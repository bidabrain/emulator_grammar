"""
Bridge PGalF galaxy catalogues -> SWIFT-Emulator GSMF input files.

For every run in a Latin-hypercube design (design_gsmf/design_parameters.json),
read that run's PGalF catalogue at the z=0 snapshot, bin the galaxy stellar
masses into a GSMF, and write one YAML per run in the exact format that
swiftemulator.io.swift.load_pipeline_outputs expects:

    stellar_mass_function:
      lines:
        mass_function:
          centers:       [...]          # M* bin centres  (LINEAR Msun)
          centers_units: Solar_Mass
          values:        [...]          # phi             (LINEAR Mpc^-3 dex^-1)
          values_units:  "1/Mpc**3"
          scatter:       [...]          # Poisson error on phi (LINEAR) -> use this
        histogram:
          values:        [...]          # raw galaxy counts per bin (reference only)

PGalF reports M* in Msun/h and phi in (Mpc/h)^-3 dex^-1 (only if NewDD found h;
see check_little_h). We convert to PHYSICAL
units (Msun, Mpc^-3 dex^-1) using h from the catalogue header, so the GSMF lines
line up directly with observational data (Baldry+2012 / Driver+2022, physical
Msun) at the calibration step -- no further h-juggling.

By default M* is the stellar mass within a 3D 50 proper-kpc aperture (COLIBRE /
observational convention; Chaikin+2026), computed from the per-particle
GALFIND.DATA. Use --aperture-kpc to change the size, or --aperture-kpc 0 to fall
back to the catalogue's total bound stellar mass (no aperture). Several values
(e.g. --aperture-kpc 0 30 50 100) write one set of YAMLs per aperture into
OUT_DIR/ap_<X>kpc/ (ap_total/ for 0), reading GALFIND.DATA once per run.

Run IDs are taken from design_parameters.json so they MATCH the emulator's
model_parameters keys. Each run is a separate SWIFT simulation, so each has its
own PGalF output directory; give the directory as a template with "{id}".

    python measure_gsmf.py \
        --pgalf-repo   /Volumes/Elements/scratch/pgalf_swift_grammar \
        --root-template "/path/to/pgalf/run_{id}" \
        --snapshot     42 \
        --design       design_gsmf/design_parameters.json \
        --out-dir      gsmf_measured

Then train the emulator on {design_parameters.json, gsmf_measured/run_*.yml}.

With --finder hbt the catalogues come from HBT-HERONS instead (run with
HBT-HERONS-GRAMMAR/grammar/run_hbt.sh; --hbt-repo points at that repo and
--root-template at the per-run HBT output dirs). Same output format.
"""

import sys
import json
import inspect
import argparse
from pathlib import Path

import numpy as np
import yaml


def h_of(hubble):
    """Catalogue header stores H0 (e.g. 69.33) or h (0.6933); return h."""
    return hubble / 100.0 if (hubble and hubble > 1.0) else hubble


_MPC_CM = 3.08567758e24


def check_little_h(root, snapshot, header):
    """Check that the PGalF catalogue really is in little-h units.

    PGalF's NewDD used to miss /Cosmology/h on SWIFT input and silently use
    h=1, leaving box, positions and masses in cMpc / Msun (see pgalf
    NewDD/NOTE_SWIFT_h_units.md). The FoF header's `hubble` is overwritten with
    a compile-time H0, so it cannot reveal this; compare the FoF box_size with
    the SWIFT snapshot that run_pgalf.sh links as <root>/snap_XXXX.hdf5.
    Returns ("ok" | "pre-fix" | "mismatch" | "unchecked", message)."""
    snap = Path(root) / f"snap_{snapshot:04d}.hdf5"
    if not snap.exists():
        return "unchecked", f"{snap} not found"
    try:
        import h5py
    except ImportError:
        return "unchecked", "h5py not available"
    with h5py.File(snap, "r") as f:
        box = float(np.atleast_1d(f["Header"].attrs["BoxSize"])[0])
        if "HubbleParam" in f["Header"].attrs:
            h = float(np.atleast_1d(f["Header"].attrs["HubbleParam"])[0])
        elif "Cosmology" in f and "h" in f["Cosmology"].attrs:
            h = float(np.atleast_1d(f["Cosmology"].attrs["h"])[0])
        else:
            return "unchecked", f"no h in {snap}"
        unit_l = _MPC_CM
        if "InternalCodeUnits" in f:
            unit_l = float(np.atleast_1d(
                f["InternalCodeUnits"].attrs["Unit length in cgs (U_L)"])[0])
    box_mpc = box * unit_l / _MPC_CM
    got = float(header["box_size"])
    if abs(got / (box_mpc * h) - 1) < 1e-3:
        return "ok", f"box_size={got:.4g} cMpc/h = {box_mpc:g} Mpc x h={h:g}"
    if abs(got / box_mpc - 1) < 1e-3:
        return "pre-fix", (f"box_size={got:.4g} equals the snapshot BoxSize in Mpc "
                           f"({box_mpc:g}), not x h={h:g}: catalogue made by the old "
                           "NewDD (h=1, units cMpc/Msun). Re-run PGalF from NewDD "
                           "(pgalf_swift_grammar >= 2570de6).")
    return "mismatch", (f"box_size={got:.4g} matches neither {box_mpc:g} Mpc nor "
                        f"{box_mpc * h:.4g} cMpc/h from {snap}")


def _is_aperture(a):
    return a is not None and a > 0


def _aperture_tag(a):
    """Sub-directory name for one aperture in multi-aperture mode."""
    return f"ap_{a:g}kpc" if _is_aperture(a) else "ap_total"


def stellar_masses_code(root, snapshot, read_fof, read_gal, apertures_kpc,
                        aperture_fn=None, centers_fn=None):
    """Per-subhalo M* [Msun/h] for each entry of apertures_kpc, plus the FoF
    header. An entry >0 is a 3D PROPER-kpc aperture (COLIBRE/obs convention)
    read from GALFIND.DATA; <=0 or None is the catalogue's total bound M*.
    GALFIND.DATA is read once for all apertures (aperture_fn must be
    read_pgalf_catalog.aperture_stellar_mass, which accepts a list of radii).
    Apertures are centred on each subhalo's centre of mass unless centers_fn
    (read_pgalf_catalog.galcenter_positions) is given, in which case the
    GalCenter stellar centres from GALFIND.CENTER are used (COM where missing)."""
    step = f"{snapshot:05d}"
    fof_dir = Path(root) / "FoF_Data" / f"FoF.{step}"
    header, _ = read_fof(fof_dir / f"FoF_halo_cat.{step}")
    _, subs, _ = read_gal(fof_dir / f"GALCATALOG.LIST.{step}")

    ap_kpc = [a for a in apertures_kpc if _is_aperture(a)]
    by_ap = {}
    if ap_kpc:
        h = h_of(header["hubble"])
        # proper kpc -> comoving Mpc/h (the catalogue's length units).
        radii = [(a / 1000.0) / header["anow"] * h for a in ap_kpc]
        data_file = fof_dir / f"GALFIND.DATA.{step}"
        if not data_file.exists():
            raise SystemExit(f"--aperture-kpc set but {data_file} is missing")
        kw = {}
        if centers_fn is not None:
            center_file = fof_dir / f"GALFIND.CENTER.{step}"
            if not center_file.exists():
                raise SystemExit(f"--center galcenter but {center_file} is missing")
            kw["centers"] = centers_fn(center_file, fof_dir / f"GALCATALOG.LIST.{step}")
            if len(kw["centers"]) != len(subs):
                raise SystemExit(f"GALFIND.CENTER aligned to {len(kw['centers'])} "
                                 f"subhalos != GALCATALOG ({len(subs)}) in {fof_dir}")
        if len(radii) == 1:
            mstar_ap = aperture_fn(data_file, radii[0], header["box_size"], **kw)[:, None]
        else:
            try:
                mstar_ap = aperture_fn(data_file, radii, header["box_size"], **kw)
            except TypeError:
                raise SystemExit("read_pgalf_catalog.aperture_stellar_mass takes "
                                 "a single radius only; git pull pgalf_swift_grammar "
                                 "(needs >= ba53e2d) for several --aperture-kpc")
        if len(mstar_ap) != len(subs):
            raise SystemExit(
                f"GALFIND.DATA subhalos ({len(mstar_ap)}) != "
                f"GALCATALOG ({len(subs)}) in {fof_dir}")
        by_ap = {a: mstar_ap[:, j] for j, a in enumerate(ap_kpc)}
    masses = [by_ap[a] if _is_aperture(a) else subs["mstar"] for a in apertures_kpc]
    return masses, header


def gsmf_physical(root, snapshot, read_fof, read_gal,
                  mass_lo=7.0, mass_hi=12.0, dlog=0.25, min_count=1,
                  aperture_kpc=None, aperture_fn=None):
    """Return (centers_Msun, phi, scatter, counts, z) in PHYSICAL units.

    centers_Msun : linear M* bin centres [Msun]
    phi          : linear number density [Mpc^-3 dex^-1]
    scatter      : linear Poisson error on phi [Mpc^-3 dex^-1]
    counts       : raw galaxy count per kept bin

    aperture_kpc : if >0, use the stellar mass within this 3D aperture (PROPER
                   kpc, COLIBRE/obs convention) read from GALFIND.DATA instead of
                   the catalogue's total bound M*. aperture_fn must then be
                   read_pgalf_catalog.aperture_stellar_mass.
    """
    (mstar_code,), header = stellar_masses_code(
        root, snapshot, read_fof, read_gal, [aperture_kpc], aperture_fn)
    return gsmf_from_masses(mstar_code, header, mass_lo, mass_hi, dlog, min_count)


def gsmf_from_masses(mstar_code, header, mass_lo=7.0, mass_hi=12.0, dlog=0.25,
                     min_count=1):
    """Bin catalogue M* [Msun/h] into a PHYSICAL-unit GSMF; same return as
    gsmf_physical."""
    h = h_of(header["hubble"])
    z = header["redshift"]
    mstar = mstar_code[mstar_code > 0] / h             # Msun/h -> Msun
    volume = (header["box_size"] / h) ** 3             # (Mpc/h)^3 -> Mpc^3

    edges = np.arange(mass_lo, mass_hi + 1e-9, dlog)
    counts, _ = np.histogram(np.log10(mstar), bins=edges)
    log_centers = 0.5 * (edges[:-1] + edges[1:])

    phi = counts / (volume * dlog)                     # Mpc^-3 dex^-1
    scatter = np.where(counts > 0, phi / np.sqrt(np.maximum(counts, 1)), 0.0)

    good = counts >= min_count
    return (
        10.0 ** log_centers[good],                     # linear Msun
        phi[good],
        scatter[good],
        counts[good],
        z,
    )


def write_emulator_yaml(path, centers, phi, scatter, counts,
                        scaling_name="stellar_mass_function"):
    """Write one GSMF file in load_pipeline_outputs format."""
    doc = {
        scaling_name: {
            "lines": {
                "mass_function": {
                    "centers": [float(x) for x in centers],
                    "centers_units": "Solar_Mass",
                    "values": [float(x) for x in phi],
                    "values_units": "1/Mpc**3",
                    "scatter": [float(x) for x in scatter],
                },
                # Raw galaxy counts per bin, kept for reference / diagnostics.
                "histogram": {
                    "values": [int(c) for c in counts],
                },
            }
        }
    }
    with open(path, "w") as fh:
        yaml.safe_dump(doc, fh, default_flow_style=False, sort_keys=False)


def write_run(args, rid, masses, header, apertures, out_dirs, done):
    """Bin one run's catalogue M* [Msun/h] per aperture and write its YAMLs."""
    for a, mstar_code, d in zip(apertures, masses, out_dirs):
        centers, phi, scatter, counts, z = gsmf_from_masses(
            mstar_code, header, mass_lo=args.mass_lo, mass_hi=args.mass_hi,
            dlog=args.dlog, min_count=args.min_count)
        write_emulator_yaml(d / f"run_{rid}.yml", centers, phi, scatter, counts,
                            scaling_name=args.scaling_name)
        done[a].append((rid, z, int(counts.sum()), len(centers)))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--finder", choices=("pgalf", "hbt"), default="pgalf",
                   help="galaxy finder whose catalogues to read: PGalF (default) or "
                        "HBT-HERONS (HBT-HERONS-GRAMMAR/grammar/run_hbt.sh outputs)")
    p.add_argument("--pgalf-repo", type=Path, default=None,
                   help="Path to pgalf_swift_grammar (for read_pgalf_catalog.py); "
                        "required for --finder pgalf")
    p.add_argument("--hbt-repo", type=Path, default=None,
                   help="Path to HBT-HERONS-GRAMMAR (for grammar/read_hbt_catalog.py); "
                        "required for --finder hbt")
    p.add_argument("--root-template", type=str, required=True,
                   help="Per-run finder output dir, with '{id}', e.g. '/data/pgalf/run_{id}' "
                        "(PGalF run_pgalf.sh -o) or '/data/hbt/run_{id}' (run_hbt.sh -o)")
    p.add_argument("--snapshot", type=int, required=True,
                   help="z=0 snapshot number (same for every run)")
    p.add_argument("--design", type=Path, default=Path("design_gsmf/design_parameters.json"),
                   help="design_parameters.json (defines the run IDs)")
    p.add_argument("--out-dir", type=Path, default=Path("gsmf_measured"))
    p.add_argument("--scaling-name", type=str, default="stellar_mass_function")
    p.add_argument("--mass-lo", type=float, default=7.0)
    p.add_argument("--mass-hi", type=float, default=12.0)
    p.add_argument("--dlog", type=float, default=0.25, help="mass bin width [dex]")
    p.add_argument("--min-count", type=int, default=1,
                   help="drop bins with fewer than this many galaxies")
    p.add_argument("--aperture-kpc", type=float, nargs="+", default=[50.0],
                   help="3D stellar aperture(s) in PROPER kpc (COLIBRE/obs use 50, "
                        "read from GALFIND.DATA); <=0 means total bound M*. "
                        "Several values (e.g. 0 30 50 100) write one set of YAMLs "
                        "per aperture into OUT_DIR/ap_<X>kpc/ (or ap_total/), "
                        "reading GALFIND.DATA once per run")
    p.add_argument("--center", choices=("com", "galcenter", "mostbound"), default=None,
                   help="aperture centre. PGalF: 'com' (default) = subhalo centre of "
                        "mass of all particle types, 'galcenter' = GalCenter stellar "
                        "centre (GALFIND.CENTER), COM where missing. HBT: 'mostbound' "
                        "(default) = most bound particle, 'com' = centre of mass of the "
                        "bound particles")
    p.add_argument("--skip-h-check", action="store_true",
                   help="do not compare the FoF box_size with <root>/snap_XXXX.hdf5 "
                        "(the check catches catalogues made by the old NewDD, which "
                        "used h=1 on SWIFT input)")
    args = p.parse_args()

    if args.center is None:
        args.center = "com" if args.finder == "pgalf" else "mostbound"
    if args.finder == "hbt":
        if args.center == "galcenter":
            p.error("--center galcenter is PGalF-only; use mostbound or com with --finder hbt")
        if args.hbt_repo is None:
            p.error("--finder hbt needs --hbt-repo")
        sys.path.insert(0, str(args.hbt_repo / "grammar"))
        try:
            from read_hbt_catalog import load_stellar_masses as hbt_masses, find_output
        except ImportError as e:
            raise SystemExit(f"Could not import grammar/read_hbt_catalog from {args.hbt_repo}: {e}")
        measure_runs(args, hbt_masses=hbt_masses, hbt_find=find_output)
        return
    if args.center == "mostbound":
        p.error("--center mostbound is HBT-only; use com or galcenter with --finder pgalf")
    if args.pgalf_repo is None:
        p.error("--finder pgalf needs --pgalf-repo")

    # Import PGalF's binary readers from its repo.
    sys.path.insert(0, str(args.pgalf_repo))
    try:
        from read_pgalf_catalog import (read_fof_catalog, read_galaxy_catalog,
                                         aperture_stellar_mass)
    except ImportError as e:
        raise SystemExit(f"Could not import read_pgalf_catalog from {args.pgalf_repo}: {e}")
    # Old readers hard-code the GALFIND.DATA record sizes (wrong for NCHEM/NDUST
    # builds -> "Truncated subhalo header"). Several apertures need >= ba53e2d.
    if any(_is_aperture(a) for a in args.aperture_kpc) and \
            "sizes" not in inspect.signature(aperture_stellar_mass).parameters:
        raise SystemExit(f"read_pgalf_catalog in {args.pgalf_repo} is too old "
                         "(hard-coded GALFIND.DATA record sizes); git pull it "
                         "(needs >= 3d0b545, or >= ba53e2d for several apertures)")
    centers_fn = None
    if args.center == "galcenter":
        try:
            from read_pgalf_catalog import galcenter_positions as centers_fn
        except ImportError:
            raise SystemExit(f"read_pgalf_catalog in {args.pgalf_repo} has no "
                             "galcenter_positions; git pull it for --center galcenter")
    measure_runs(args, read_fof_catalog=read_fof_catalog,
                 read_galaxy_catalog=read_galaxy_catalog,
                 aperture_stellar_mass=aperture_stellar_mass, centers_fn=centers_fn)


def measure_runs(args, read_fof_catalog=None, read_galaxy_catalog=None,
                 aperture_stellar_mass=None, centers_fn=None, hbt_masses=None,
                 hbt_find=None):
    """Measure and write the GSMF of every design run with the finder chosen in
    args (PGalF readers or HBT-HERONS reader functions passed in)."""
    with open(args.design) as fh:
        run_ids = list(json.load(fh).keys())

    apertures = args.aperture_kpc
    multi = len(apertures) > 1
    out_dirs = [args.out_dir / _aperture_tag(a) if multi else args.out_dir
                for a in apertures]
    for d in out_dirs:
        d.mkdir(parents=True, exist_ok=True)

    done, missing, unchecked = {a: [] for a in apertures}, [], []
    for rid in run_ids:
        root = args.root_template.format(id=rid)
        if args.finder == "hbt":
            # HBT-HERONS reads units, h and box from the SWIFT metadata, so the
            # PGalF little-h check does not apply.
            try:
                hbt_find(root, args.snapshot)
            except FileNotFoundError as e:
                missing.append((rid, str(e)))
                continue
            masses, header = hbt_masses(root, args.snapshot, apertures, center=args.center)
            write_run(args, rid, masses, header, apertures, out_dirs, done)
            continue
        step = f"{args.snapshot:05d}"
        cat = Path(root) / "FoF_Data" / f"FoF.{step}" / f"GALCATALOG.LIST.{step}"
        if not cat.exists():
            missing.append((rid, str(cat)))
            continue

        if not args.skip_h_check:
            fof = Path(root) / "FoF_Data" / f"FoF.{step}" / f"FoF_halo_cat.{step}"
            status, msg = check_little_h(root, args.snapshot, read_fof_catalog(fof)[0])
            if status in ("pre-fix", "mismatch"):
                raise SystemExit(f"run {rid}: {msg}\n(use --skip-h-check to override)")
            if status == "unchecked":
                unchecked.append((rid, msg))

        masses, header = stellar_masses_code(
            root, args.snapshot, read_fof_catalog, read_galaxy_catalog,
            apertures, aperture_fn=aperture_stellar_mass, centers_fn=centers_fn)
        write_run(args, rid, masses, header, apertures, out_dirs, done)

    for a, d in zip(apertures, out_dirs):
        ap = (f"{a:g} pkpc 3D aperture, centre={args.center}" if _is_aperture(a)
              else "total bound M*")
        print(f"Wrote {len(done[a])} GSMF files to {d}/  "
              f"(finder={args.finder}, scaling='{args.scaling_name}', M*={ap})")
        for rid, z, ntot, nbin in done[a]:
            print(f"  run {rid:>3}: z={z:.3f}  {ntot:>7d} galaxies  {nbin} bins")
    if unchecked:
        print(f"\nlittle-h check SKIPPED for {len(unchecked)} runs "
              "(cannot confirm the catalogue is in h units):")
        for rid, msg in unchecked:
            print(f"  run {rid:>3}: {msg}")
    if missing:
        print(f"\nMISSING catalogues for {len(missing)} runs (not yet processed by {args.finder}):")
        for rid, path in missing:
            print(f"  run {rid:>3}: {path}")
    if any(done.values()):
        # NOTE: use the written 'scatter' field + sanitize_output=True. Do NOT pass
        # histogram_for_error together with log_dependent -- swiftemulator leaves the
        # values linear in that path (swift.py:155 is skipped), giving log errors on
        # linear values. The 'histogram' line is kept only for reference.
        print("\nNext: train the emulator, e.g.")
        print("  from swiftemulator.io.swift import load_pipeline_outputs")
        print(f"  values, units = load_pipeline_outputs(")
        print(f"      filenames={{rid: '{out_dirs[0]}/run_'+str(rid)+'.yml' for rid in ...}},")
        print(f"      scaling_relations=['{args.scaling_name}'],")
        print(f"      log_independent=['{args.scaling_name}'], log_dependent=['{args.scaling_name}'],")
        print(f"      sanitize_output=True)")


if __name__ == "__main__":
    main()
