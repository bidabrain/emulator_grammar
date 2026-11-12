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

PGalF reports M* in Msun/h and phi in (Mpc/h)^-3 dex^-1. We convert to PHYSICAL
units (Msun, Mpc^-3 dex^-1) using h from the catalogue header, so the GSMF lines
line up directly with observational data (Baldry+2012 / Driver+2022, physical
Msun) at the calibration step -- no further h-juggling.

By default M* is the stellar mass within a 3D 50 proper-kpc aperture (COLIBRE /
observational convention; Chaikin+2026), computed from the per-particle
GALFIND.DATA. Use --aperture-kpc to change the size, or --aperture-kpc 0 to fall
back to the catalogue's total bound stellar mass (no aperture).

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
"""

import sys
import json
import argparse
from pathlib import Path

import numpy as np
import yaml


def h_of(hubble):
    """Catalogue header stores H0 (e.g. 69.33) or h (0.6933); return h."""
    return hubble / 100.0 if (hubble and hubble > 1.0) else hubble


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
    step = f"{snapshot:05d}"
    fof_dir = Path(root) / "FoF_Data" / f"FoF.{step}"
    header, _ = read_fof(fof_dir / f"FoF_halo_cat.{step}")
    _, subs, _ = read_gal(fof_dir / f"GALCATALOG.LIST.{step}")

    h = h_of(header["hubble"])
    z = header["redshift"]

    if aperture_kpc is None or aperture_kpc <= 0:
        mstar_code = subs["mstar"]                     # Msun/h, total bound
    else:
        # proper kpc -> comoving Mpc/h (the catalogue's length units).
        radius = (aperture_kpc / 1000.0) / header["anow"] * h
        data_file = fof_dir / f"GALFIND.DATA.{step}"
        if not data_file.exists():
            raise SystemExit(f"--aperture-kpc set but {data_file} is missing")
        mstar_code = aperture_fn(data_file, radius, header["box_size"])
        if len(mstar_code) != len(subs):
            raise SystemExit(
                f"GALFIND.DATA subhalos ({len(mstar_code)}) != "
                f"GALCATALOG ({len(subs)}) in {fof_dir}")

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


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pgalf-repo", type=Path, required=True,
                   help="Path to pgalf_swift_grammar (for read_pgalf_catalog.py)")
    p.add_argument("--root-template", type=str, required=True,
                   help="Per-run PGalF output dir, with '{id}', e.g. '/data/pgalf/run_{id}'")
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
    p.add_argument("--aperture-kpc", type=float, default=50.0,
                   help="3D stellar aperture in PROPER kpc (COLIBRE/obs use 50, "
                        "read from GALFIND.DATA); pass <=0 for total bound M*")
    args = p.parse_args()

    # Import PGalF's binary readers from its repo.
    sys.path.insert(0, str(args.pgalf_repo))
    try:
        from read_pgalf_catalog import (read_fof_catalog, read_galaxy_catalog,
                                         aperture_stellar_mass)
    except ImportError as e:
        raise SystemExit(f"Could not import read_pgalf_catalog from {args.pgalf_repo}: {e}")

    with open(args.design) as fh:
        run_ids = list(json.load(fh).keys())

    args.out_dir.mkdir(parents=True, exist_ok=True)

    done, missing = [], []
    for rid in run_ids:
        root = args.root_template.format(id=rid)
        step = f"{args.snapshot:05d}"
        cat = Path(root) / "FoF_Data" / f"FoF.{step}" / f"GALCATALOG.LIST.{step}"
        if not cat.exists():
            missing.append((rid, str(cat)))
            continue

        centers, phi, scatter, counts, z = gsmf_physical(
            root, args.snapshot, read_fof_catalog, read_galaxy_catalog,
            mass_lo=args.mass_lo, mass_hi=args.mass_hi,
            dlog=args.dlog, min_count=args.min_count,
            aperture_kpc=args.aperture_kpc, aperture_fn=aperture_stellar_mass,
        )
        out = args.out_dir / f"run_{rid}.yml"
        write_emulator_yaml(out, centers, phi, scatter, counts,
                            scaling_name=args.scaling_name)
        done.append((rid, z, int(counts.sum()), len(centers)))

    ap = (f"{args.aperture_kpc:g} pkpc 3D aperture"
          if args.aperture_kpc and args.aperture_kpc > 0 else "total bound M*")
    print(f"Wrote {len(done)} GSMF files to {args.out_dir}/  "
          f"(scaling='{args.scaling_name}', M*={ap})")
    for rid, z, ntot, nbin in done:
        print(f"  run {rid:>3}: z={z:.3f}  {ntot:>7d} galaxies  {nbin} bins")
    if missing:
        print(f"\nMISSING catalogues for {len(missing)} runs (not yet processed by PGalF):")
        for rid, path in missing:
            print(f"  run {rid:>3}: {path}")
    if done:
        # NOTE: use the written 'scatter' field + sanitize_output=True. Do NOT pass
        # histogram_for_error together with log_dependent -- swiftemulator leaves the
        # values linear in that path (swift.py:155 is skipped), giving log errors on
        # linear values. The 'histogram' line is kept only for reference.
        print("\nNext: train the emulator, e.g.")
        print("  from swiftemulator.io.swift import load_pipeline_outputs")
        print(f"  values, units = load_pipeline_outputs(")
        print(f"      filenames={{rid: '{args.out_dir}/run_'+str(rid)+'.yml' for rid in ...}},")
        print(f"      scaling_relations=['{args.scaling_name}'],")
        print(f"      log_independent=['{args.scaling_name}'], log_dependent=['{args.scaling_name}'],")
        print(f"      sanitize_output=True)")


if __name__ == "__main__":
    main()
