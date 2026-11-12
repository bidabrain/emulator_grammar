"""Standalone cosmic star-formation-history plot for one or more SWIFT runs.

Reads each run directory's SFR.txt (SWIFT's star-formation log) and overplots
the cosmic SFR density vs redshift, compared to three analytic fits (Madau &
Dickinson 2014, Behroozi+2013, Harikane+2022; same as the swift-pipeline
general SFH script). Independent of the swift-pipeline framework; it only needs
`swiftsimio` (to read units + box size from one snapshot per run) and `unyt`.

The SFR curve comes entirely from SFR.txt (columns a, z, SFR in internal units,
as in the swift-pipeline script); one snapshot per run is loaded ONLY for its
units and box size, so any available snapshot in the directory works.

    python plot_sfh.py /path/run_0 /path/run_1 \
        --labels run0 run1 --output sfh.png

    # all design runs at once:
    python plot_sfh.py /scratch/gsmf_runs/run_* --output sfh_all.png
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import unyt

from swiftsimio import load

SFR_OUTPUT_UNITS = unyt.msun / (unyt.year * unyt.Mpc ** 3)

# Columns in SWIFT's SFR.txt (0-indexed), matching the swift-pipeline script.
_COL_A, _COL_Z, _COL_SFR = 2, 3, 7


# Analytic cosmic SFRD fits [Msun/yr/Mpc^3], h=0.7, matching pipeline-configs
# general/star_formation_history.py. The SWIFT models use a Chabrier IMF, so
# MD14 and Harikane+22 (native Salpeter) are divided by 1.7; Behroozi+13 is
# already Chabrier. No h-rescaling (fits assume h=0.7, ~few-percent difference).
def madau_dickinson_2014(z):
    """Madau & Dickinson (2014), ARA&A 52, Eq. 15 (Salpeter->Chabrier /1.7)."""
    z = np.asarray(z, dtype=float)
    return (0.015 * (1.0 + z) ** 2.7 / (1.0 + ((1.0 + z) / 2.9) ** 5.6)) / 1.7


def behroozi_2013(z):
    """Behroozi, Wechsler & Conroy (2013), ApJ 770, 57 (native Chabrier)."""
    z = np.asarray(z, dtype=float)
    return 0.180 / (10.0 ** (-0.997 * (z - 1.243))
                    + 10.0 ** (0.241 * (z - 1.243)))


def harikane_2022(z):
    """Harikane et al. (2022), ApJS 259, 20, Eq. 60 (Salpeter->Chabrier /1.7)."""
    z = np.asarray(z, dtype=float)
    return (1.0 / (61.7 * (1.0 + z) ** -3.13
                   + 10.0 ** (0.22 * (1.0 + z))
                   + 2.4 * 10.0 ** (0.50 * (1.0 + z) - 3.0))) / 1.7


_SFRD_FITS = [
    ("Madau & Dickinson (2014)", madau_dickinson_2014, "dashed"),
    ("Behroozi et al. (2013)", behroozi_2013, "dotted"),
    ("Harikane et al. (2022)", harikane_2022, "dashdot"),
]


def find_snapshot(run_dir, name=None):
    """A snapshot .hdf5 in run_dir for metadata (any one works). Skips ICs/FoF."""
    d = Path(run_dir)
    if name:
        p = d / name
        if not p.exists():
            raise SystemExit(f"snapshot {p} not found")
        return p
    snaps = [c for c in sorted(d.glob("*.hdf5"))
             if "ics" not in c.name.lower() and "fof" not in c.name.lower()]
    if not snaps:
        raise SystemExit(f"no snapshot .hdf5 found in {d} (needed for units/box)")
    return snaps[-1]


def run_sfh(run_dir, snapshot_name=None):
    """Return (scale_factor, redshift, sfr_density[Msun/yr/Mpc^3]) for a run."""
    d = Path(run_dir)
    sfr_file = d / "SFR.txt"
    if not sfr_file.exists():
        raise SystemExit(f"{sfr_file} missing")

    data = np.genfromtxt(sfr_file).T
    if data.ndim != 2 or data.shape[0] <= _COL_SFR:
        raise SystemExit(f"{sfr_file}: expected >= {_COL_SFR + 1} columns")

    a = data[_COL_A]
    z = data[_COL_Z]
    sfr_internal = data[_COL_SFR]

    snap = load(str(find_snapshot(d, snapshot_name)))
    box = snap.metadata.boxsize
    volume = box[0] * box[1] * box[2]
    sfr_units = snap.gas.star_formation_rates.units
    sfr_density = (sfr_internal * sfr_units / volume).to(SFR_OUTPUT_UNITS)

    return a, z, sfr_density.value


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("run_dirs", nargs="+", type=Path,
                   help="one or more SWIFT run directories (each with SFR.txt)")
    p.add_argument("--labels", nargs="+", default=None,
                   help="legend labels (default: directory names)")
    p.add_argument("--snapshot-name", default=None,
                   help="snapshot filename inside each run dir (default: auto-detect)")
    p.add_argument("--no-obs", action="store_true", help="skip the observational fit")
    p.add_argument("--output", type=Path, default=Path("star_formation_history.png"))
    args = p.parse_args()

    if args.labels and len(args.labels) != len(args.run_dirs):
        raise SystemExit("--labels count must match the number of run dirs")

    fig, ax = plt.subplots(figsize=(7, 5.5))
    ax.loglog()

    for i, d in enumerate(args.run_dirs):
        a, z, sfr = run_sfh(d, args.snapshot_name)
        label = args.labels[i] if args.labels else Path(d).name
        ax.plot(a, sfr, lw=1.3, zorder=100, label=label)

    if not args.no_obs:
        z_obs = np.logspace(-3, np.log10(20.0), 300)
        a_obs = 1.0 / (1.0 + z_obs)
        for label, fit, ls in _SFRD_FITS:
            ax.plot(a_obs, fit(z_obs), color="black", ls=ls, lw=1.1,
                    zorder=1, label=label)

    # Redshift ticks on the (scale-factor) x-axis.
    z_ticks = np.array([0.0, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0])
    ax.set_xticks(1.0 / (1.0 + z_ticks))
    ax.set_xticklabels([f"${z:g}$" for z in z_ticks])
    ax.set_xlim(1.02, 0.055)
    ax.set_ylim(1.8e-4, 1.7)
    ax.tick_params(axis="x", which="minor", bottom=False)

    ax.set_xlabel("Redshift $z$")
    ax.set_ylabel(r"SFR density $\dot{\rho}_\star$ [M$_\odot$ yr$^{-1}$ Mpc$^{-3}$]")
    ax.legend(frameon=False, fontsize=8, loc="lower center", ncol=2)
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=160)
    plt.close(fig)
    print(f"wrote {args.output.resolve()}  ({len(args.run_dirs)} runs)")


if __name__ == "__main__":
    main()
