"""
Stage 2 helper: turn a Latin-hypercube design into ready-to-submit SWIFT runs.

Reads
  --design    the directory produced by design_gsmf.py (run_<id>.yml files,
              plus design_parameters.json)
  --template  a "working-directory template": every file a SWIFT run needs,
              i.e. symlinks to the ICs / cooling / yield / photometry / X-ray
              tables, output_list.txt, the swift executable, etc.
  --sbatch    a SLURM job-script template with {{...}} placeholders (below)

and, for every run in the design, creates   <out>/run_<id>/   containing
  * a copy of the whole template directory (symlinks kept AS symlinks, so the
    big ICs/tables are not duplicated 48x),
  * the run's SWIFT parameter file (design_gsmf/run_<id>.yml -> <param-name>),
  * a rendered SLURM script (the --sbatch template with placeholders filled).

So after running this you can submit each simulation with its own sbatch.

Placeholders understood in the --sbatch template (use {{ }}, which never
clashes with bash ${VAR} / ${arr[@]}):
    {{RUN_ID}}      run number, e.g. 0
    {{RUN_NAME}}    run_<id>, e.g. run_0
    {{RUN_DIR}}     absolute path of the run directory
    {{PARAM_FILE}}  parameter-file name inside the run dir (--param-name)

    python setup_runs.py \
        --design   design_gsmf \
        --template run_template \
        --sbatch   mpirun.sbatch \
        --out      /scratch/gsmf_runs \
        [--param-name params.yml] [--sbatch-name job.sbatch] \
        [--submit] [--dry-run] [--force]

Without --submit it only prepares the directories and prints the submit
commands. With --submit it runs `sbatch job.sbatch` from inside each run dir.
"""

import re
import sys
import json
import shutil
import argparse
import subprocess
from pathlib import Path


def discover_runs(design_dir):
    """Return [(run_id:int, param_path:Path), ...] sorted by id, from run_*.yml."""
    runs = []
    for p in design_dir.glob("run_*.yml"):
        m = re.fullmatch(r"run_(\d+)\.yml", p.name)
        if m:
            runs.append((int(m.group(1)), p))
    runs.sort(key=lambda t: t[0])
    return runs


def render_sbatch(template_text, *, run_id, run_name, run_dir, param_file):
    """Fill the {{...}} placeholders. Returns (text, used_any:bool)."""
    subs = {
        "{{RUN_ID}}": str(run_id),
        "{{RUN_NAME}}": run_name,
        "{{RUN_DIR}}": str(run_dir),
        "{{PARAM_FILE}}": param_file,
    }
    used = any(k in template_text for k in subs)
    text = template_text
    for key, val in subs.items():
        text = text.replace(key, val)
    return text, used


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--design", type=Path, required=True,
                   help="design_gsmf.py output dir (run_<id>.yml + design_parameters.json)")
    p.add_argument("--template", type=Path, required=True,
                   help="working-dir template copied into every run folder")
    p.add_argument("--sbatch", type=Path, required=True,
                   help="SLURM template with {{RUN_ID}}/{{RUN_NAME}}/{{RUN_DIR}}/{{PARAM_FILE}}")
    p.add_argument("--out", type=Path, required=True,
                   help="base dir under which run_<id>/ folders are created")
    p.add_argument("--param-name", default="flamingo.yml",
                   help="name of the SWIFT parameter file inside each run dir "
                        "(default: flamingo.yml -- must match the name your sbatch passes to swift)")
    p.add_argument("--sbatch-name", default=None,
                   help="name of the rendered SLURM script inside each run dir "
                        "(default: same basename as --sbatch)")
    p.add_argument("--only", type=str, default=None,
                   help="comma-separated run ids to process (default: all)")
    p.add_argument("--force", action="store_true",
                   help="overwrite an existing run dir's template files / scripts")
    p.add_argument("--submit", action="store_true",
                   help="run `sbatch <sbatch-name>` from inside each run dir")
    p.add_argument("--dry-run", action="store_true",
                   help="show what would happen, create nothing")
    args = p.parse_args()

    if args.sbatch_name is None:
        args.sbatch_name = args.sbatch.name

    for label, path in (("design", args.design), ("template", args.template),
                        ("sbatch", args.sbatch)):
        if not path.exists():
            raise SystemExit(f"--{label} not found: {path}")

    runs = discover_runs(args.design)
    if not runs:
        raise SystemExit(f"No run_*.yml found in {args.design}")

    if args.only is not None:
        wanted = {int(x) for x in args.only.split(",") if x.strip() != ""}
        runs = [(i, p) for (i, p) in runs if i in wanted]
        if not runs:
            raise SystemExit(f"--only {args.only} matched none of the design runs")

    # Optional: sanity-log parameter values alongside each run.
    design_params = {}
    dp_file = args.design / "design_parameters.json"
    if dp_file.is_file():
        with open(dp_file) as fh:
            design_params = json.load(fh)

    template_text = args.sbatch.read_text()
    probe, used = render_sbatch(template_text, run_id=0, run_name="run_0",
                                run_dir=Path("/x"), param_file=args.param_name)
    if not used:
        print(f"note: {args.sbatch} has no {{{{...}}}} placeholders -> every run gets an "
              "identical sbatch. That is fine if it uses RELATIVE paths and is submitted "
              "from inside each run dir (e.g. `./swift_mpi ... flamingo.yml`).",
              file=sys.stderr)

    args.out.mkdir(parents=True, exist_ok=True)
    submitted = []
    for run_id, param_src in runs:
        run_name = f"run_{run_id}"
        run_dir = (args.out / run_name).resolve()
        print(f"[{run_name}] -> {run_dir}")

        if args.dry_run:
            if str(run_id) in design_params:
                print(f"    params: {design_params[str(run_id)]}")
            continue

        # 1) copy the template working dir (keep symlinks as symlinks).
        try:
            shutil.copytree(args.template, run_dir, symlinks=True, dirs_exist_ok=args.force)
        except FileExistsError:
            raise SystemExit(f"{run_dir} already exists; pass --force to overwrite, "
                             "or --out a fresh directory")

        # 2) drop in this run's SWIFT parameter file.
        shutil.copyfile(param_src, run_dir / args.param_name)

        # 3) render + write the per-run SLURM script.
        text, _ = render_sbatch(template_text, run_id=run_id, run_name=run_name,
                                run_dir=run_dir, param_file=args.param_name)
        sbatch_path = run_dir / args.sbatch_name
        sbatch_path.write_text(text)

        # 4) optionally submit.
        if args.submit:
            r = subprocess.run(["sbatch", args.sbatch_name], cwd=run_dir,
                               capture_output=True, text=True)
            if r.returncode == 0:
                print(f"    submitted: {r.stdout.strip()}")
                submitted.append(run_name)
            else:
                print(f"    sbatch FAILED: {r.stderr.strip()}", file=sys.stderr)

    if args.dry_run:
        print(f"\n[dry-run] {len(runs)} runs would be prepared under {args.out}/")
        return

    print(f"\nPrepared {len(runs)} run dirs under {args.out}/")
    if args.submit:
        print(f"Submitted {len(submitted)}/{len(runs)} jobs.")
    else:
        print("To submit, either re-run with --submit, or:")
        print(f"  for d in {args.out}/run_*; do (cd $d && sbatch {args.sbatch_name}); done")


if __name__ == "__main__":
    main()
