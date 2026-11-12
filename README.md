SWIFT-Emulator
==============

[![Documentation Status](https://readthedocs.org/projects/swiftemulator/badge/?version=latest)](https://swiftemulator.readthedocs.io/en/latest/?badge=latest)
![Test Status](https://github.com/swiftsim/emulator/actions/workflows/pytest.yml/badge.svg)
[![PyPI version](https://badge.fury.io/py/swiftemulator.svg)](https://badge.fury.io/py/swiftemulator)
[![status](https://joss.theoj.org/papers/61d082196ef861cc0b612486c1fa6d40/status.svg)](https://joss.theoj.org/papers/61d082196ef861cc0b612486c1fa6d40)

The SWIFT-emulator (henceforth 'the emulator') was initially designed for [SWIFT](http://swift.dur.ac.uk)
outputs, and includes utilities to read and write SWIFT data.

The emulator can be used used to predict 
outputs of simulations without having to run them, by employing Gaussian Process
Regression with `george` and sensitivity analysis with `SALib`.

Dcumentation is available at [ReadTheDocs](https://swiftemulator.readthedocs.io/).

GSMF Calibration Pipeline (this fork)
-------------------------------------

This fork adds a pipeline that calibrates the 6 SWIFT feedback parameters
(dual-channel thermal+kinetic SN + spin-jet AGN, `ColibrePressure` f_E model)
against the **z=0 galaxy stellar mass function (GSMF)** of a 25 Mpc box, using
the same emulator method as COLIBRE (Chaikin et al. 2026) / FLAMINGO
(Kugel et al. 2023). Full details in [`GSMF_CALIBRATION.md`](GSMF_CALIBRATION.md).

```
design_gsmf.py  ->  setup_runs.py  ->  SWIFT  ->  pgalf  ->  measure_gsmf.py  ->  train + calibrate
 Latin hypercube    per-run dirs       48 sims   galaxy      catalogue -> GSMF    GP emulator ->
 48 x run_X.yml     + sbatch                      catalogue   emulator YAML        fit to obs GSMF
```

The 6 calibrated knobs (see `design_gsmf.py`): `SNII_f_kinetic`, `SNII_delta_T_K`,
`SNII_energy_fraction_min/max`, `SNII_energy_fraction_P_pivot_K_p_cm3`,
`SPINJETAGN:subgrid_seed_mass_Msun`.

**Stage 0 — environment** (once)

Run the scripts from a Python env with this repo (editable) installed. The
scripts import from source, so install the package and its dependencies:

```bash
# optional but recommended: a dedicated env instead of conda (base)
conda create -n emulator python=3.11 -y && conda activate emulator

cd /path/to/emulator_grammar
pip install -e .                 # swiftemulator + deps (attrs, george, SALib, corner, ...)
pip install "pyDOE==0.3.8"       # Latin hypercube; NOT plain `pip install pyDOE` (wrong 'pydoe' 1.5.0)
pip install emcee                # only needed for the Stage-5 MCMC calibration
```

A bare `ModuleNotFoundError: No module named 'attr'` means the deps aren't
installed yet — run `pip install -e .` above (`attr` is the import name of `attrs`).

**Stage 1 — design the parameter sweep** (`design_gsmf.py`)

```bash
python design_gsmf.py            # base template defaults to ./templet.yml
```
Writes `design_gsmf/run_<id>.yml` (48 complete SWIFT parameter files),
`design_parameters.json` (the sampled values, used later to train the emulator)
and a corner plot. Edit `parameter_limits` / `number_of_simulations` at the top
to change the ranges or sample count. Re-running RESAMPLES and overwrites, so
don't re-run once you've started simulating.

**Stage 2 — build ready-to-submit run directories** (`setup_runs.py`)

```bash
python setup_runs.py \
    --design   design_gsmf \
    --template run_template \     # a working dir of symlinks to swift_mpi, ICs,
    --sbatch   mpirun.sbatch \    #   cooling/yield/xray tables, output_list.txt
    --out      /scratch/gsmf_runs \
    [--dry-run] [--submit] [--force] [--only 0,1,2]
```
For every run it creates `<out>/run_<id>/` containing a copy of `run_template/`
(symlinks are kept AS symlinks, so the big ICs/tables are not duplicated), the
run's parameter file as `flamingo.yml`, and a rendered SLURM script. The
`--sbatch` template may use `{{RUN_ID}}`, `{{RUN_NAME}}`, `{{RUN_DIR}}`,
`{{PARAM_FILE}}` placeholders. With `--submit` it runs `sbatch` from each dir.

> Point `templet.yml`'s input paths at `./<file>` and place matching (absolute)
> symlinks in `run_template/`, so they resolve inside every run directory.

Submit all runs (from the `--out` directory). The subshell `cd` matters: the
sbatch uses relative paths (`./swift_mpi`, `flamingo.yml`, `output-%j.txt`) and
must be submitted from inside each run dir.

```bash
cd /scratch/gsmf_runs
for d in run_*/; do ( cd "$d" && sbatch mpirun.sbatch ); done
```

Equivalently, add `--submit` to the `setup_runs.py` call above to submit as the
directories are created.

Resume after a SLURM time limit. SWIFT writes restart files; add its `-r` flag
before `flamingo.yml` and resubmit every run. This one-liner is idempotent (it
will not add a second `-r` on later resumes):

```bash
cd /scratch/gsmf_runs
for d in run_*/; do grep -q 'swift_mpi.*-r.*flamingo.yml' "$d/mpirun.sbatch" || sed -i 's#\(swift_mpi.*\) flamingo.yml#\1 -r flamingo.yml#' "$d/mpirun.sbatch"; ( cd "$d" && sbatch mpirun.sbatch ); done
```

**Stage 3 — run the galaxy finder** (pgalf, on the cluster)

For each run's z=0 snapshot run pgalf to get a galaxy catalogue
(`FoF_Data/FoF.<NNNNN>/GALCATALOG.LIST.<NNNNN>`). The catalogue is in little-h
units (Msun/h, Mpc/h).

**Stage 4 — measure the GSMF** (`measure_gsmf.py`)

```bash
python measure_gsmf.py \
    --pgalf-repo    /path/to/pgalf_swift_grammar \
    --root-template "/path/to/pgalf/run_{id}" \
    --snapshot      <z=0 snapshot> \
    --design        design_gsmf/design_parameters.json \
    --out-dir       gsmf_measured
```
Bins each catalogue into a GSMF (converting little-h -> physical Msun, Mpc^-3)
and writes one `gsmf_measured/run_<id>.yml` per run in the format
`swiftemulator.io.swift.load_pipeline_outputs` expects. By default M* is the
stellar mass within a **3D 50 proper-kpc aperture** (COLIBRE / observational
convention, read from `GALFIND.DATA`); set `--aperture-kpc` to change the size,
or `--aperture-kpc 0` for the catalogue's total bound stellar mass.

**Stage 5 — train the emulator and calibrate** (to be written)

Load `design_parameters.json` + `gsmf_measured/*.yml`, train a Gaussian process
emulator, then fit to an observed z=0 GSMF (Driver+2022 / Baldry+2012) via an
L2 penalty / MCMC. See the skeleton in [`GSMF_CALIBRATION.md`](GSMF_CALIBRATION.md).

Predicting Simulations
----------------------

The emulator can predict a given scaling relation
(the relationship between two output variables in a simulation, for instance the
masses of galaxies and their size) when varying the underlying physical model
simulated in a continuous way.

As an example from cosmological simulations, imagine varing the energy that supernovae
release when they explode as a parameter `x`. This affects both the sizes and masses of galaxies.
The emulator, using a few 'base' simulations, performed with the real code,
at various values of `x` spaced evenly throughout the viable region, can predict
what the shape of the relationship between mass and size would be at other values
of `x`, given that it has been trained on the base simulation outputs.

Why SWIFT Emulator?
-------------------

The emulator works at a much higher level than other Gaussian Process Regression
libraries, such as `george` (which it is built on).

Working with base machine learning libraries can be tricky, as it typically
requries knowledge of how to structure input arrays in a specific way (both for
training and prediction). They also rarely also include model design routines.
Additionally, validation and visualisation routines
are typically missing.

The SWIFT Emulator package provides a one-stop solution, with a consistent API,
for developing a model design, running it (if using SWIFT), reading in data (
again if using SWIFT), building an emulation model, validating said model,
comparing the model against ground-truth data _across parameter space_ 
(e.g. observations), and visualising the results.

Installation
------------

The package can be installed easily from PyPI under the name `swiftemulator`,
so:

```
pip3 install swiftemulator
```

This will install all necessary dependencies.

The package can be installed from source, by cloning the repository and
then using `pip install -e .` for development purposes.


Requirements
------------

The package requires a number of numerical and experimental design packages.
These have been tested (and are continuously tested) using GitHub actions CI
to use the latest versions available on PyPI. See `requirements.txt` for
details for the packages required to develop SWIFT-Emulator. The packages
will be installed automatically by `pip` when installing from PyPI.


Authors
-------

+ Roi Kugel (@Moyoxkit)
+ Josh Borrow (@jborrow)
