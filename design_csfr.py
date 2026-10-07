"""
Design a Latin-hypercube parameter sweep for CSFR calibration.

Produces N SWIFT parameter files (one per run) under ./design_csfr/,
each = your base .yml with the swept subgrid parameters overwritten,
plus a corner plot of the design.

Usage:
    pip install swiftemulator          # once
    python design_csfr.py /path/to/base_eagle_spin_jet.yml

Then run SWIFT once per generated .yml, keep each run's SFR.txt.
"""

import sys
from pathlib import Path

from swiftemulator.design import latin
from swiftemulator.io.swift import write_parameter_files
from swiftemulator import ModelSpecification

# ---------------------------------------------------------------
# 1) base parameter file (your working eagle_spin_jet .yml)
# ---------------------------------------------------------------
base_parameter_file = Path(sys.argv[1]) if len(sys.argv) > 1 else None

# ---------------------------------------------------------------
# 2) how many simulations in the sweep (compute-budget dependent)
# ---------------------------------------------------------------
number_of_simulations = 24        # TODO: 16-32 typical; >= ~5x number_of_parameters
output_path = Path("design_csfr")
output_path.mkdir(exist_ok=True)

# ---------------------------------------------------------------
# 3) parameters to vary + ranges.
#    For CSFR the dominant knobs are STELLAR (SNII) feedback;
#    AGN matters for the low-z decline. Names MUST match keys in
#    your SWIFT .yml exactly (section:key). Edit to your model.
# ---------------------------------------------------------------
spec = ModelSpecification(
    number_of_parameters=4,
    parameter_names=[
        "EAGLEFeedback:SNII_energy_fraction_min",       # TODO: confirm vs your .yml
        "EAGLEFeedback:SNII_energy_fraction_max",
        "EAGLEFeedback:SNII_energy_fraction_n_0_H_p_cm3",  # log-sampled (see transform)
        "EAGLEAGN:coupling_efficiency",                 # TODO: AGN knob name in spin_jet model
    ],
    parameter_printable_names=[
        r"$f_{\rm E,min}$",
        r"$f_{\rm E,max}$",
        r"$\log_{10} n_{\rm H,0}$",
        r"$\epsilon_{\rm AGN}$",
    ],
    parameter_limits=[
        [0.0, 1.0],     # TODO: set sensible ranges for your model
        [1.0, 7.0],
        [-1.0, 1.5],    # sampled in log10, converted below
        [0.0, 0.1],
    ],
)

# parameters sampled in log and written as linear (key -> transform)
parameter_transforms = {
    "SNII_energy_fraction_n_0_H_p_cm3": lambda x: 10.0**x,
}

# ---------------------------------------------------------------
# 4) build the Latin hypercube and write one .yml per run
# ---------------------------------------------------------------
model_parameters = latin.create_hypercube(
    model_specification=spec,
    number_of_samples=number_of_simulations,
)

model_parameters.plot_model(
    model_specification=spec,
    filename=output_path / "design_csfr_corner.png",
    corner_kwargs=dict(bins=5),
)

write_parameter_files(
    filenames={
        key: output_path / f"run_{key}.yml"
        for key in model_parameters.model_parameters.keys()
    },
    model_parameters=model_parameters,
    parameter_transforms=parameter_transforms,
    base_parameter_file=base_parameter_file,
)

# save the parameter table for later (emulator training needs it)
import json
with open(output_path / "design_parameters.json", "w") as fh:
    json.dump(model_parameters.model_parameters, fh, indent=2)

print(f"Wrote {number_of_simulations} parameter files + corner plot to {output_path}/")
print("Each run_<id>.yml -> run SWIFT, keep its SFR.txt. Params saved in design_parameters.json")
