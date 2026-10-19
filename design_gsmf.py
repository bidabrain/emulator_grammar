"""
Latin-hypercube parameter sweep for z=0 GSMF calibration of the dual-channel
(thermal + kinetic) SN feedback + spin-jet AGN model (base: flamingo.yml,
25 Mpc box, SNII_energy_fraction_model = ColibrePressure).

Varies the 6 calibration knobs of Table 2 (Chaikin+2026 ColibrePressure f_E):

    EAGLEFeedback:SNII_f_kinetic                        f_kin     kinetic/thermal split
    EAGLEFeedback:SNII_delta_T_K                        dT_SNe    thermal-channel heating dT  (log-sampled)
    EAGLEFeedback:SNII_energy_fraction_min              f_E,min   floor  (low-P stars)
    EAGLEFeedback:SNII_energy_fraction_max              f_E,max   ceiling(high-P stars)
    EAGLEFeedback:SNII_energy_fraction_P_pivot_K_p_cm3  P_pivot   logistic pivot          (log-sampled)
    SPINJETAGN:subgrid_seed_mass_Msun                  M_seed    BH seed mass            (log-sampled)

Everything else in flamingo.yml is held FIXED (sigma_P=0.3, delta_v=50 km/s,
coupling_efficiency, etc.). The f_E model is already ColibrePressure in the
base file, so the PowerLaw knobs (n_0, n_n, Z_0, ...) are inert and untouched.

    pip install swiftemulator          # once
    python design_gsmf.py [/path/to/flamingo.yml]

Produces N SWIFT .yml files under ./design_gsmf/ (base .yml with the 6 swept
keys overwritten), a corner plot, and design_parameters.json (needed later to
train the emulator). Run SWIFT once per file, then measure the z=0 GSMF of each
run and train the emulator on {design_parameters.json, per-run GSMF}.
"""

import sys
import json
from pathlib import Path

from swiftemulator.design import latin
from swiftemulator.io.swift import write_parameter_files
from swiftemulator import ModelSpecification

# Default to the real 25 Mpc config; override with argv[1] if given.
_default_base = Path("./templet.yml")
base_parameter_file = Path(sys.argv[1]) if len(sys.argv) > 1 else _default_base

number_of_simulations = 48          # ~8 x number_of_parameters (6)
output_path = Path("design_gsmf")
output_path.mkdir(exist_ok=True)

# Keys MUST match flamingo.yml (section:key).
#
# f_E,min and f_E,max ranges are kept NON-OVERLAPPING ([.., 0.4] and [0.4, ..])
# so that SNII_energy_fraction_min <= max for EVERY Latin-hypercube sample
# (the ColibrePressure logistic requires min <= max).
#
# delta_T_K, P_pivot and M_seed span orders of magnitude -> sampled uniformly in
# log10 (parameter_limits are log10 values) and transformed to linear on write.
spec = ModelSpecification(
    number_of_parameters=6,
    parameter_names=[
        "EAGLEFeedback:SNII_f_kinetic",
        "EAGLEFeedback:SNII_delta_T_K",                       # log10-sampled
        "EAGLEFeedback:SNII_energy_fraction_min",
        "EAGLEFeedback:SNII_energy_fraction_max",
        "EAGLEFeedback:SNII_energy_fraction_P_pivot_K_p_cm3",  # log10-sampled
        "SPINJETAGN:subgrid_seed_mass_Msun",                  # log10-sampled
    ],
    parameter_printable_names=[
        r"$f_{\rm kin}$",
        r"$\log_{10}\,\Delta T_{\rm SNe}\,[{\rm K}]$",
        r"$f_{\rm E,min}$",
        r"$f_{\rm E,max}$",
        r"$\log_{10}\,P_{\rm pivot}\,[{\rm K\,cm^{-3}}]$",
        r"$\log_{10}\,M_{\rm seed}\,[M_\odot]$",
    ],
    parameter_limits=[
        [0.0, 1.0],      # f_kin        (current 0.3)      0=pure thermal, 1=pure kinetic
        [6.5, 8.0],      # log10 dT_SNe (current ~10^6.54; EAGLE 10^7.5)
        [0.05, 0.4],     # f_E,min      (current 0.1)      floor,   <= f_E,max
        [0.4, 8.0],      # f_E,max      (current 0.4)      ceiling, >= f_E,min
        [2.5, 4.5],      # log10 P_pivot(current 8e3 -> 10^3.9; COLIBRE ~10^3.4)
        [4.0, 6.0],      # log10 M_seed (current 1e5 -> 10^5)
    ],
)

# Transform keys are the FULL "section:key" (write_parameter_files matches the
# full name, not the bare key).
parameter_transforms = {
    "EAGLEFeedback:SNII_delta_T_K": lambda x: 10.0 ** x,
    "EAGLEFeedback:SNII_energy_fraction_P_pivot_K_p_cm3": lambda x: 10.0 ** x,
    "SPINJETAGN:subgrid_seed_mass_Msun": lambda x: 10.0 ** x,
}

model_parameters = latin.create_hypercube(
    model_specification=spec,
    number_of_samples=number_of_simulations,
)

model_parameters.plot_model(
    model_specification=spec,
    filename=output_path / "design_gsmf_corner.png",
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

# Saved in SAMPLING space (log10 for the 3 log-sampled params). Keep the SAME
# space when training the emulator so predictions map back correctly.
with open(output_path / "design_parameters.json", "w") as fh:
    json.dump(model_parameters.model_parameters, fh, indent=2)

print(f"Base config : {base_parameter_file}")
print(f"Wrote {number_of_simulations} parameter files + corner plot to {output_path}/")
print("Next: run SWIFT per run_<id>.yml; params saved in design_parameters.json")
