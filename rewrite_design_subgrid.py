"""Run from emulator_grammar/: python rewrite_design_subgrid.py

Re-write the SAME 6 design points (design_gsmf/design_parameters.json) on top of
the current templet.yml (Subgrid SF threshold) into design_gsmf_subgrid/."""
import json, shutil
from pathlib import Path
from swiftemulator import ModelParameters
from swiftemulator.io.swift import write_parameter_files
src, out = Path("design_gsmf"), Path("design_gsmf_subgrid")
out.mkdir(exist_ok=False)
params = json.load(open(src / "design_parameters.json"))
tr = {k: (lambda x: 10.0 ** x) for k in ["EAGLEFeedback:SNII_delta_T_K",
      "EAGLEFeedback:SNII_energy_fraction_P_pivot_K_p_cm3", "SPINJETAGN:subgrid_seed_mass_Msun"]}
write_parameter_files(filenames={k: out / f"run_{k}.yml" for k in params},
                      model_parameters=ModelParameters(model_parameters=params),
                      parameter_transforms=tr, base_parameter_file=Path("templet.yml"))
shutil.copyfile(src / "design_parameters.json", out / "design_parameters.json")
print("wrote", sorted(p.name for p in out.iterdir()))
