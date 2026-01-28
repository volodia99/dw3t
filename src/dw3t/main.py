import os
import argparse
import json
import importlib
from pathlib import Path

import numpy as np

import toml
from nonos.api import GasDataSet

from dw3t.model import load_model, Opacity
from dw3t._typing import F, FArray2D

def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dw3t",
        description=__doc__,
    )
    parser.suggest_on_error = True  # type: ignore [attr-defined]

    parser.add_argument(
        "parameter_file", 
        type=str, 
        help="work on dw3t parameter file"
    )

    return parser

def main(argv: list[str] | None = None) -> int:
    #TODO: think if CLI or minimalist parameter file? Other option?
    parser = get_parser()
    args = parser.parse_args(argv)
    config = toml.load(args.parameter_file)

    file = config["simulation"]["on"]
    input_dir = config["simulation"]["input_dir"]
    ds = GasDataSet(file, directory=input_dir)

    component = config["simulation"]["component"]
    component = np.atleast_1d(component).tolist()
    if not (set(component) & set(["gas","dust"])):
        raise ValueError(
            f"{component=} should be 'dust', 'gas' or ['dust', 'gas']."
        )

    model = load_model(
        ds=ds,
        unit_length_au=config["simulation"]["unit_length_au"],
        unit_mass_msun=config["simulation"]["unit_mass_msun"],
        component=component,
        config=config,
    )

    if config["simulation"]["processing"] is not None:
        processing_category = config["simulation"]["processing"].split(":")
        if processing_category[0]=="template" and len(processing_category)>1:
            template_modules = []
            for ii in range(1, len(processing_category)):
                try:
                    template_modules.append(importlib.import_module(f"dw3t.template.{processing_category[ii]}", package=None))
                except ModuleNotFoundError:
                    raise ModuleNotFoundError(
                        f"No module named 'template.{processing_category[ii]}'."
                    )
        elif processing_category[0].endswith(".py") and len(processing_category)==1:
            processing_file = processing_category[0]
            if not(os.path.isfile(processing_file)):
                raise FileNotFoundError(f"absolute path of the file '{processing_file}' must exist.")
            spec = importlib.util.spec_from_file_location(Path(processing_file).stem, processing_file)
            template_modules = [importlib.util.module_from_spec(spec)]
            spec.loader.exec_module(template_modules[0])
        else:
            raise ValueError(
                f"{processing_category=} should have the form 'template:func' or path_to_file/file.py."
            )
        for tm in template_modules:
            model = tm.processing(model=model)

    if model.dimension!=3:
        raise ValueError(f"{model.dimension=}D, should be 3D for RADMC3D. Try to post-process the data for it to be 3D.")

    #TODO: Need to improve opacity from .lnk file
    #TODO: test smoothing opacities
    #TODO: Add default config (ChainMap?) + improve config
    model.write_files(
        directory=config["simulation"]["output_dir"],
        write_opacities=True,
        opacity=Opacity(opacity=config["opacity"]["mix"], rho=config["opacity"]["rho"]),
        smoothing=config["opacity"]["smoothing"],
        simulation_files_only=config["simulation"]["simulation_files_only"],
        config=config,
    )

    # Writing config to output directory
    with open(os.path.join(config["simulation"]["output_dir"], "dw3t.full.toml"), 'w') as f:
        output_toml = toml.dump(config, f)

    return 0