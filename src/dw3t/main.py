import os
import argparse
import json
import importlib
from pathlib import Path

import numpy as np

import tomli
import tomli_w
from deep_chainmap import DeepChainMap
from nonos.api import GasDataSet

from dw3t.model import load_model, Opacity
from dw3t._typing import F, FArray2D
from dw3t._parsing import is_set, list_of_middle_keys
from dw3t.default import DEFAULT_LAYER, MANDATORY_SET

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
    parser = get_parser()
    args = parser.parse_args(argv)
    with open(args.parameter_file, "rb") as f:
        config_file_layer = tomli.load(f)

    #TODO: improve handling of mandatory parameters
    # ensures that we do not add mandatory parameters without knowing it
    expected_length_mandatory = 5
    if len(MANDATORY_SET)!=expected_length_mandatory:
        raise ValueError(
            f"expecting {expected_length_mandatory} mandatory parameters, not {len(MANDATORY_SET)}."
        )
    if "gas" in config_file_layer["simulation"]["component"]:
        MANDATORY_SET.update(["species", "abundance"])
    if "dust" in config_file_layer["simulation"]["component"]:
        MANDATORY_SET.update(["opacity"])
    # ensures that all mandatory parameters are defined in toml file
    if MANDATORY_SET.difference(set(list_of_middle_keys(config_file_layer))):
        raise ValueError(
            f"at least one mandatory parameter is missing: {MANDATORY_SET.difference(set(list_of_middle_keys(config_file_layer)))}"
        )

    config = DeepChainMap(config_file_layer, DEFAULT_LAYER)
    if not is_set(config["stars"]["M_star"]):
        config["stars"]["M_star"] = config["simulation"]["unit_mass_msun"]

    component = config["simulation"]["component"]
    component = np.atleast_1d(component).tolist()
    if not (set(component) & set(["gas","dust"])):
        raise ValueError(
            f"{component=} should be 'dust', 'gas' or ['dust', 'gas']."
        )

    model = load_model(
        unit_length_au=config["simulation"]["unit_length_au"],
        unit_mass_msun=config["simulation"]["unit_mass_msun"],
        component=component,
        config=config,
    )

    if "processing" in config["simulation"]:
        processing_dict = np.atleast_1d(config["simulation"]["processing"]).tolist()
        processing_category = [d.get("mode") for d in processing_dict]
        if "userdef" in processing_category:
            if len(processing_category)!=1:
                raise NotImplementedError(
                    "mode='userdef' cannot yet be combined with other processing modes."
                )
            processing_file = processing_dict[0].get("file")
            if not processing_file.endswith(".py"):
                raise ValueError(
                    f"{proceprocessing_file=} should be a .py file. See example in dw3t/src/dw3t/template/userdef.py."
                )
            if not(os.path.isfile(processing_file)):
                raise FileNotFoundError(f"absolute path of the file '{processing_file}' must exist.")
            spec = importlib.util.spec_from_file_location(Path(processing_file).stem, processing_file)
            template_modules = [importlib.util.module_from_spec(spec)]
            spec.loader.exec_module(template_modules[0])
            kwargs = [processing_dict[0].copy()]
            del kwargs[0]["mode"]
        else: 
            if "nonos" in processing_category:
                if not is_set(config["simulation"]["output_dir"]):
                    config["simulation"]["output_dir"] = os.path.join(
                        processing_dict[processing_category.index("nonos")]["input_dir"], 
                        "radmc3d"
                    )
            template_modules = []
            kwargs = []
            for ii in range(len(processing_category)):
                try:
                    template_modules.append(importlib.import_module(f"dw3t.template.{processing_category[ii]}", package=None))
                    kwargs.append(processing_dict[ii].copy())
                    del kwargs[ii]["mode"]
                except ModuleNotFoundError:
                    raise ModuleNotFoundError(
                        f"No module named 'template.{processing_category[ii]}'."
                    )
        for tm, kw in zip(template_modules, kwargs):
            model = tm.processing(model=model, kwargs=kw)

    if not is_set(config["simulation"]["output_dir"]):
        MANDATORY_SET.update(["output_dir"])

    if MANDATORY_SET.difference(set(list_of_middle_keys(config_file_layer))):
        raise ValueError(
            f"at least one mandatory parameter is missing: {MANDATORY_SET.difference(set(list_of_middle_keys(config_file_layer)))}"
        )

    if model.dimension!=3:
        raise ValueError(f"{model.dimension=}D, should be 3D for RADMC3D. Try to post-process the data for it to be 3D.")

    #TODO: test smoothing opacities
    if write_opacities:=("dust" in config["simulation"]["component"]):
        opacity = Opacity(
            mix=config["dust"]["opacity"]["mix"], 
            rho=config["dust"]["opacity"]["rho"]
        )
    else:
        opacity = None

    model.write_files(
        directory=config["simulation"]["output_dir"],
        write_opacities=write_opacities,
        opacity=opacity,
        smoothing=config["dust"]["opacity"]["smoothing"],
        simulation_files_only=config["simulation"]["simulation_files_only"],
        config=config,
    )

    # Writing config to output directory
    with open(os.path.join(config["simulation"]["output_dir"], "dw3t.full.toml"), "wb") as f:
        tomli_w.dump(config, f)

    return 0