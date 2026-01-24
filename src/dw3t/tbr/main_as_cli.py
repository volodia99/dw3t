import os
import argparse
import json
import importlib

import numpy as np

from nonos.api import GasDataSet

from dw3t.model import load_model
from dw3t._typing import F, FArray2D

def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dw3t",
        description=__doc__,
    )
    parser.suggest_on_error = True  # type: ignore [attr-defined]

    subparsers = parser.add_subparsers(
        help="Choice between from_on (from output number) and from_path (from absolute path to file name)", 
        dest="input_processing"
    )

    ## Create the parser for the "on" command
    parser_on = subparsers.add_parser(
        "from_on", 
        help="Work with output number approach. Try 'dw3t from_on -h' for more info."
    )
    parser_on.add_argument(
        "file_on", 
        type=int, 
        help="Get simulation output from its output number"
    )

    parser_on.add_argument(
        "-dir",
        type=str,
        required=True,
        dest="directory",
        help="required: location of the simulated output, if described by its output number.",
    )

    ## Create the parser for the "name" command
    parser_file = subparsers.add_parser(
        "from_path", 
        help="Work with file name approach. Try 'dw3t from_path -h' for more info."
    )
    parser_file.add_argument(
        "file_path", 
        type=str, 
        help="Get simulation output from its name"
    )

    for subparser in [parser_on, parser_file]:
        subparser.add_argument(
            "-unit_length_au",
            type=float,
            default=None,
            required=True,
            help="required: code unit of length [au].",
        )

        subparser.add_argument(
            "-unit_mass_msun",
            type=float,
            default=None,
            required=True,
            help="required: code unit of mass [solMass].",
        )

        subparser.add_argument(
            "-processing",
            type=str,
            default=None,
        )

        subparser.add_argument(
            "-component",
            type=str,
            required=True,
            default=None,
        )

        subparser.add_argument(
            "-simulation_files_only",
            action="store_true",
        )

    return parser

def main(argv: list[str] | None = None) -> int:
    #TODO: think if CLI or minimalist parameter file? Other option?
    parser = get_parser()
    args = parser.parse_args(argv)

    if args.input_processing=="from_path":
        file = args.file_path
        directory = None
        if not(os.path.isfile(file)):
            raise FileNotFoundError(f"absolute path of the file '{file}' must exist.")
        ds = GasDataSet(file)
    elif args.input_processing=="from_on":
        file = args.file_on
        directory = args.directory
        ds = GasDataSet(file, directory=directory)

    component = args.component
    if component not in ("dust", "gas", "dustgas", "gasdust"):
        raise ValueError(
            f"{component=} should be 'dust', 'gas', 'dustgas' (or equivalently 'gasdust')."
        )
    if component=="gasdust":
        component="dustgas"

    model = load_model(
        ds=ds,
        unit_length_au=args.unit_length_au,
        unit_mass_msun=args.unit_mass_msun,
        component=args.component,
    )

    if args.processing is not None:
        processing_category = args.processing.split(":")
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
            #TODO: implement userdef processing
            raise NotImplementedError(
                f"userdef processing not yet implemented."
            )
            processing_file = processing_category[0]
            if not(os.path.isfile(processing_file)):
                raise FileNotFoundError(f"absolute path of the file '{processing_file}' must exist.")
            template_modules = [importlib.import_module(f"{processing_file}", package=None)]
        else:
            raise ValueError(
                f"{processing_category=} should have the form 'template:func' or path_to_file/file.py."
            )
        for tm in template_modules:
            model = tm.processing(model=model)

    # dargs = vars(args)
    # # Writing CLI args to a JSON file
    # with open(os.path.join(args.prodimo_model_directory, "dw3t_args.json"), "w") as outfile:
    #     json.dump(dargs, outfile, indent=1)

    if model.dimension!=3:
        raise ValueError(f"{model.dimension=}D, should be 3D for RADMC3D. Try to post-process the data for it to be 3D.")

    #TODO: Need to check how to activate smoothing
    #TODO: Need to check opacity from .lnk file
    #TODO: Add RT directory 
    model.write_files(
        directory="/home/gwf/projects/dw3t/radmc3d_test",
        write_opacities=True,
        opacity=None,
        smoothing=False,
        simulation_files_only=args.simulation_files_only,
    )

    #TODO: test RADMC3D (mctherm + dust RT + gas RT)
    return 0