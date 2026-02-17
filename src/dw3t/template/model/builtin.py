import os

import numpy as np
import astropy.units as u
import astropy.constants as uc

import inifix
from nonos.api import GasDataSet
from dw3t._typing import FArray1D
from dw3t.model import Model, Grid, Gas, Dust

def computeSizeMM(betai:FArray1D, *, rhoint:float, unit_length_au:float, unit_mass_msun:float) -> u.Quantity:
    sizeMM_k = betai/(np.sqrt(np.pi/8.0)*rhoint.to(u.g/u.cm/u.cm/u.cm)*(unit_length_au.to(u.cm))**2/unit_mass_msun.to(u.g))
    return sizeMM_k.to(u.mm)

def processing(*, model:"Model", kwargs:dict) -> "Model":
    file = kwargs["input_number"]
    input_dir = kwargs["input_dir"]
    ds = GasDataSet(file, directory=input_dir)

    if ds.native_geometry!=model.geometry:
        raise NotImplementedError(f"native_geometry='{ds.native_geometry}' should be {model.geometry}")

    UNIT_DENSITY = model.unit_mass_msun / model.unit_length_au**3
    UNIT_VELOCITY = np.sqrt(uc.G*model.unit_mass_msun/model.unit_length_au).to(u.m/u.s)

    grid = Grid(
            x1 = ((ds.coords.get_axis_array("r") * model.unit_length_au).to(u.cm)),#.value,
            x2 = ds.coords.get_axis_array("theta") * u.radian,
            x3 = ds.coords.get_axis_array("phi") * u.radian,
        )

    if "gas" in model.component:
        #TODO: add flexibility
        print(f"WARNING: Assuming no omegraframe.")
        gas = Gas(
            rho = ((ds["RHO"].data * UNIT_DENSITY).to(u.g / u.cm**3)),#.value,
            v1 = ((ds["VX1"].data * UNIT_VELOCITY).to(u.cm / u.s)),#.value,
            v2 = ((ds["VX2"].data * UNIT_VELOCITY).to(u.cm / u.s)),#.value,
            v3 = ((ds["VX3"].data * UNIT_VELOCITY).to(u.cm / u.s)),#.value,
        )
    if "dust" in model.component:
        print(f"WARNING: 'dust' not implemented in a general way with nonos. Implementation specific to IDEFIX.")
        directory = ds._parameters_input["directory"]
        rhoint_csg = kwargs["internal_rho"]*(u.g/u.cm/u.cm/u.cm)
        inifile = inifix.load(os.path.join(directory, "idefix.ini"))
        dragType = inifile["Dust"]["drag"][0]
        if dragType!="size":
            raise ValueError(f"{dragType=} should be 'size'.")
        dustBeta = np.array(inifile["Dust"]["drag"][1:])
        dustSize = computeSizeMM(
            dustBeta, 
            rhoint=rhoint_csg,
            unit_length_au=model.unit_length_au,
            unit_mass_msun=model.unit_mass_msun,
            )#.value
        dustSize_ascending_indices = np.argsort(dustSize)
        dustrho = np.empty(grid.shape+(len(dustSize),))
        for kk in dustSize_ascending_indices:
            dustrho[..., kk] = ds[f"DUST{kk}_RHO"].data
        dust = Dust(
            rho = ((dustrho * UNIT_DENSITY).to(u.g / u.cm**3)),#.value,
            size = dustSize,
        )

    updated_model = model
    updated_model.grid = grid
    updated_model.gas = gas
    updated_model.dust = dust

    return updated_model
