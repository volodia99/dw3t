import os
import sys

import numpy as np
import astropy.units as u
import astropy.constants as uc

import inifix
from nonos.api import GasDataSet
from dw3t._typing import FArray1D
from dw3t.model import Model, Grid, Gas, Dust

if sys.version_info >= (3, 13):
    from copy import replace
else:
    from dataclasses import replace

def computeSizeMM(betai:FArray1D, *, rhoint:float, unit_length_au:float, unit_mass_msun:float) -> u.Quantity:
    sizeMM_k = betai/(np.sqrt(np.pi/8.0)*rhoint.to(u.g/u.cm/u.cm/u.cm)*(unit_length_au.to(u.cm))**2/unit_mass_msun.to(u.g))
    return sizeMM_k.to(u.mm)

def processing(*, model:"Model", **kwargs) -> "Model":
    file = kwargs["input_number"]
    input_dir = kwargs["input_dir"]
    ds = GasDataSet(file, directory=input_dir)

    unit_length_au = model.unit_length_au
    unit_mass_msun = model.unit_mass_msun
    UNIT_DENSITY = unit_mass_msun / unit_length_au**3
    UNIT_VELOCITY = np.sqrt(uc.G*unit_mass_msun/unit_length_au).to(u.m/u.s)

    model = replace(
        model,
        grid=Grid(
            x1=((ds.coords.get_axis_array("r") * unit_length_au).to(u.cm)),#.value,
            x2=ds.coords.get_axis_array("theta") * u.radian,
            x3=ds.coords.get_axis_array("phi") * u.radian,
            geometry=ds.native_geometry,
        ),
    )

    nphi = 128
    phi = np.linspace(0, 2.*np.pi, nphi+1) * u.radian

    gas = None
    dust = None
    if "gas" in model.component:
        #TODO: add flexibility
        print("WARNING: Assuming no omegraframe.")
        gas = Gas(
            rho = ((np.repeat(ds["RHO"].data, nphi, axis=2) * UNIT_DENSITY).to(u.g / u.cm**3)),
            v1 = ((np.repeat(ds["VX1"].data, nphi, axis=2) * UNIT_VELOCITY).to(u.cm / u.s)),
            v2 = ((np.repeat(ds["VX2"].data, nphi, axis=2) * UNIT_VELOCITY).to(u.cm / u.s)),
            v3 = ((np.repeat(ds["VX3"].data, nphi, axis=2) * UNIT_VELOCITY).to(u.cm / u.s)),
        )

    if "dust" in model.component:
        print("WARNING: 'dust' not implemented in a general way with nonos. Implementation specific to IDEFIX.")
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
            unit_length_au=unit_length_au,
            unit_mass_msun=unit_mass_msun,
            )#.value
        dustSize_ascending_indices = np.argsort(dustSize)
        dustrho = np.empty((*model.grid.shape, len(dustSize)))
        for kk in dustSize_ascending_indices:
            dustrho[..., kk] = ds[f"DUST{kk}_RHO"].data
        dust = Dust(
            rho = ((np.repeat(dustrho, nphi, axis=2) * UNIT_DENSITY).to(u.g / u.cm**3)),
            size = dustSize,
        )

    return replace(
        model,
        grid=replace(
            model.grid,
            x3=phi,
        ),
        gas=gas,
        dust=dust,
    )
