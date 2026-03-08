from dataclasses import replace

import numpy as np
import astropy.units as u

from dw3t.model import Model
from nonos._geometry import Axis

def processing(*, model:"Model", kwargs:dict) -> "Model":
    #TODO: to be improved?
    if model.dimension!=2:
        raise ValueError(f"phi_expansion works in 2D, model is {model.dimension}D.")
    if model.grid.geometry!="spherical":
        raise ValueError("phi_expansion works if the model geometry is spherical.")
    if Axis.AZIMUTH not in model.reduced_axes:
        raise ValueError("phi_expansion works if there is ONE cell only in the azimuthal direction.")

    nphi = kwargs["nphi"]
    phi = np.linspace(0, 2.*np.pi, nphi+1) * u.radian

    grid = model.grid
    gas = model.gas
    dust = model.dust

    grid = replace(
        grid,
        x1=model.grid.x1,
        x2=model.grid.x2,
        x3=phi,
        geometry=model.grid.geometry,
    )

    if "gas" in model.component:
        gas = replace(
            gas,
            rho = np.repeat(model.gas.rho, nphi, axis=2),
            v1 = np.repeat(model.gas.v1, nphi, axis=2),
            v2 = np.repeat(model.gas.v2, nphi, axis=2),
            v3 = np.repeat(model.gas.v3, nphi, axis=2),
        )
    if "dust" in model.component:
        dust = replace(
            dust,
            rho = np.repeat(model.dust.rho, nphi, axis=2),
            size = model.dust.size,
        )

    return replace(
        model,
        grid=grid,
        gas=gas,
        dust=dust,
    )
