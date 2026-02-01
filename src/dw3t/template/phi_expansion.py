import numpy as np
import astropy.units as u

from dw3t.model import Model, Grid, Gas, Dust
from nonos._geometry import Axis

def processing(*, model:"Model", kwargs:dict) -> "Model":
    #TODO: to be improved?
    if model.dimension!=2:
        raise ValueError(f"phi_expansion works in 2D, model is {model.dimension}D.")
    if model.geometry!="spherical":
        raise ValueError(f"phi_expansion works if the model geometry is spherical.")
    if Axis.AZIMUTH in model.reduced_axes:
        raise ValueError(f"phi_expansion works if there is ONE cell only in the azimuthal direction.")

    nphi = kwargs["nphi"]
    phi = np.linspace(0, 2.*np.pi, nphi+1) * u.radian

    grid=Grid(
        x1 = model.grid.x1,
        x2 = model.grid.x2,
        x3 = phi,
    )

    gas = None
    dust = None
    if "gas" in model.component:
        gas = Gas(
            rho = np.repeat(model.gas.rho, nphi, axis=2),
            v1 = np.repeat(model.gas.v1, nphi, axis=2),
            v2 = np.repeat(model.gas.v2, nphi, axis=2),
            v3 = np.repeat(model.gas.v3, nphi, axis=2),
        )
    if "dust" in model.component:
        dust = Dust(
            rho = np.repeat(model.dust.rho, nphi, axis=2),
            size = model.dust.size,
        )

    updated_model = model
    updated_model.grid = grid
    updated_model.gas = gas
    updated_model.dust = dust
    
    return updated_model
