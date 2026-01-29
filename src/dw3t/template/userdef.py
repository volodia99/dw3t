import numpy as np
import astropy.units as u

from dw3t.model import Model, Grid, Gas, Dust
from nonos._geometry import Axis

def processing(*, model:"Model", kwargs:dict) -> "Model":
    if model.dimension!=2:
        raise ValueError(f"phi_expansion works in 2D, model is {model.dimension}D.")
    if model.geometry!="spherical":
        raise ValueError(f"phi_expansion works if the model geometry is spherical.")
    if Axis.AZIMUTH in model.reduced_axes:
        raise ValueError(f"phi_expansion works if there is ONE cell only in the azimuthal direction.")

    nphi = 128
    phi = np.linspace(0, 2.*np.pi, nphi+1) * u.radian

    gas = None
    dust = None
    if model.gas is not None:
        gas = Gas(
            rho = np.repeat(model.gas.rho, nphi, axis=2),
            v1 = np.repeat(model.gas.v1, nphi, axis=2),
            v2 = np.repeat(model.gas.v2, nphi, axis=2),
            v3 = np.repeat(model.gas.v3, nphi, axis=2),
        )
    if model.dust is not None:
        mask_inner_edge = (model.grid.x1c < 1.5*model.grid.x1c.min())
        dustrhog = model.dust.rho
        dustrhog[mask_inner_edge] = 0.0#np.nan
        dust = Dust(
            rho = np.repeat(dustrhog, nphi, axis=2),
            size = model.dust.size,
        )

    model = Model(
        grid=Grid(
            x1 = model.grid.x1,
            x2 = model.grid.x2,
            x3 = phi,
        ),
        gas=gas,
        dust=dust,
        unit_length_au=model.unit_length_au,
        unit_mass_msun=model.unit_mass_msun,
        component=model.component,
        geometry=model.geometry,
    )
    return model
