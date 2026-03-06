from dataclasses import replace
import numpy as np
from dw3t.model import Model

def processing(*, model:"Model", kwargs:dict) -> "Model":
    inner_factor = kwargs["inner_factor"]
    r_spherical = model.grid.x1c
    imask, = np.where(r_spherical < inner_factor * r_spherical.min())

    gas = model.gas
    dust = model.dust
    if "gas" in model.component:
        rhogas, v1, v2, v3 = (gas.rho, gas.v1, gas.v2, gas.v3)
        for dfield in (rhogas, v1, v2, v3):
            dfield[imask,:,:] = 0.0
        gas = replace(
            gas,
            rho=rhogas,
            v1=v1,
            v2=v2,
            v3=v3,
        )
    if "dust" in model.component:
        rhodust = dust.rho
        rhodust[imask,:,:] = 0.0
        dust = replace(
            dust, 
            rho=rhodust,
        )
    return replace(
        model,
        gas=gas,
        dust=dust,
    )
