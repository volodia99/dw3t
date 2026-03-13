from dw3t.model import Model

def processing(*, model:"Model", **_kwargs) -> "Model":
    model = Model(
        grid=model.grid,
        gas=model.gas,
        dust=model.dust,
        unit_length_au=model.unit_length_au,
        unit_mass_msun=model.unit_mass_msun,
        component=model.component,
    )
    return model