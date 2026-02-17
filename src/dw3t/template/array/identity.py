from dw3t.model import Array

def processing(*, array:"Array", kwargs:dict) -> "Array":
    array = Array(
        value=array.value,
        unit_length_au=array.unit_length_au,
    )
    return array