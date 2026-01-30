from typing import Any

def is_set(x: Any) -> bool:
    """
    Adapted from nonos_cli (https://github.com/la-niche/nonos/tree/main/cli)
    """
    return x!="unset"

def list_of_middle_keys(dictionary: dict) -> list:
    list_all_middle_keys = []
    for key, value in dictionary.items():
        if isinstance(value, dict):
            for element in list(value.keys()):
                list_all_middle_keys.append(element)
    return list_all_middle_keys