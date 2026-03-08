from typing import Any

def is_set(x: Any) -> bool:
    """
    Adapted from nonos_cli (https://github.com/la-niche/nonos/tree/main/cli)
    """
    return x!="unset"

def list_of_middle_keys(dictionary: dict) -> list:
    list_all_middle_keys = []
    for value in dictionary.values():
        if isinstance(value, dict):
            list_all_middle_keys.extend(value.keys())
    return list_all_middle_keys