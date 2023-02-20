from typing import Union, Tuple


def scale_si(num: Union[int, float]) -> Tuple[Union[float, int], str]:
    prefix = ""
    if num > 1024:
        num = num / 1024
        prefix = "Ki"
    if num > 1024:
        num = num / 1024
        prefix = "Mi"
    if num > 1024:
        num = num / 1024
        prefix = "Gi"
    if num > 1024:
        num = num / 1024
        prefix = "Ti"
    if num > 1024:
        num = num / 1024
        prefix = "Pi"
    if num > 1024:
        num = num / 1024
        prefix = "Ei"
    return num, prefix
