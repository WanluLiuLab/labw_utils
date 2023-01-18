from __future__ import annotations
import base64
import random

from typing import Any, Dict, Iterable, Optional


class MyItem:
    __slots__ = (
        "_field_i1",
        "_field_i2",
        "_field_op_b",
        "_field_s",
        "_field_op_d",
    )

    _field_i1:int
    _field_i2:int
    _field_op_b:Optional[bool]
    _field_s:str
    _field_op_d:Optional[Dict[Any, Any]]

    @property
    def i1(self) -> int:
        return self._field_i1
    
    @property
    def i2(self) -> int:
        return self._field_i2

    @property
    def b(self) -> Optional[bool]:
        return self._field_op_b
    
    @property
    def s(self) -> str:
        return self._field_s
    
    @property
    def d(self) -> Optional[Dict[Any, Any]]:
        return self._field_op_d
    
    def __init__(
        self,
        i1:int,
        i2:int, 
        b:Optional[bool],
        s:str,
        d:Optional[Dict[Any, Any]]
    ) -> None:
        self._field_i1 = i1
        self._field_i2 = i2
        self._field_op_b = b
        self._field_s = s
        self._field_op_d = d

    @classmethod
    def blank_none(cls) -> MyItem:
        return cls(0, 0, None, "", None)

    @classmethod
    def blank(cls) -> MyItem:
        return cls(0, 0, False, "", {})


def generate_random_items(num_items:int) -> Iterable[MyItem]:
    for _ in range(num_items):
        yield MyItem(
            random.randint(1, 100),
            random.randint(120, 12000),
            False,
            str(base64.b64encode(random.randbytes(1200)), encoding="UTF-8"),
            {i:random.random() for i in range(1024)}
        )
