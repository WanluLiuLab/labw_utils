"""
TODO docs

.. versionadded:: 1.0.3
"""

import enum

from labw_utils.typing_importer import Any


class EnhancedEnum(enum.Enum):
    """
    TODO docs

    .. versionadded:: 1.0.3
    """

    @classmethod
    def from_name(cls, in_name: str):
        for choices in cls:
            if choices.name == in_name:
                return choices
        raise TypeError

    @classmethod
    def from_value(cls, in_value: Any):
        for choices in cls:
            if choices.value == in_value:
                return choices
        raise TypeError
