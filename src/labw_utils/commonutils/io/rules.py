"""
rules.py -- Rule definition for Rule-based IO proxy.
"""

from __future__ import annotations

__all__ = (
    "RuleType",
    "OpenerType",
    "BaseFileRule",
    "FileRuleRing"
)

import bz2
import gzip
import lzma

from labw_utils.commonutils.io import PathType, FDType, convert_path_to_str
from labw_utils.typing_importer import Callable, Type, List


RuleType: Type[Callable]
"""
Type definition of rules.

Is a function which accepts a path as input and give result (True/False) as output.
"""

OpenerType:Type[Callable]
"""
Type definition of openers.

Is a function which accepts a path and (optionally) kwyword arguments as input and give
file descriptor as output.
"""

# This goes wrong in Python 3.7
try:
    RuleType = Callable[[PathType, ...], bool]
except TypeError:
    RuleType = Callable[..., bool]
try:
    OpenerType = Callable[[PathType, ...], FDType]
except TypeError:
    OpenerType = Callable[..., FDType]


class BaseFileRule:
    """
    Base file rule type. Is an association between rule and opener.
    """
    _rule: RuleType
    _opener: OpenerType
    _rule_name: str

    @property
    def rule_name(self) -> str:
        """Name of the rule"""
        return self._rule_name

    def check_applicability(self, path: PathType, *args, **kwargs) -> bool:
        """Application of rule which checks whether a file can be opened by the associated opener"""
        return self._rule(path, *args, **kwargs)

    def open(self, path: PathType, *args, **kwargs) -> FDType:
        """Open that file."""
        return self._opener(path, *args, **kwargs)

    def __init__(self, rule: RuleType, opener: OpenerType, rule_name: str = "unnames_rule"):
        self._rule = rule
        self._opener = opener
        self._rule_name = rule_name

    @classmethod
    def from_extension(cls, *extensions: str, opener: OpenerType, rule_name: str = "unnamed_rule"):
        """
        Initialize a new rule from extension name.

        :param extensions: List of valid extensions.
        :param opener: Associated opener.
        :param rule_name: Name of the rule.
        """

        def extension_rule(path: PathType, *args, **kwargs):
            _, _ = args, kwargs
            path = convert_path_to_str(path)
            for extension in extensions:
                if path.endswith(extension):
                    return True
            return False

        new_instance = cls(extension_rule, opener=opener, rule_name=rule_name)
        return new_instance


class FileRuleRing:
    """
    The rule ring.
    """
    _registered_rules: List[BaseFileRule] = []

    @staticmethod
    def register(rule: BaseFileRule):
        """
        Add a new rule.

        :param rule: The rule you wish to add.
        """
        FileRuleRing._registered_rules.append(rule)

    @staticmethod
    def open(path: str, *args, **kwargs) -> FDType:
        """
        Open a file. If no rule was identified, would fallback to :py:func:`open`.
        """
        for rule in FileRuleRing._registered_rules:
            if rule.check_applicability(path, *args, **kwargs):
                return rule.open(path, *args, **kwargs)
        return open(path, *args, **kwargs)


FileRuleRing.register(
    BaseFileRule.from_extension(".gz", ".GZ", opener=gzip.open, rule_name="gzip_suffix_rule")
)
FileRuleRing.register(
    BaseFileRule.from_extension(".xz", ".lzma", opener=lzma.open, rule_name="lzma_suffix_rule")
)
FileRuleRing.register(
    BaseFileRule.from_extension(".bz2", opener=bz2.open, rule_name="bz2_suffix_rule")
)
