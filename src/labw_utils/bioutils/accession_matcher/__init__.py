from __future__ import annotations

from abc import abstractmethod, ABC
from typing import Mapping, Union, List, Optional, Type

_AccessionMatchResultValueType = Union[str, Mapping[str, str]]


class AccessionMatchResult:
    _toplevel: str
    _details: Mapping[str, AccessionMatchResult]

    def __init__(self, toplevel: str, details: Mapping[str, Union[str, AccessionMatchResult]]):
        self._toplevel = toplevel
        self._details = details

    def as_dict(self) -> Mapping[str, _AccessionMatchResultValueType]:
        return {
            "toplevel": self._toplevel,
            "details": {
                k: v.as_dict() if isinstance(v, AccessionMatchResult) else v
                for k, v in self._details.items()
            }
        }

    def __repr__(self) -> str:
        return f"{self._toplevel}, with {self._details}"


class AccessionMatcherRuleType(ABC):

    @abstractmethod
    def match(self, accession: str) -> Optional[AccessionMatchResult]:
        raise NotImplementedError


class ChainAccessionMatcherRuleType(AccessionMatcherRuleType):
    _rule_chain: List[Type[AccessionMatcherRuleType]]

    def match(self, accession: str) -> Optional[AccessionMatchResult]:
        for rule in self._rule_chain:
            match_result = rule().match(accession)
            if match_result is not None:
                return match_result
        return None
