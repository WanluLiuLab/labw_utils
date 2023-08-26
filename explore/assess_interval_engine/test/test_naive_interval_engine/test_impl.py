import os

import pytest

from labw_utils.typing_importer import Type
from naive_interval_engine import IntervalEngineType
from naive_interval_engine.interval_tree_impl import IntervalTreeIntervalEngine
from naive_interval_engine.ncls_impl import NclsIntervalEngine
from naive_interval_engine.ne_impl import NumExprIntervalEngine
from naive_interval_engine.np_impl import NumpyIntervalEngine
from naive_interval_engine.pd_impl import PandasIntervalEngine

argvalues = (
    {"engine": PandasIntervalEngine},
    {"engine": NumpyIntervalEngine},
    {"engine": IntervalTreeIntervalEngine},
    {"engine": NumExprIntervalEngine},
    {"engine": NclsIntervalEngine},
)

test_filename = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_data.tsv")
test_match_filename = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_match.tsv")


@pytest.mark.parametrize(
    argnames="kwargs", argvalues=argvalues, ids=list(map(lambda d: d["engine"].__name__, argvalues))
)
def test_match(kwargs):
    engine_type: Type[IntervalEngineType] = kwargs["engine"]
    engine_test = engine_type(test_filename)
    engine_match = engine_type(test_match_filename)

    l_engine_match = list(engine_match)

    assert list(engine_test.match(l_engine_match[0])) == []
    assert list(engine_test.match(l_engine_match[1])) == [4, 6]
    assert list(engine_test.matches(engine_match)) == [[], [4, 6]]
    assert list(engine_match.matches(engine_test)) == [[0, 1], [], [], [], [], [], [], [], []]

    assert list(engine_test.match(("2", 0, 1))) == []
    assert list(engine_test.matches((("2", 0, 1), ("2", 13000, 14000)))) == [[], []]


@pytest.mark.parametrize(
    argnames="kwargs", argvalues=argvalues, ids=list(map(lambda d: d["engine"].__name__, argvalues))
)
def test_overlap(kwargs):
    engine_type: Type[IntervalEngineType] = kwargs["engine"]
    engine_test = engine_type(test_filename)
    engine_match = engine_type(test_match_filename)

    l_engine_match = list(engine_match)

    assert list(engine_test.overlap(l_engine_match[0])) == [0, 1, 2, 3]
    assert list(engine_test.overlap(l_engine_match[1])) == [0, 1, 2, 3, 4, 5, 6, 7, 8]
    assert list(engine_test.overlaps(engine_match)) == [[0, 1, 2, 3], [0, 1, 2, 3, 4, 5, 6, 7, 8]]
    assert list(engine_match.overlaps(engine_test)) == [[0, 1], [0, 1], [0, 1], [0, 1], [1], [1], [1], [1], [1]]

    assert list(engine_test.overlap(("2", 0, 1))) == []
    assert list(engine_test.overlaps((("2", 0, 1), ("2", 13000, 14000)))) == [[], []]


@pytest.mark.parametrize(
    argnames="kwargs", argvalues=argvalues, ids=list(map(lambda d: d["engine"].__name__, argvalues))
)
def test_io(kwargs):
    engine_type: Type[IntervalEngineType] = kwargs["engine"]
    engine_test = engine_type(test_filename)
    engine_match = engine_type(test_match_filename)

    l_engine_test = list(engine_test)
    l_engine_match = list(engine_match)
    assert l_engine_test[0] == ("1", 12006, 184926)
    assert len(l_engine_match) == 2
