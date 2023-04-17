from labw_utils import UnmetDependenciesError

_REQUIRED_MODNAMES = ("flask", "sqlalchemy", "psutil", "gevent", "tomli_w")

try:
    import pytest

    for modname in _REQUIRED_MODNAMES:
        _ = pytest.importorskip(modname)
except ImportError:
    pytest = None
    for modname in _REQUIRED_MODNAMES:
        try:
            __import__(modname)
        except ImportError:
            raise UnmetDependenciesError(modname)
