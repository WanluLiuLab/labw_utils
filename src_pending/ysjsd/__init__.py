from labw_utils import UnmetDependenciesError

_REQUIRED_MODNAMES = ("flask", "sqlalchemy", "psutil", "gevent", "tomli_w")

if os.environ.get("LABW_UTILS_UNDER_PYTEST", None) is not None:
    import pytest

    for modname in _REQUIRED_MODNAMES:
        _ = pytest.importorskip(modname)
else:
    pytest = None
    for modname in _REQUIRED_MODNAMES:
        try:
            __import__(modname)
        except ImportError:
            raise UnmetDependenciesError(modname)
