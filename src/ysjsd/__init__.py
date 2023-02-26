import os

from labw_utils import UnmetDependenciesError

try:
    import flask
except ImportError:
    raise UnmetDependenciesError("flask")

try:
    import sqlalchemy
except ImportError:
    raise UnmetDependenciesError("sqlalchemy")

try:
    import psutil
except ImportError:
    raise UnmetDependenciesError("psutil")

try:
    import gevent
except ImportError:
    raise UnmetDependenciesError("gevent")

try:
    import tomli
except ImportError:
    raise UnmetDependenciesError("tomli")

try:
    import tomli_w
except ImportError:
    raise UnmetDependenciesError("tomli_w")
APP_DIR = os.path.dirname(os.path.abspath(__file__))
APP_NAME = "YSJSD BACKEND"
NOT_READY = ("Not ready\n", 500)
