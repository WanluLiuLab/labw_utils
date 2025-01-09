from labw_utils import __version__
from labw_utils.commonutils import __doc__ as doc
from labw_utils.commonutils.libfrontend import setup_frontend

if __name__ == "__main__":
    setup_frontend(f"{__package__}._main", doc.splitlines()[1], __version__)
