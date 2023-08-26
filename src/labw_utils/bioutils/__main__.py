from labw_utils import __version__
from labw_utils.commonutils.libfrontend import setup_frontend
from labw_utils.bioutils import __doc__ as doc

if __name__ == "__main__":
    setup_frontend(f"{__package__}._main", doc.splitlines()[1], __version__)
