from labw_utils.commonutils.libfrontend import setup_frontend
from labw_utils import __version__

if __name__ == '__main__':
    setup_frontend(
        "labw_utils.bioutils.main",
        "DESC TODO",
        __version__
    )
