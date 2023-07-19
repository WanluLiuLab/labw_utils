from labw_utils.commonutils.libfrontend import setup_frontend
from libysjs import __version__

if __name__ == '__main__':
    setup_frontend(
        "ysjs._main",
        "ysjs -- Commandline Interface of YSJS",
        __version__
    )
