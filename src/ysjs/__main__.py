from libysjs import __version__
from labw_utils.commonutils.libfrontend import setup_frontend

if __name__ == '__main__':
    setup_frontend(
        "ysjs._main",
        "ysjs -- Commandline Interface of YSJS",
        __version__
    )
