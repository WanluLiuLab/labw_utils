from labw_utils import __version__
from labw_utils.commonutils.libfrontend import setup_frontend

if __name__ == '__main__':
    setup_frontend(
        "test_labw_utils.commonutils.libfrontend_test_main",
        "",
        __version__
    )
