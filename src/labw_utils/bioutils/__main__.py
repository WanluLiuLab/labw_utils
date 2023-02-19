from labw_utils import __version__
from labw_utils.commonutils.libfrontend import setup_frontend

if __name__ == '__main__':
    setup_frontend(
        "labw_utils.bioutils._main",
        "labw_utils.bioutils -- Biological Utilities used in LabW projects",
        __version__
    )
