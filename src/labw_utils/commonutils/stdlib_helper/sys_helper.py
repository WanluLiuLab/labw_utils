"""
sys_helper.py -- A rewritten of various sysadmin commands.

Due to the issues caused by default :py:mod:`sys` and :py:mod:`platform`,
we rewrote some functions inside.

This module is intended to work on Microsoft Windows.
"""

import ctypes
import os


def is_user_admin() -> int:
    """
    To detect whether the user is root (admin) or not.
    A part of this function comes from
    <https://stackoverflow.com/questions/19672352/how-to-run-script-with-elevated-privilege-on-windows>

    :return: 0=no, 1=yes, -1=error
    """
    if os.name == 'nt':
        # WARNING: requires Windows XP SP2 or higher!
        try:
            return ctypes.windll.shell32.IsUserAnAdmin()
        except AttributeError:
            return -1
    elif os.name == 'posix':
        if os.getuid() == 0:  # root
            return 1
        else:
            return 0
    else:
        return -1