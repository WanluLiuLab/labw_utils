"""
l``abw_utils.stdlib_helper.sys_helper`` -- A rewritten of various sysadmin commands.

.. versionadded:: 1.0.2
"""

__all__ = ("is_user_admin",)

import ctypes
import os

if os.name == "nt":

    def is_user_admin() -> int:
        try:
            return ctypes.windll.shell32.IsUserAnAdmin()  # type: ignore
        except AttributeError:
            return -1

elif os.name == "posix":

    def is_user_admin() -> int:
        if os.getuid() == 0:
            return 1
        else:
            return 0

else:

    def is_user_admin() -> int:
        return -1


is_user_admin.__doc__ = """
To detect whether the user is root (admin) or not.

A part of this function comes from
<https://stackoverflow.com/questions/19672352/how-to-run-script-with-elevated-privilege-on-windows>

.. warning ::
    For Microsoft Windows users: requires Windows XP SP2 or higher!

:return: 0=no, 1=yes, -1=error

.. versionadded:: 1.0.2
"""
