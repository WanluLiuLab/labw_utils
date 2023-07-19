"""
``labw_utils.stdlib.cpy311`` -- Code copied from CPython 3.11 standard library for compatibility.

This mainly contains following modules:

- ``tomllib`` -- TOML Parser.
  Please note that this can only parse TOML files. It cannot be used as a serializer.
  Signature of :py:func:`labw_utils.stdlib.cpy311.tomllib.load` and
  :py:func:`labw_utils.stdlib.cpy311.tomllib.loads` was modified for Python 3.6/3.7 compatibility.
  
  .. versionadded:: 1.0.2
"""
