import importlib
import logging
import logging.handlers
import os
import shutil
import sys
import tempfile
import unittest

from labw_utils.stdlib.cpy310.pkgutil import resolve_name


class PkgutilTests(unittest.TestCase):
    def setUp(self):
        self.dirname = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.dirname)
        sys.path.insert(0, self.dirname)

    def tearDown(self):
        del sys.path[0]

    def test_name_resolution(self):
        success_cases = (
            ("os", os),
            ("os.path", os.path),
            ("os.path:pathsep", os.path.pathsep),
            ("logging", logging),
            ("logging:", logging),
            ("logging.handlers", logging.handlers),
            ("logging.handlers:", logging.handlers),
            ("logging.handlers:SysLogHandler", logging.handlers.SysLogHandler),
            ("logging.handlers.SysLogHandler", logging.handlers.SysLogHandler),
            ("logging.handlers:SysLogHandler.LOG_ALERT", logging.handlers.SysLogHandler.LOG_ALERT),
            ("logging.handlers.SysLogHandler.LOG_ALERT", logging.handlers.SysLogHandler.LOG_ALERT),
            ("builtins.int", int),
            ("builtins:int", int),
            ("builtins.int.from_bytes", int.from_bytes),
            ("builtins:int.from_bytes", int.from_bytes),
            ("builtins.ZeroDivisionError", ZeroDivisionError),
            ("builtins:ZeroDivisionError", ZeroDivisionError),
            ("os:path", os.path),
        )

        failure_cases = (
            (None, TypeError),
            (1, TypeError),
            (2.0, TypeError),
            (True, TypeError),
            ("", ValueError),
            ("?abc", ValueError),
            ("abc/foo", ValueError),
            ("foo", ImportError),
            ("os.foo", AttributeError),
            ("os.foo:", ImportError),
            ("os.pth:pathsep", ImportError),
            ("logging.handlers:NoSuchHandler", AttributeError),
            ("logging.handlers:SysLogHandler.NO_SUCH_VALUE", AttributeError),
            ("logging.handlers.SysLogHandler.NO_SUCH_VALUE", AttributeError),
            ("ZeroDivisionError", ImportError),
            ("os.path.9abc", ValueError),
            ("9abc", ValueError),
        )

        # add some Unicode package names to the mix.

        unicode_words = (
            "\u0935\u092e\u0938",
            "\xe9",
            "\xc8",
            "\uc548\ub155\ud558\uc138\uc694",
            "\u3055\u3088\u306a\u3089",
            "\u3042\u308a\u304c\u3068\u3046",
            "\u0425\u043e\u0440\u043e\u0448\u043e",
            "\u0441\u043f\u0430\u0441\u0438\u0431\u043e",
            "\u73b0\u4ee3\u6c49\u8bed\u5e38\u7528\u5b57\u8868",
        )
        for uw in unicode_words:
            d = os.path.join(self.dirname, uw)
            try:
                os.makedirs(d, exist_ok=True)
            except UnicodeEncodeError:
                # When filesystem encoding cannot encode uw: skip this test
                continue
            # make an empty __init__.py file
            f = os.path.join(d, "__init__.py")
            with open(f, "w") as f:
                f.write("")
                f.flush()
            # now import the package we just created; clearing the caches is
            # needed, otherwise the newly created package isn't found
            importlib.invalidate_caches()
            mod = importlib.import_module(uw)
            success_cases += ((uw, mod),)
            if len(uw) > 1:
                failure_cases += ((uw[:-1], ImportError),)

        # add an example with a Unicode digit at the start
        failure_cases += (("\u0966\u0935\u092e\u0938", ValueError),)

        for s, expected in success_cases:
            with self.subTest(s=s):
                o = resolve_name(s)
                self.assertEqual(o, expected)

        for s, exc in failure_cases:
            with self.subTest(s=s):
                with self.assertRaises(exc):
                    resolve_name(s)
