import logging
import os
import types
import uuid
from typing import Callable, Any


def copy_doc(source: Any) -> Callable:
    """
    The following piece of code is from
    https://stackoverflow.com/questions/68901049/copying-the-docstring-of-function-onto-another-function-by-name
    by Iced Chai at Aug 24, 2021 at 2:56

    This wrapper copies docstring from one function to another.

    Use Example: copy_doc(self.copy_func)(self.func) or used as deco

    >>> class Test:
    ...     def foo(self) -> None:
    ...         \"\"\"Woa\"\"\"
    ...         ...
    ...
    ...     @copy_doc(foo)
    ...     def this(self) -> None:
    ...         pass
    >>> Test.this.__doc__
    'Woa'

    This function should be used on so-called "proxy" classes. For example,

    >>> class A:
    ...     def foo(self) -> None:
    ...         \"\"\"Woa\"\"\"
    ...         ...
    ...
    >>> class AProxy:
    ...     _A: A
    ...     @copy_doc(A.foo)
    ...     def foo(self) -> None:
    ...         self._A.foo()
    """

    def wrapper(func: Any) -> Callable:
        func.__doc__ = source.__doc__
        return func

    return wrapper


def chronolog(display_time: bool = False, log_error: bool = False):
    """
    The logging decorator, will inject a logger variable named _lh to the code.
    From <https://stackoverflow.com/questions/17862185/how-to-inject-variable-into-scope-with-a-decorator>

    .. note::
        The :py:func:`error` (or :py:func:`exception`, :py:func:`critical`, :py:func:`fatal`
        functions DO NOT exit the program! You have to exit the program by yourself!

    .. warning::
        Call this function, do NOT call functions inside this function!

    :param display_time: Whether to display calling time, arguments and return value in log level.
    :param log_error: Whether add error captured
    """

    def msg_decorator(f: Callable) -> Callable:
        f: types.FunctionType
        if os.environ.get('SPHINX_BUILD') is not None:
            return f  # To make Sphinx get the right result.

        def inner_dec(*args, **kwargs):
            """
            Decorator which performs the logging and do the work.

            :param args: Unnamed arguments of the decorated function call.
            :param kwargs: Named arguments of the decorated function call.
            :return: The return value of the decorated function call.
            :raise: The return value of the decorated function call.
            """
            call_id = f"CHRONOLOG CID={uuid.uuid4()}"
            try:
                _ = f.__globals__
            except AttributeError:
                return f(*args, **kwargs)
            lh = logging.getLogger(f.__module__)
            if display_time:
                args_repr = [repr(a) for a in args]
                kwargs_repr = [f"{k}={v!r}" for k, v in kwargs.items()]
                signature = ", ".join(args_repr + kwargs_repr)
                # FIXME: Function name incorrect!
                lh.debug("%s %s(%s)", call_id, f.__name__, signature)
            res = None
            try:
                res = f(*args, **kwargs)
            except Exception as e:
                if log_error:
                    lh.exception("%s exception inside func: %s", call_id, str(e), stack_info=True, exc_info=True)
                raise e
            finally:
                lh.debug("%s returns %s", f.__name__, repr(res))
            return res

        return inner_dec

    return msg_decorator
