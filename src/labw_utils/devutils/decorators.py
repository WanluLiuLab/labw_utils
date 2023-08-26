"""
``labw_utils.devutils.decorators`` -- Decorators for miscellaneous features.

.. versionadded:: 1.0.2
"""

__all__ = ("copy_doc", "chronolog", "create_class_init_doc_from_property", "supress_inherited_doc")

import inspect
import logging
import os
import types
import uuid

from labw_utils.stdlib.cpy310.pkgutil import resolve_name
from labw_utils.typing_importer import Any, List, Mapping, TypeVar, Callable

_InType = TypeVar("_InType")


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
    >>> AProxy.foo.__doc__
    'Woa'

    .. versionadded:: 1.0.2
    """
    if isinstance(source, str):
        source = resolve_name(source)

    def wrapper(func: Any) -> Callable:
        func.__doc__ = source.__doc__
        return func

    return wrapper


def supress_inherited_doc(force: bool = False, modify_overwritten: bool = False):
    """
    Stripping inhereited documentations from an object.

    >>> class A:
    ...     def f(self):
    ...         \"\"\"Woa\"\"\"
    ...         ...
    ...
    >>> @supress_inherited_doc(force=True)
    ... class B(A):
    ...     ...
    ...
    >>> class C(A):
    ...     ...
    ...
    >>> class D(A):
    ...     def f(self):
    ...         \"\"\"Woee\"\"\"
    ...         ...
    ...
    >>> @supress_inherited_doc(force=True, modify_overwritten=False)
    ... class E(A):
    ...     def f(self):
    ...         \"\"\"Woee\"\"\"
    ...         ...
    ...
    >>> @supress_inherited_doc(force=True, modify_overwritten=True)
    ... class F(A):
    ...     def f(self):
    ...         \"\"\"Woee\"\"\"
    ...         ...
    ...
    >>> print(A.f.__doc__)
    Woa
    >>> print(B.f.__doc__)
    Inherited from :py:mod:`labw_utils.devutils.decorators.A.f`
    >>> print(C.f.__doc__)
    Woa
    >>> print(D.f.__doc__)
    Woee
    >>> print(E.f.__doc__)
    Woee
    >>> print(F.f.__doc__)
    Inherited from :py:mod:`labw_utils.devutils.decorators.A.f`

    :param force: Proceed even if not in Sphinx build environment.
    :param modify_overwritten: Proceed even if this method was overewritten.
    :return: The object without stripped documentation.

    .. versionadded:: 1.0.0
    .. versionchanged:: 1.0.2
        The decorator becomes parameterized.
    .. warning::
        This method is well-tested. Use with care!
        Currently it works with methods and :py:class:`builtin,property` objects.
        Usage on other member types (e.g., constants of arbitrary type) not tested.
    """

    def empty_function(orig_func: Callable):
        def f(*args, **kwargs):
            return orig_func(*args, **kwargs)

        return f

    if os.getenv("SPHINX_BUILD") is not None or force:

        def real_decorator(obj: _InType) -> _InType:
            def _perform(require_del: bool):
                old_member = getattr(obj, member_name)
                if inspect.isfunction(old_member) or inspect.ismethod(old_member):
                    _new_attr = empty_function(old_member)
                elif isinstance(old_member, property):
                    _new_attr = property()  # FIXME: Function loss!
                else:
                    _new_attr = type(old_member.__class__.__name__, (object,), {})()  # FIXME: Add unit tests.
                if require_del:
                    delattr(obj, member_name)
                setattr(obj, member_name, _new_attr)
                setattr(getattr(obj, member_name), "__doc__", f"Inherited from :py:mod:`{mro_fullname}.{member_name}`")

            mro_name_member: Mapping[str : List[str]] = {}
            mro_name_object: Mapping[str:str] = {}
            for mro in obj.__mro__[1:]:  # Exclude myself
                mro_mod = inspect.getmodule(mro)
                if mro_mod is None:
                    raise TypeError
                mro_fullname = ".".join(
                    (mro_mod.__name__, getattr(mro, "__name__", ""))  # FIXME: bugs! Not fully qualified name!
                )
                mro_name_member[mro_fullname] = set()
                mro_name_object[mro_fullname] = mro
                for member_name in dir(mro):
                    if not member_name.startswith("_"):
                        mro_name_member[mro_fullname].add(member_name)
            obj_dir = filter(lambda x: not x.startswith("_"), dir(obj))
            for member_name in obj_dir:
                for mro_fullname, mro_member_names in mro_name_member.items():
                    mro = mro_name_object[mro_fullname]
                    if member_name in mro_member_names:
                        if getattr(obj, member_name) is not getattr(mro, member_name):
                            # Overwritten
                            if modify_overwritten:
                                _perform(require_del=True)
                        else:
                            _perform(require_del=False)
                        break
            return obj

    else:

        def real_decorator(obj: _InType) -> _InType:
            return obj

    return real_decorator


def doc_del_attr(
    attr_names: List[str],
    force: bool = False,
):
    """
    Delete attribute from object.

    :param force: Proceed even if not in Sphinx build environment.
    :param attr_names: Attributes to be removed.
    :return: The object with attributes deleted

    .. versionadded:: 1.0.3
    """
    if os.getenv("SPHINX_BUILD") is not None or force:

        def real_decorator(obj: _InType) -> _InType:
            for attr_name in attr_names:
                try:
                    delattr(obj, attr_name)
                except AttributeError:
                    pass  # TODO: Inherited not overwritten
            return obj

    else:

        def real_decorator(obj: _InType) -> _InType:
            return obj

    return real_decorator


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

    .. versionadded:: 1.0.2
    """

    def msg_decorator(f: types.FunctionType) -> Callable:
        if os.environ.get("SPHINX_BUILD") is not None:
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


def create_class_init_doc_from_property(
    text_before: str = "",
    text_after: str = "",
):
    """
    Place documentations at attributes to ``__init__`` function of a class.

    :param text_before: Text placed before parameters.
    :param text_after: Text placed after parameters.

    Example:

    >>> @create_class_init_doc_from_property()
    ... class TestInitDoc:
    ...     _a: int
    ...     _b: int
    ...     def __init__(self, a: int, b: int):
    ...         ...
    ...
    ...     @property
    ...     def a(self) -> int:
    ...         \"\"\"Some A value\"\"\"
    ...         return self._a
    ...
    ...     @property
    ...     def b(self) -> int:
    ...         \"\"\"Some B value\"\"\"
    ...         return self._b
    >>> print(TestInitDoc.__init__.__doc__)
    :param a: Some A value
    :param b: Some B value
    <BLANKLINE>

    Note that this example would NOT work:

    >>> @create_class_init_doc_from_property()
    ... class TestInitDoc:
    ...     a: int
    ...     \"\"\"Some A value\"\"\"
    ...
    ...     b: int
    ...     \"\"\"Some B value\"\"\"
    ...
    ...     def __init__(self, a: int, b: int):
    ...         ...
    >>> print(TestInitDoc.__init__.__doc__)
    <BLANKLINE>

    .. versionadded:: 1.0.2
    """

    def inner_dec(cls: _InType) -> _InType:
        init_func = cls.__init__
        mro = list(cls.__mro__)
        sig = inspect.signature(init_func)
        result_doc = ""
        for argname in sig.parameters.keys():
            curr_mro = list(mro)
            while curr_mro:
                curr_class = curr_mro.pop(0)
                try:
                    doc = getattr(curr_class, argname).__doc__
                except AttributeError:
                    continue
                if doc is None:
                    continue
                result_doc += f":param {argname}: {doc}\n"
                break

        init_func.__doc__ = text_before + result_doc + text_after
        return cls

    return inner_dec
