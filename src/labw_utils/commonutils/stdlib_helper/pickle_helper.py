"""
``labw_utils.stdlib_helper.pickle_helper`` -- Pickle helper with compression and progress-bar

We may firstly generate some random sequence

>>> import random
>>> import tempfile
>>> import os
>>> from labw_utils.commonutils.stdlib_helper.shutil_helper import rm_rf
>>> random_arr = [random.random() for _ in range(1000)]
>>> test_path = tempfile.mkdtemp()

Test with a normal file

>>> pickle_fn = os.path.join(test_path, 'rd.pickle')
>>> dump(random_arr, pickle_fn)
>>> unpickle_obj = load(pickle_fn)
>>> assert unpickle_obj == random_arr

This module can also handle compressed pickle

>>> pickle_fn = os.path.join(test_path, 'rd.pickle.xz')
>>> dump(random_arr, pickle_fn)
>>> unpickle_obj = load(pickle_fn)
>>> assert unpickle_obj == random_arr

Without a progress bar

>>> unpickle_obj = load(pickle_fn, with_tqdm=False)
>>> assert unpickle_obj == random_arr

Clean up the environment.

>>> rm_rf(test_path)

.. versionadded:: 1.0.2

.. todo::
    Support StringIO and friends
"""

__all__ = ("load", "dump")

import pickle
from pickle import Unpickler

from labw_utils.commonutils.lwio.safe_io import get_writer, get_reader
from labw_utils.commonutils.lwio.tqdm_reader import get_tqdm_reader
from labw_utils.typing_importer import Any


def load(filename: str, with_tqdm: bool = True) -> Any:
    """
    Unpickle a file with tqdm.

    :param filename: Filename to be load from
    :param with_tqdm: Whether to display a progress bar.
    :return: Picked object.

    .. versionadded:: 1.0.2
    """
    if with_tqdm:
        fd = get_tqdm_reader(filename, is_binary=True)
    else:
        fd = get_reader(filename, is_binary=True)
    up = Unpickler(fd)
    obj = up.load()
    fd.close()
    return obj


def dump(obj: Any, filename: str) -> None:
    """
    Pickle an object into a file.

    :param obj: The object to be pickled.
    :param filename: The filename to be written to.
        If the filename have compressed suffixes like ``xz``, it will be compressed.

    .. versionadded:: 1.0.2
    """
    with get_writer(filename, is_binary=True) as writer:
        pickle.dump(obj, writer)
