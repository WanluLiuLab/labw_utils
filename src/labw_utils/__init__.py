"""
labw_utils -- Utility Python functions & classes used in LabW

This is the top-level package of LabW Utils.
It also defines some commonly-used dependencies.

Import of this module may raise following errors & warnings:

- :py:obj:`RuntimeError`: If Python version lower than or equal to 3.6.
- :py:obj:`UserWarning`: If Python version is 3.7.
"""

from __future__ import annotations

__all__ = (
    "PackageSpec",
    "PackageSpecs",
    "UnmetDependenciesError",
    "__version__"
)

__version__ = "1.0.1"

import sys
import warnings

from typing import Optional

if sys.version_info <= (3, 6):
    raise RuntimeError("Python version <= 3.6, refuse to work.")
elif sys.version_info < (3, 8):
    warnings.warn("Python version == 3.7, not recommended.")

if sys.version_info >= (3, 9):
    from collections.abc import Iterable

    Dict = dict
else:
    from typing import Iterable, Dict


class PackageSpec:
    """
    Basic package specification. Can be used to specify packages in following package registry:

    - PYPI: https://pypi.org/
    - Conda: http://anaconda.org/

    .. versionadded:: 1.0.0
    """
    _name: str
    _conda_channel: Optional[str]
    _conda_name: Optional[str]
    _pypi_name: Optional[str]

    def __repr__(self):
        if self._conda_name is not None:
            if self._conda_channel is not None:
                conda_str = f"Use ``conda install -c {self._conda_channel} {self._conda_name}``"
            else:
                conda_str = f"Use ``conda install {self._conda_name}``"
        else:
            conda_str = ""
        if self._pypi_name is not None:
            pypi_str = f"Use ``pip install {self._pypi_name}``"
        else:
            pypi_str = ""
        return "; ".join((conda_str, pypi_str))

    def __init__(
            self,
            name: str,
            conda_channel: Optional[str],
            conda_name: Optional[str],
            pypi_name: Optional[str]
    ):
        self._name = name
        self._conda_name = conda_name
        self._pypi_name = pypi_name
        self._conda_channel = conda_channel

    @property
    def name(self) -> str:
        """
        Commonly-used name of that package.
        """
        return self._name

    @property
    def conda_name(self) -> Optional[str]:
        """
        Name as-is in Conda
        """
        return self._conda_name

    @property
    def conda_channel(self) -> Optional[str]:
        """
        Channel name as-is in Conda
        """
        return self._conda_channel

    @property
    def pypi_name(self) -> Optional[str]:
        """
        Name as-is in PyPI
        """
        return self._pypi_name


class PackageSpecs:
    """
    Package specifications.
    Maintains a list of :py:class:`PackageSpec`.
    Used in :py:class:`UnmetDependenciesError`.

    Current recognized optional dependencies:
    """
    _deps: Dict[str, PackageSpec] = {}

    @staticmethod
    def get(name: str) -> PackageSpec:
        """
        Get package specification.

        :param name: Name of package.
        :return: Specification of that package.
        :raises KeyError: If the package was not found.
        """
        return PackageSpecs._deps[name]

    @staticmethod
    def add(item: PackageSpec) -> None:
        """
        Add a package into the list.
        """
        PackageSpecs._deps[item.name] = item
        if PackageSpecs.__doc__ is None:  # Supress mypy
            PackageSpecs.__doc__ = ""
        PackageSpecs.__doc__ = PackageSpecs.__doc__ + f"\n    - ``{item.name}``: {item}"

    @staticmethod
    def iter_names() -> Iterable[str]:
        """
        Iterate known package names.
        """
        return iter(PackageSpecs._deps.keys())


PackageSpecs.add(PackageSpec(
    name="pandas",
    conda_name="pandas",
    pypi_name="pandas",
    conda_channel="conda-forge"
))
PackageSpecs.add(PackageSpec(
    name="numpy",
    conda_name="numpy",
    pypi_name="numpy",
    conda_channel="conda-forge"
))
PackageSpecs.add(PackageSpec(
    name="torch",
    conda_name="pytorch",
    pypi_name="torch",
    conda_channel="pytorch"
))
PackageSpecs.add(PackageSpec(
    name="pytables",
    conda_name="pytables",
    pypi_name="pytables",
    conda_channel="conda-forge"
))
PackageSpecs.add(PackageSpec(
    name="pysam",
    conda_name="pysam",
    pypi_name="pysam",
    conda_channel="bioconda"
))
PackageSpecs.add(PackageSpec(
    name="pyarrow",
    conda_name="pyarrow",
    pypi_name="pyarrow",
    conda_channel="conda-forge"
))
PackageSpecs.add(PackageSpec(
    name="fastparquet",
    conda_name="fastparquet",
    pypi_name="fastparquet",
    conda_channel="conda-forge"
))
PackageSpecs.add(PackageSpec(
    name="flask",
    conda_name="flask",
    pypi_name="flask",
    conda_channel="conda-forge"
))
PackageSpecs.add(PackageSpec(
    name="sqlalchemy",
    conda_name="sqlalchemy",
    pypi_name="sqlalchemy",
    conda_channel="conda-forge"
))
PackageSpecs.add(PackageSpec(
    name="psutil",
    conda_name="psutil",
    pypi_name="psutil",
    conda_channel="conda-forge"
))
PackageSpecs.add(PackageSpec(
    name="gevent",
    conda_name="gevent",
    pypi_name="gevent",
    conda_channel="conda-forge"
))
PackageSpecs.add(PackageSpec(
    name="tomli_w",
    conda_name="tomli-w",
    pypi_name="tomli-w",
    conda_channel="conda-forge"
))
PackageSpecs.add(PackageSpec(
    name="requests",
    conda_name="requests",
    pypi_name="requests",
    conda_channel="conda-forge"
))
PackageSpecs.add(PackageSpec(
    name="joblib",
    conda_name="joblib",
    pypi_name="joblib",
    conda_channel="conda-forge"
))
PackageSpecs.add(PackageSpec(
    name="jinja2",
    conda_name="jinja2",
    pypi_name="jinja2",
    conda_channel="conda-forge"
))
PackageSpecs.add(PackageSpec(
    name="matplotlib",
    conda_name="matplotlib",
    pypi_name="matplotlib",
    conda_channel="conda-forge"
))
PackageSpecs.add(PackageSpec(
    name="scipy",
    conda_name="scipy",
    pypi_name="scipy",
    conda_channel="conda-forge"
))
PackageSpecs.add(PackageSpec(
    name="snappy",
    conda_name="python-snappy",
    pypi_name="python-snappy",
    conda_channel="conda-forge"
))


class UnmetDependenciesError(RuntimeError):
    """
    An error indicating some additional packages should be installed.
    
    .. versionadded:: 1.0.0
    """
    _package_name: str

    def __init__(self, package_name: str):
        self._package_name = package_name
        err_message = f"{package_name} not installed; " + repr(PackageSpecs.get(package_name))
        super().__init__(err_message)

    @property
    def package_name(self) -> str:
        """Name of the missing package"""
        return self._package_name
