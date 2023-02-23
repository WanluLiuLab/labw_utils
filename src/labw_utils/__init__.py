__version__ = "0.1.14"

_deps = {
    "pandas": ("pandas", "pandas"),
    "pytables": ("pytables", "pytables"),
    "pysam": ("pysam", "pysam"),
    "fastparquet": ("fastparquet", "fastparquet"),
    "flask": ("flask", "flask"),
    "sqlalchemy": ("sqlalchemy", "sqlalchemy"),
    "psutil": ("psutil", "psutil"),
    "gevent": ("gevent", "gevent"),
    "tomli": ("tomli", "tomli"),
    "tomli_w": ("tomli-w", "tomli-w")
}


class UnmetDependenciesError(RuntimeError):
    err_message: str

    def __init__(self, package_name: str):
        err_message = f"{package_name} not installed; " \
                      f"Use conda install -c conda-forge -c bioconda {_deps[package_name][0]}, " \
                      f"or pip install {_deps[package_name][1]}."
        super().__init__(err_message)
        self.err_message = err_message
