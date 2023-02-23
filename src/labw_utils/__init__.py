__version__ = "0.1.14"


_deps = {
    "pandas": "",
    "pytables": "",
}

class UnmetDependenciesError(RuntimeError):
    ...  # TODO

    def __init__(self, package_name: str):
        ...
