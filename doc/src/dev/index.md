# Developing `labw_utils`

## Environment Variables

- `LABW_UTILS_SPHINX_BUILD`: Should be set while building Sphinx documents. This environment variable may supress Chronolog and change `tqdm` implementation to the silent one.
- `LABW_UTILS_UNDER_PYTEST`: Should be set while executing automatic tests. This environment variable will change the behavior of importing third-party packages. If an import is failed, it would raise `pytest.Skipped` instead of `ImportError` or `UnmetDependenciesError`.
