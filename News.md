---
orphan: true
---

# NEWS

The initially released version is [`1.0.2`](https://pypi.org/project/labw-utils/1.0.2/) at 2023-07-17 CST.

## Version 1.X

- `1.0.3`:
  - Documentation enhancements for `labw_utils.PackageSpecs`.
  - A command-line interface added to `labw_utils.commonutils`
  - Enum support added to `labw_utils.commonutils.stdlib_helper.argparse_helper`; Setting both `required` and `default` would generate `ValueError`.
  - Mainstream supported Python version shifted from 3.8 to 3.9; Python 3.7 excluded from automated tests with Python 3.11 and 3.12 added.
  - [`pytype`](https://google.github.io/pytype/) static linter removed.
  - Miscellaneous bug fixes.
