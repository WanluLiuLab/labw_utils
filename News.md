---
orphan: true
---

# NEWS

The initially released version is [`1.0.2`](https://pypi.org/project/labw-utils/1.0.2/) at 2023-07-17 CST.

## Version 1.X

- [`1.0.3`](https://pypi.org/project/labw-utils/1.0.3/) at 2023-12-29 CST:
  - Documentation enhancements for `labw_utils.PackageSpecs`.
  - A command-line interface added to `labw_utils.commonutils`.
  - Enum support added to `labw_utils.commonutils.stdlib_helper.argparse_helper`; Setting both `required` and `default` would generate `ValueError`.
  - Mainstream supported Python version shifted from 3.8 to 3.9; Python 3.7 excluded from automated tests with Python 3.11 and 3.12 added.
  - [`pytype`](https://google.github.io/pytype/) static linter removed.
  - Miscellaneous bug fixes.
- 1.0.4 (current):
  - `labw_utils` now comes with no dependencies; Use `labw_utils[defaults]` for default dependencies.
  - Support for reading [RepeatMasker](https://www.repeatmasker.org) output.
  - Support for parsing [RepBase](https://www.girinst.org/repbase) and [Dfam](https://www.dfam.org/home) database in [EMBL](https://raw.githubusercontent.com/enasequence/read_docs/master/submit/fileprep/flatfile_user_manual.txt) format.
  - Miscellaneous bug fixes.
