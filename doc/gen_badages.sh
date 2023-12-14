#!/usr/bin/env bash

python -m pybadges \
    --left-text "Build" \
    --right-text "Failed" \
    --right-color "red" \
    > src/installation/build_failed.svg
python -m pybadges \
    --left-text "Build" \
    --right-text "Passed" \
    --right-color "greem" \
    > src/installation/build_passed.svg
python -m pybadges \
    --left-text "Conda" \
    --right-text "Unmet Dependencies" \
    --right-color "red" \
    > src/installation/unmet_dependencies.svg
python -m pybadges \
    --left-text "Build" \
    --right-text "Ignored" \
    --right-color "orange" \
    > src/installation/build_ignored.svg
