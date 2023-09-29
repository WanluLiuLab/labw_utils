#!/usr/bin/env bash
# shellcheck disable=SC2002
set -ue

git ls-files |
    grep \.py$ |
    while read -r fn; do
        cat "${fn}" |
            sed 's;\#.*;;' |
            sed 's;^[[:blank:]]*;;' |
            sed 's;[[:blank:]]*$;;' |
            sed 's;$;'" # ${fn}"';'
    done |
    grep 'import ' |
    grep -e '^from ' -e '^import' |
    grep -v '^\#' |
    grep -v 'typing' |
    grep -v '^>>>' |
    sort | uniq
