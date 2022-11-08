#!/usr/bin/env bash
printf "chr\ts\te\n" > "${2}"
pv "${1}" | grep -v '^#' | cut -f 1,3-5 | grep transcript | cut -f 1,3,4 >> "${2}"
