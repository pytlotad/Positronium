#!/usr/bin/env bash

set -euo pipefail

if (( $# != 1 )); then
    echo "usage: $0 PATH_TO_SANITIZED_BINARY" >&2
    exit 64
fi

binary=$1
if [[ $binary != /* ]]; then
    binary="$PWD/${binary#./}"
fi
if [[ ! -x $binary ]]; then
    echo "sanitized production binary is not executable: $binary" >&2
    exit 66
fi

smoke_dir=$(mktemp -d "${TMPDIR:-/tmp}/positronium-sanitizer.XXXXXX")
cleanup() {
    status=$?
    trap - EXIT
    rm -rf -- "$smoke_dir"
    exit "$status"
}
trap cleanup EXIT

cd "$smoke_dir"

# CREM collapse reaches the large production-only header and its censored
# statistical report.  A short budget is intentional: sanitizer CI verifies
# memory/UB safety, not the publication-quality completion fraction.
set +e
"$binary" --mode statistical --phenomenon 1 --runs 1 --seed 42 \
    --level 1 --crem-wallclock-budget-s 5 >crem-smoke.log 2>&1
crem_status=$?
set -e
if (( crem_status != 0 )); then
    if (( crem_status != 2 )) \
        || ! grep -q "No CREM trajectory reached the collision boundary" \
            crem-smoke.log; then
        echo "CREM sanitizer smoke failed with status $crem_status" >&2
        sed -n '1,240p' crem-smoke.log >&2
        exit "$crem_status"
    fi
    echo "CREM sanitizer smoke: expected right-censoring after the short budget"
fi

# Experiment 5 exercises collision classification, censoring regressions,
# statistical analysis and the multi-page PDF export path.
"$binary" --mode statistical --phenomenon 5 --runs 1 --seed 42 \
    --interaction-energy-ev 0.6

first_pdf=$(find distributions -type f -name '*.pdf' -print -quit)
if [[ -z $first_pdf ]]; then
    echo "production sanitizer smoke did not export any PDF" >&2
    exit 1
fi
