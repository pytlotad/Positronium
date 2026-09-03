#!/bin/bash
# Compile one modules/*.hpp in isolation; with no argument, sweeps all 36.
#
# Three headers exist only in the validation configuration: positronium.cpp
# includes them inside #ifdef POSITRONIUM_ENABLE_FIELD_VALIDATION, and their
# bodies are either wholly inside that guard or use types electrodynamics.hpp
# defines only there.  They are compiled here with that macro, which is what
# "compiles on its own" means for a conditional module -- without it the probe
# would compile an empty file and report a PASS that tested nothing.
cd /home/teddy/Projekty/Positronium
SP="$(mktemp -d)"
trap 'rm -rf "$SP"' EXIT
one() {
    local extra=""
    case "$1" in
        *maxwell_validation.hpp|*maxwell_validation_backend.hpp|*legendre_anisotropy_fit.hpp)
            extra="-DPOSITRONIUM_ENABLE_FIELD_VALIDATION" ;;
    esac
    echo "#include \"$1\"" > "$SP/probe.cpp"
    g++ -fsyntax-only -std=c++20 -I. $(root-config --cflags) $extra \
        "$SP/probe.cpp" 2> "$SP/probe_err.txt"
}
if [ -n "$1" ]; then
    if one "$1"; then echo "PASS $1"; else
        echo "FAIL $1 ($(grep -c 'error:' "$SP/probe_err.txt") errors)"
        grep 'error:' "$SP/probe_err.txt" | head -"${2:-12}"
    fi
else
    pass=0; fail=0
    for h in modules/*.hpp; do
        if one "$h"; then pass=$((pass+1)); else
            fail=$((fail+1)); echo "FAIL $h ($(grep -c 'error:' "$SP/probe_err.txt") errors)"
        fi
    done
    echo "--- self-contained: $pass/36, failing: $fail ---"
fi
