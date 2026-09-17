#!/bin/bash
# Build a tool against a temporary copy of HEAD with one of the variants the
# audit used, without touching the repository (audit sections 96-98, 100).
#
#   --max-jump X          maximumJumpParameter = X instead of 0.30 (a constexpr,
#                         so a checkpoint-size scan needs a separate build)
#   --freeze-spin-orbit   skip the secular spin-orbit half-step when
#                         CREM_FREEZE_SPIN_ORBIT is set at run time (section 100)
#   --revert COMMIT       reverse that commit's changes to modules/ on top of
#                         HEAD, leaving every later change in place (section 103)
#
# Usage: tools/build_variant.sh <tools/tool.cpp> <output binary> [options]
# The tool source is taken from the working tree; the model from HEAD.
set -euo pipefail
[ $# -ge 2 ] || { echo "usage: $0 <tools/tool.cpp> <output> [--max-jump X] [--freeze-spin-orbit]"; exit 2; }
tool=$1; out=$(realpath -m "$2"); shift 2
jump=""; freeze=0; revert=""
while [ $# -gt 0 ]; do
  case $1 in
    --max-jump) [ $# -ge 2 ] || { echo "--max-jump needs a value"; exit 2; }
                jump=$2; shift 2;;
    --freeze-spin-orbit) freeze=1; shift;;
    --revert) [ $# -ge 2 ] || { echo "--revert needs a commit"; exit 2; }
              revert=$2; shift 2;;
    *) echo "unknown option $1"; exit 2;;
  esac
done
root=$(git rev-parse --show-toplevel)
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
git -C "$root" archive HEAD | tar -x -C "$tmp"
mkdir -p "$tmp/$(dirname "$tool")"; cp "$root/$tool" "$tmp/$tool"
if [ -n "$revert" ]; then
  git -C "$root" diff "$revert^" "$revert" -- modules > "$tmp/revert.patch"
  (cd "$tmp" && git apply -R --whitespace=nowarn revert.patch) \
    || { echo "reverse patch of $revert does not apply to HEAD"; exit 1; }
fi
header="$tmp/modules/crem_collapse.hpp"
if [ -n "$jump" ]; then
  [ "$(grep -c 'constexpr double maximumJumpParameter=0.30;' "$header")" = 1 ] \
    || { echo "maximumJumpParameter anchor not found exactly once"; exit 1; }
  sed -i "s/constexpr double maximumJumpParameter=0.30;/constexpr double maximumJumpParameter=$jump;/" "$header"
fi
if [ $freeze = 1 ]; then
  python3 - "$header" <<'PY'
import sys
p=sys.argv[1]; s=open(p).read()
anchor="""        const auto advanceSpinOrbitHalf=[&](double semiMajorAxis,
                                            double elapsedTime) {
"""
assert s.count(anchor)==1, "advanceSpinOrbitHalf anchor not found exactly once"
s=s.replace(anchor,anchor+'            if(std::getenv("CREM_FREEZE_SPIN_ORBIT")) return true;\n')
open(p,'w').write(s)
PY
fi
(cd "$tmp" && g++ -std=c++20 -O2 -I . $(root-config --cflags) "$tool" -o "$out" $(root-config --libs))
echo "built $out (max-jump ${jump:-0.30}, freeze-spin-orbit $freeze, revert ${revert:-none})"
