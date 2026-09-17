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
#   --probe-switches      add run-time switches that remove single parts of the
#                         retarded dipole sector (section 105): CREM_PROBE_NO_POLES,
#                         _NO_MAGNETIZATION, _NO_LW_IN_MOMENT, _NO_MATERIAL,
#                         _NO_DIPOLE_ON_CHARGE, _NO_MOMENT_FORCE
#
# Usage: tools/build_variant.sh <tools/tool.cpp> <output binary> [options]
# The tool source is taken from the working tree; the model from HEAD.
set -euo pipefail
[ $# -ge 2 ] || { echo "usage: $0 <tools/tool.cpp> <output> [--max-jump X] [--freeze-spin-orbit]"; exit 2; }
tool=$1; out=$(realpath -m "$2"); shift 2
jump=""; freeze=0; revert=""; probe=0
while [ $# -gt 0 ]; do
  case $1 in
    --max-jump) [ $# -ge 2 ] || { echo "--max-jump needs a value"; exit 2; }
                jump=$2; shift 2;;
    --freeze-spin-orbit) freeze=1; shift;;
    --probe-switches) probe=1; shift;;
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
if [ $probe = 1 ]; then
  python3 - "$tmp/modules/electrodynamics.hpp" <<'PY2'
import sys
p=sys.argv[1]; s=open(p).read()
def rep(old,new):
    global s
    assert s.count(old)==1, "probe anchor not found exactly once: "+old[:50]
    s=s.replace(old,new)
rep("""    ElectromagneticField field{dual.magnetic*(-c*c),dual.electric};""",
    """    ElectromagneticField field{dual.magnetic*(-c*c),dual.electric};
    if(std::getenv("CREM_PROBE_NO_POLES")) field=ElectromagneticField{};""")
rep("""    if(const double floor=separationFloor(); floor>0.0) {
        double magnetizationTime=observationTime;""",
    """    if(const double floor=separationFloor(); floor>0.0&&!std::getenv("CREM_PROBE_NO_MAGNETIZATION")) {
        double magnetizationTime=observationTime;""")
rep("""        sourceCharge);
    const ElectromagneticField magneticDipole=retardedMagneticDipoleField(""",
    """        sourceCharge);
    if(std::getenv("CREM_PROBE_NO_LW_IN_MOMENT")) field=ElectromagneticField{};
    const ElectromagneticField magneticDipole=retardedMagneticDipoleField(""")
rep("""    const double materialRate=dipoleCouplingMaterialRate(
        state,history,targetIsFirst);""",
    """    if(std::getenv("CREM_PROBE_NO_MATERIAL")) return gradient/targetGamma;
    const double materialRate=dipoleCouplingMaterialRate(
        state,history,targetIsFirst);""")
rep("""            {secondField.electric+secondDipoleField.electric,
             secondField.magnetic+secondDipoleField.magnetic}),
        lorentzForce(secondCharge,s.secondVelocity,
            {firstField.electric+firstDipoleField.electric,
             firstField.magnetic+firstDipoleField.magnetic})};""",
    """            {secondField.electric+(std::getenv("CREM_PROBE_NO_DIPOLE_ON_CHARGE")?Vec3{}:secondDipoleField.electric),
             secondField.magnetic+(std::getenv("CREM_PROBE_NO_DIPOLE_ON_CHARGE")?Vec3{}:secondDipoleField.magnetic)}),
        lorentzForce(secondCharge,s.secondVelocity,
            {firstField.electric+(std::getenv("CREM_PROBE_NO_DIPOLE_ON_CHARGE")?Vec3{}:firstDipoleField.electric),
             firstField.magnetic+(std::getenv("CREM_PROBE_NO_DIPOLE_ON_CHARGE")?Vec3{}:firstDipoleField.magnetic)})};""")
rep("""    MutualForces tensorGradient=gDipoleForceEnabled
        ?MutualForces{covariantDipoleGradientForce(s,history,true),""",
    """    MutualForces tensorGradient=(gDipoleForceEnabled&&!std::getenv("CREM_PROBE_NO_MOMENT_FORCE"))
        ?MutualForces{covariantDipoleGradientForce(s,history,true),""")
open(p,'w').write(s)
PY2
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
echo "built $out (max-jump ${jump:-0.30}, freeze-spin-orbit $freeze, revert ${revert:-none}, probe-switches $probe)"
