#!/bin/bash
cd "$(dirname "$0")"
for s in 0.05 0.2417 0.4834 0.9668 1.9336 3.8672 9.668 452; do
  CREM_MAGNETIC_RADIUS_SCALE=$s CREM_SKIP_CENSUS=1 \
    /tmp/radius 300 42 > radius_$s.txt 2>&1
done
run() { s=$1
  CREM_MAGNETIC_RADIUS_SCALE=$s /tmp/engdrift 1.0 0 0 1 1e-8 256 \
    > orbit_$s.txt 2>&1
}
export -f run
printf '%s\n' 0.05 0.2417 0.4834 0.9668 1.9336 3.8672 9.668 452 \
  | xargs -P 4 -I{} bash -c 'run {}'
echo SCAN139_DONE
