#!/usr/bin/env python3
"""Aggregate CREM_TILT probe output into the ensemble statistic that decides
whether the para/ortho dipole-dipole sign survives averaging.

The probe (modules/crem_collapse.hpp, gated on CREM_TILT) prints one line per
secular spin-orbit half-step:

    CREM_TILT <time_s> <a_m> <cos_after_1> <cos_after_2> <cos_before_1>

cos is the cosine of the angle between a particle's magnetic moment and the
orbital angular momentum -- the tilt.  "before" is the state the half-step was
handed, "after" is what it left behind, so the pair brackets one application of
the coupled transport and makes the comparison paired.

Why the tilt is the quantity.  The orbit-averaged RADIAL dipole-dipole force
carries the second Legendre polynomial of the tilt,

    F_r  =  +2.06e-13 N * P2(cos tilt)     for para (aligned moments)
    F_r  =  -2.06e-13 N * P2(cos tilt)     for ortho (anti-aligned)

at the pair Bohr radius, so P2 = +1 (moments along the orbit normal) is
repulsive for para, P2 = -1/2 (moments in the orbital plane) is attractive for
it, and the sign crosses zero at the magic angle, 54.7 degrees.  "Para attracts
and ortho repels" is therefore a statement about the tilt distribution, not
about the spin state.  Isotropic tilts give <P2> = 0 exactly, with per-sample
sd(P2) = sqrt(1/5) and hence a standard error of 0.447/sqrt(N).

Usage:
    CREM_TILT=1 ./positronium --mode statistical --phenomenon 1 \
        --runs 200 --seed 101 --crem-wallclock-budget-s 8 > para.log
    CREM_TILT=1 ./positronium --mode statistical --phenomenon 2 \
        --runs 200 --seed 101 --crem-wallclock-budget-s 8 > ortho.log
    python3 tools/tilt_ensemble.py para.log ortho.log
"""
import math
import sys


def read(path):
    """Rows of (time, a, cosAfter1, cosAfter2, cosBefore1)."""
    rows = []
    with open(path) as handle:
        for line in handle:
            if not line.startswith("CREM_TILT "):
                continue
            parts = line.split()
            if len(parts) != 6:
                continue
            try:
                rows.append(tuple(float(value) for value in parts[1:]))
            except ValueError:
                continue
    return rows


def legendre2(cosine):
    return 0.5 * (3.0 * cosine * cosine - 1.0)


def moments(values):
    """Mean and standard error of the mean."""
    count = len(values)
    if count == 0:
        return float("nan"), float("nan"), 0
    mean = sum(values) / count
    if count == 1:
        return mean, float("nan"), 1
    variance = sum((v - mean) ** 2 for v in values) / (count - 1)
    return mean, math.sqrt(variance / count), count


def verdict(mean, error):
    if error != error or error == 0.0:
        return "exactly unchanged" if mean == 0.0 else "no error estimate"
    sigmas = abs(mean) / error
    return "consistent with zero" if sigmas < 3.0 else f"NONZERO ({sigmas:.1f} sigma)"


def line(label, values):
    mean, error, count = moments(values)
    print(f"  {label:<34} {mean:>+9.4f} +/- {error:<8.4f} N={count:<5} "
          f"{verdict(mean, error)}")


def summarize(label, rows):
    if not rows:
        print(f"\n=== {label} ===\n  no CREM_TILT samples -- was the run started "
              f"with CREM_TILT=1?")
        return
    # The first half-step of each trajectory carries the tilt exactly as the
    # sampler drew it, which is the isotropic null this is tested against.
    first = [row for row in rows if row[0] == 0.0]
    print(f"\n=== {label} ===")
    print(f"  trajectories {len(first)}, samples {len(rows)}, "
          f"time span 0 .. {max(r[0] for r in rows):.4g} s")
    line("P2 as sampled, before transport", [legendre2(r[4]) for r in first])
    line("P2 after one half-step", [legendre2(r[2]) for r in first])
    line("P2 over every sample", [legendre2(r[2]) for r in rows])
    line("d|cos tilt| per half-step, paired",
         [abs(r[2]) - abs(r[4]) for r in first])
    print(f"  <cos^2 tilt> over every sample      "
          f"{sum(r[2]**2 for r in rows)/len(rows):.5f}   (isotropic 0.33333)")

    span = max(r[0] for r in rows)
    if span > 0.0:
        print("\n  P2 against time (late windows are few surviving trajectories,"
              "\n  so read them as survivorship, not as a trend):")
        for index in range(4):
            lo, hi = span * index / 4, span * (index + 1) / 4
            window = [legendre2(r[2]) for r in rows
                      if lo <= r[0] < hi or (index == 3 and r[0] == hi)]
            if window:
                print(f"    {lo*1e12:>8.4g}..{hi*1e12:<9.4g} ps  "
                      f"{sum(window)/len(window):>+8.4f}  ({len(window)} samples)")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    labels = ["PARA (aligned moments)", "ORTHO (anti-aligned moments)"]
    for index, path in enumerate(sys.argv[1:]):
        summarize(labels[index] if index < len(labels) else path, read(path))
    print("\nSamples inside one trajectory are correlated; the per-trajectory\n"
          "rows above (everything but 'over every sample') avoid that.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
