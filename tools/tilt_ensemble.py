#!/usr/bin/env python3
"""Aggregate CREM_TILT probe output into the ensemble average that decides
whether the para/ortho dipole-dipole sign survives averaging.

The probe (modules/crem_collapse.hpp, gated on CREM_TILT) prints one line per
secular spin-orbit half-step:

    CREM_TILT <time_s> <semiMajorAxis_m> <cos_tilt_1> <cos_tilt_2>

cos_tilt is the cosine of the angle between each particle's magnetic moment and
the orbital angular momentum.  The orbit-averaged RADIAL dipole-dipole force
carries the second Legendre polynomial P2(cos tilt) = (3cos^2-1)/2, so the
ensemble mean <P2> is the quantity of interest:

    <P2> = 0  -> the attraction/repulsion cancels over the ensemble, and
                 "para attracts, ortho repels" is a statement about one
                 geometry rather than about the state.
    <P2> != 0 -> a systematic sign survives, and the channels really do differ
                 in the sign of their dipole-dipole shift.

Isotropic tilts give <P2> = 0 exactly, with per-sample sd(P2) = sqrt(1/5),
hence a standard error of 0.447/sqrt(N).

Usage:
    CREM_TILT=1 ./positronium --mode statistical --phenomenon 1 \
        --runs 40 --seed 7 --crem-wallclock-budget-s 25 > para.log
    CREM_TILT=1 ./positronium --mode statistical --phenomenon 2 \
        --runs 40 --seed 7 --crem-wallclock-budget-s 25 > ortho.log
    python3 tools/tilt_ensemble.py para.log ortho.log
"""
import math
import sys


def read(path):
    """Rows of (time, semiMajorAxis, cosTilt1, cosTilt2) from a run log."""
    rows = []
    with open(path) as handle:
        for line in handle:
            if not line.startswith("CREM_TILT "):
                continue
            parts = line.split()
            if len(parts) != 5:
                continue
            try:
                rows.append(tuple(float(value) for value in parts[1:]))
            except ValueError:
                continue
    return rows


def split_trajectories(rows):
    """The probe restarts each trajectory at t = 0."""
    runs, current = [], []
    for row in rows:
        if row[0] == 0.0 and current:
            runs.append(current)
            current = []
        current.append(row)
    if current:
        runs.append(current)
    return runs


def legendre2(cosine):
    return 0.5 * (3.0 * cosine * cosine - 1.0)


def summarize(label, rows):
    if not rows:
        print(f"{label}: no CREM_TILT samples -- was the run built with the "
              f"probe and started with CREM_TILT=1?")
        return
    runs = split_trajectories(rows)
    values = [legendre2(row[2]) for row in rows]
    count = len(values)
    mean = sum(values) / count
    # Samples inside one trajectory are correlated, so quote the error two ways:
    # naive per-sample, and per-trajectory, which is the honest one.
    naive = 0.4472135955 / math.sqrt(count)
    per_run = [sum(legendre2(r[2]) for r in run) / len(run) for run in runs]
    if len(per_run) > 1:
        run_mean = sum(per_run) / len(per_run)
        variance = sum((v - run_mean) ** 2 for v in per_run) / (len(per_run) - 1)
        run_error = math.sqrt(variance / len(per_run))
    else:
        run_mean, run_error = per_run[0], float("nan")

    print(f"\n=== {label} ===")
    print(f"trajectories                 {len(runs)}")
    print(f"samples                      {count}")
    print(f"time span                    0 .. {max(r[0] for r in rows):.4g} s")
    print(f"<cos^2 tilt>                 {sum(r[2]**2 for r in rows)/count:.5f}"
          f"   (isotropic 0.33333)")
    print(f"<P2> per sample              {mean:+.5f} +/- {naive:.5f} (naive)")
    print(f"<P2> per trajectory          {run_mean:+.5f} +/- {run_error:.5f}"
          f"   <- use this one")
    if run_error == run_error and run_error > 0.0:
        sigmas = abs(run_mean) / run_error
        print(f"distance from zero           {sigmas:.2f} sigma -> "
              f"{'CONSISTENT WITH ZERO' if sigmas < 3.0 else 'NONZERO'}")

    print("\ncos(tilt) histogram (isotropic is flat at 10% per bin):")
    bins = [0] * 10
    for row in rows:
        bins[min(9, int((row[2] + 1.0) / 2.0 * 10))] += 1
    for index, hits in enumerate(bins):
        share = 100.0 * hits / count
        print(f"  [{-1.0+0.2*index:+.1f},{-1.0+0.2*(index+1):+.1f})"
              f"  {share:5.1f}%  {'#' * int(round(share * 1.5))}")

    print("\n<P2> against time, to show whether the run drives it away from 0:")
    span = max(r[0] for r in rows)
    if span > 0.0:
        print(f"  {'t window [ps]':<22} {'<P2>':<12} samples")
        for bin_index in range(5):
            lo, hi = span * bin_index / 5, span * (bin_index + 1) / 5
            window = [legendre2(r[2]) for r in rows if lo <= r[0] < hi
                      or (bin_index == 4 and r[0] == hi)]
            if not window:
                continue
            print(f"  {lo*1e12:>9.4g}..{hi*1e12:<10.4g} "
                  f"{sum(window)/len(window):<+12.5f} {len(window)}")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    labels = ["PARA (aligned moments)", "ORTHO (anti-aligned moments)"]
    for index, path in enumerate(sys.argv[1:]):
        summarize(labels[index] if index < len(labels) else path, read(path))
    print("\nReminder: samples within one trajectory are correlated, which is "
          "why\nthe per-trajectory error is the one to quote.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
