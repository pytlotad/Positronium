#!/usr/bin/env python3
"""Cost-versus-gain split for a coarser photon quantum.

A larger quantum Q buys a longer wait for each photon (the hazard goes as
P/Q) and a bigger step down (each photon removes Q).  Because the emission
power runs away as P ~ a^-4 ~ |E|^4, the bigger step reaches the fast region
sooner.  Which of the two wins decides the SIGN of the collapse-time shift.

Model: circular orbit, |E| after a photon of energy Q is |E| + Q, waiting
time for that photon is Q/P(|E|), and P = A |E|^4.  A cancels in every ratio
below, so it is set to 1.
"""
import sys

def descend(startE, endE, quantum_of):
    """Total waiting time from binding startE to endE, in units of 1/A."""
    energy, total, steps = startE, 0.0, 0
    while energy < endE and steps < 100000:
        q = quantum_of(energy)
        if q <= 0.0:
            return None, steps
        total += q / energy**4          # wait = Q / P
        energy += q
        steps += 1
    return total, steps

def report(label, startE, endE, hbar_omega, ladder):
    t_hw, n_hw = descend(startE, endE, hbar_omega)
    t_lad, n_lad = descend(startE, endE, ladder)
    print(f"  {label:<34} hbar*omega {t_hw:9.4f} ({n_hw:3d} photons)"
          f"   ladder  {t_lad:9.4f} ({n_lad:3d})   ratio {t_lad/t_hw:6.3f}")

def main():
    # |E| in units of the Rydberg R, so level n means |E| = 1/n^2.
    # hbar*omega = 2R/n^3 = 2 |E|^{3/2};  ladder n -> n-1 at n = |E|^{-1/2}.
    hw = lambda E: 2.0 * E**1.5
    def ladder(E):
        n = E**-0.5
        if n >= 2.0:                   # n -> n-1
            lo = n - 1.0
            return 1.0/lo**2 - 1.0/n**2
        if n > 1.0:                    # only rung left is the ground state
            return 1.0 - 1.0/n**2
        return 2.0 * E**1.5            # below n=1 there is no rung at all

    print("Descent time in units of 1/A, |E| in rydbergs.  A ratio above 1")
    print("means the ladder is SLOWER than hbar*omega, below 1 faster.\n")
    report("n=2 -> n=1     (ladder 3x COARSER)", 0.25, 1.0, hw, ladder)
    report("n=1.05 -> n=0.5 (ladder 19x finer)", 1.0/1.05**2, 4.0, hw, ladder)
    report("n=1.02 -> n=0.5 (ladder 49x finer)", 1.0/1.02**2, 4.0, hw, ladder)
    report("n=1.01 -> n=0.5 (ladder 99x finer)", 1.0/1.01**2, 4.0, hw, ladder)
    print()
    print("One mechanism, both signs.  Which way it goes is set by whether the")
    print("substitution makes the quantum coarser or finer WHERE THE RUN SPENDS")
    print("ITS TIME -- near n=2 the ladder is 3x coarser and the collapse slows")
    print("(2.55 predicted here against 2.65 measured on the real code); packed")
    print("against n=1, where 6.6% of checkpoints sit, it is tens of times finer")
    print("and the collapse speeds up.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
