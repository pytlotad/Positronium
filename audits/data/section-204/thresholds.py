# 204c/204e: the disagreement between the para/ortho LABEL
# (dipoleAlignment >= t) and the drawn MULTIPLICITY (Bernoulli against
# (1+cos)/2), as a function of the threshold, under the model's default
# uniform-cosine preparation.  Exact integrals, no simulation.
#
# Run: python3 thresholds.py
from fractions import Fraction as F


def disagreement(t):
    """P(label para and draw 3g) + P(label ortho and draw 2g)."""
    def para_side(c):
        return F(1, 4) * (c - c * c / 2)

    def ortho_side(c):
        return F(1, 4) * (c + c * c / 2)

    return ((para_side(F(1)) - para_side(t))
            + (ortho_side(t) - ortho_side(F(-1))))


def label_share(t):
    """P(cos >= t) for a uniform cosine."""
    return (F(1) - t) / 2


print('  %10s %12s %12s %14s' % ('threshold', 'disagree', 'exact',
                                 'para label'))
for t in (F(-1), F(-1, 2), F(0), F(1, 4), F(1, 2), F(3, 4), F(1)):
    d = disagreement(t)
    print('  %10s %12.6f %12s %14s' % (t, float(d), d, label_share(t)))

print('\nd/dt of the disagreement is t/2, so the minimum sits at t = 0')
print('and equals 1/4.  No threshold does better under this preparation.')
print('\nreconciliation with audit 197f:')
print('  P(2g and ortho label) = %s' % (F(1, 2) * F(9, 16)))
print('  P(3g and para label)  = %s' % (F(1, 4) - F(1, 2) * F(7, 16)))
print('  sum                   = %s = disagreement at t = 1/2'
      % (F(1, 2) * F(9, 16) + F(1, 4) - F(1, 2) * F(7, 16)))
