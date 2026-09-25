# 197b/197f: exact rational spin algebra -- the matrix elements of
# S1.S2 and (S1-S2)_z in the total-spin basis, and the two divisions
# the model carries, integrated under its own default preparation.
#
# Run: python3 spinalgebra.py
from fractions import Fraction as F


def kron(A, B):
    n, m = len(A), len(B)
    return [[A[i // m][j // m] * B[i % m][j % m]
             for j in range(n * m)] for i in range(n * m)]


def zip4(A, B, f):
    return [[f(A[i][j], B[i][j]) for j in range(4)] for i in range(4)]


def mul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(4))
             for j in range(4)] for i in range(4)]


def scale(c, A):
    return [[c * A[i][j] for j in range(4)] for i in range(4)]


add = lambda a, b: a + b
sub = lambda a, b: a - b

h = F(1, 2)
I = [[F(1), F(0)], [F(0), F(1)]]
sz = [[h, F(0)], [F(0), -h]]
sp = [[F(0), F(1)], [F(0), F(0)]]
sm = [[F(0), F(0)], [F(1), F(0)]]
S1z, S2z = kron(sz, I), kron(I, sz)
S1p, S2p = kron(sp, I), kron(I, sp)
S1m, S2m = kron(sm, I), kron(I, sm)

flip = zip4(mul(S1p, S2m), mul(S1m, S2p), add)
S1S2 = zip4(mul(S1z, S2z), scale(F(1, 2), flip), add)
D = zip4(S1z, S2z, sub)

states = {'singlet': [F(0), F(1), F(-1), F(0)],
          't+': [F(1), F(0), F(0), F(0)],
          't0': [F(0), F(1), F(1), F(0)],
          't-': [F(0), F(0), F(0), F(1)]}
norm = {'singlet': F(2), 't+': F(1), 't0': F(2), 't-': F(1)}


def elem(M, a, b):
    va, vb = states[a], states[b]
    return sum(va[i] * sum(M[i][j] * vb[j] for j in range(4))
               for i in range(4))


print('<S1.S2> diagonal, units hbar^2:')
for s in states:
    print('  %-8s %s' % (s, elem(S1S2, s, s) / norm[s]))

print('\n(S1-S2)_z, units hbar:')
for a in states:
    row = []
    for b in states:
        v = elem(D, a, b)
        row.append('%-6s' % (v / F(2) if v else 0))
    print('  %-8s %s' % (a, '  '.join(row)))


def integ(lo, hi):
    def g(c):
        return F(1, 4) * (c + c * c / 2)
    return g(hi) - g(lo)


tot = integ(F(-1), F(1))
para = integ(F(1, 2), F(1))
print('\nP(2g) = %s   P(2g|para) = %s   P(2g|ortho) = %s'
      % (tot, para / F(1, 4), (tot - para) / F(3, 4)))
print('share of 2-gamma events from ortho-labelled draws = %s'
      % ((tot - para) / tot))
