# 200b/200c: the Kepler harmonic amplitudes in closed form, checked first
# against a pure Kepler ellipse (which caught an error in the first form I
# wrote) and then against the model's measured spectrum.
#
# Run: python3 besselcheck.py
import math
import cmath


def J(n, x, terms=80):
    return sum(((-1) ** m) / (math.factorial(m) * math.factorial(m + n))
               * (x / 2.0) ** (2 * m + n) for m in range(terms))


def Jp(n, x):
    return (J(n - 1, x) - J(n + 1, x)) / 2.0


def amp(n, e):
    # sin(xi) = (2/e) sum (1/n) J_n(ne) sin(nM) follows from Kepler's
    # equation directly, so y_n/a = (2 sqrt(1-e^2)/(e n)) J_n(ne) and the
    # weight below is (1-e^2)/e^2 -- NOT its square.
    w = (1.0 - e * e) / (e * e)
    return 2.0 * n * math.sqrt(Jp(n, n * e) ** 2 + w * J(n, n * e) ** 2)


def kepler_dft(e, N=8192, K=7):
    xs, ys = [], []
    for i in range(N):
        M = 2 * math.pi * i / N
        xi = M
        for _ in range(80):
            xi -= (xi - e * math.sin(xi) - M) / (1 - e * math.cos(xi))
        xs.append(math.cos(xi) - e)
        ys.append(math.sqrt(1 - e * e) * math.sin(xi))
    out = []
    for k in range(1, K):
        cx = sum(xs[i] * cmath.exp(-2j * math.pi * k * i / N)
                 for i in range(N)) / N
        cy = sum(ys[i] * cmath.exp(-2j * math.pi * k * i / N)
                 for i in range(N)) / N
        out.append(k * k * math.sqrt(abs(cx) ** 2 + abs(cy) ** 2))
    return out


measured = {0.10: [1.0, 1.9957e-01, 3.3571e-02, 5.3029e-03,
                   7.9780e-04, 1.2593e-04],
            0.30: [1.0, 5.86583e-01, 2.90254e-01, 1.34527e-01,
                   6.0181e-02, 2.6359e-02],
            0.50: [1.0, 9.3686e-01, 7.4075e-01, 5.4860e-01,
                   3.9230e-01, 2.7452e-01],
            0.70: [1.0, 1.22215e+00, 1.26290e+00, 1.22326e+00,
                   1.1445e+00, 1.0481e+00]}

print('closed form against a pure Kepler ellipse:')
for e in (0.3, 0.7):
    d = kepler_dft(e)
    a1 = amp(1, e)
    row = ['%.6f' % ((d[k] / d[0]) / (amp(k + 1, e) / a1)) for k in range(6)]
    print('  e = %.2f  ' % e + '  '.join(row))

print('\nthe model against the closed form:')
for e in sorted(measured):
    a1 = amp(1, e)
    row = ['%.6f' % (measured[e][k] / (amp(k + 1, e) / a1))
           for k in range(6)]
    print('  e = %.2f  ' % e + '  '.join(row))
