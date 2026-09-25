# 202c/202d: the inspiral time integral in coordinate time and in the
# particles' proper time, and where along the radius it accrues.
# The E1 powers are audit 201b's measurements.
#
# Run: python3 frames.py
import math

c = 2.99792458e8
k = 8.9875517873681764e9 * (1.602176634e-19) ** 2
mu = 9.1093837139e-31 / 2.0
a = 1.0583544190745552e-10

# (r/a_pair, E1 para, E1 ortho) from audit 201b
rows = sorted([(1.0, 1.1673121110e-08, 1.1673009177e-08),
               (0.3, 1.4411508358e-06, 1.4409972454e-06),
               (0.1, 1.1677263065e-04, 1.1666036128e-04),
               (0.03, 1.4485722905e-02, 1.4327070408e-02),
               (0.01, 1.2473157007e+00, 1.1052770488e+00),
               (0.003, 3.7360634350e+02, 7.0996330697e+01),
               (0.001, 3.4764632528e+04, 1.7626547589e+03)])


def gamma(rf):
    beta = math.sqrt(k / (mu * rf * a)) / 2.0 / c
    return 1.0 / math.sqrt(1.0 - beta * beta), beta


def integrate(idx, proper):
    total = 0.0
    for i in range(len(rows) - 1):
        seg = []
        for r, power in ((rows[i][0], rows[i][idx]),
                         (rows[i + 1][0], rows[i + 1][idx])):
            w = 1.0 / (2 * r * r * power)
            if proper:
                w /= gamma(r)[0]
            seg.append(w)
        total += 0.5 * (seg[0] + seg[1]) * (rows[i + 1][0] - rows[i][0])
    return total


print('  %10s %12s %14s' % ('r/a_pair', 'beta', 'gamma-1'))
for rf, _, _ in sorted(rows, reverse=True):
    g, b = gamma(rf)
    print('  %10.3g %12.5e %14.5e' % (rf, b, g - 1.0))

print()
for proper in (False, True):
    tp, to = integrate(1, proper), integrate(2, proper)
    label = 'proper time' if proper else 'coordinate time'
    print('%-18s t_ortho/t_para - 1 = %.6e' % (label, to / tp - 1.0))

print('\nwhere the coordinate time accrues (para):')
total = integrate(1, False)
for i in range(len(rows) - 1):
    r1, r2 = rows[i][0], rows[i + 1][0]
    seg = 0.5 * (1.0 / (2 * r1 * r1 * rows[i][1])
                 + 1.0 / (2 * r2 * r2 * rows[i + 1][1])) * (r2 - r1)
    print('  %7.3g -> %-7.3g %9.4f%%' % (r1, r2, 100 * seg / total))
