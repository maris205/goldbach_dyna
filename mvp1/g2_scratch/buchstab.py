# Heuristic only: F_N(y) = e^{2g} * (1/N) int_0^N w(ln t/ln y) w(ln(N-t)/ln y) dt,  w = Buchstab omega
# compared with D2 count/(N*Haar) = 1 + E/(N*Haar).  Not a proof of anything.
import json, math, sys
import numpy as np
g = 0.5772156649015329
h = 1e-4
U = np.arange(1.0, 12.0 + h, h)
w = np.empty_like(U)
m2 = U <= 2.0
w[m2] = 1.0 / U[m2]
i0 = int(round(1.0 / h))
for i in range(len(U)):
    if U[i] <= 2.0: continue
    # (u w)' = w(u-1); trapezoid
    uw_prev = U[i-1] * w[i-1]
    uw = uw_prev + 0.5 * h * (w[i-1-i0] + w[i-i0])
    w[i] = uw / U[i]
def omega(u):
    u = np.asarray(u, float)
    out = np.interp(u, U, w)
    out[u < 1.0] = 0.0
    out[u > U[-1]] = math.exp(-g)
    return out
print('omega(2)=%.6f omega(2.5)=%.6f omega(3)=%.6f (1+ln2)/3=%.6f e^-g=%.6f' %
      (omega([2])[0], omega([2.5])[0], omega([3])[0], (1 + math.log(2)) / 3, math.exp(-g)))
for s in (2.0, 2.05, 2.1, 2.15, 2.165, 2.2, 2.3, 2.5, 3.0):
    print('  s=%.3f (e^g w(s))^2 = %.4f' % (s, (math.exp(g) * omega([s])[0]) ** 2))

def F(N, y):
    L = math.log(N); ly = math.log(y)
    # u in (0,1): t = N u ; use variable v = ln u on halves with symmetry
    v = np.linspace(-L, math.log(0.5), 200001)
    u = np.exp(v)
    f = omega((L + v) / ly) * omega((L + np.log1p(-u)) / ly) * u
    return 2 * math.exp(2 * g) * np.trapezoid(f, v)

d = json.load(open(sys.argv[1]))
for c in d['cases']:
    if c['N'] not in (10000, 1000000, 100000000): continue
    N = c['N']; print('N =', N)
    L = c['layers']; n = len(L)
    for r in L:
        if r['k'] in (n // 4, n // 2, int(0.7 * n), int(0.85 * n), n) or (c['critical'] and r['k'] == c['critical']['k']):
            y = r['p']
            print('  k=%d p=%d s=%.3f  D2 count/(N mu)=%.4f   F_N=%.4f   limit (e^g w(s))^2=%.4f' %
                  (r['k'], y, math.log(N) / math.log(y), 1 + r['E'] / (N * r['Haar']), F(N, y),
                   (math.exp(g) * omega([math.log(N) / math.log(y)])[0]) ** 2))
