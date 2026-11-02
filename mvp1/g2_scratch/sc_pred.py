# Heuristic only: first crossing of 1-F_N(y) > 0.1 when y increases (s = ln N/ln y decreases).
import math, numpy as np, importlib.util, sys
spec = importlib.util.spec_from_file_location('b', 'g2_scratch/buchstab.py')
src = open('g2_scratch/buchstab.py').read().split('d = json.load')[0]
ns = {}; exec(src, ns)
F = ns['F']; omega = ns['omega']; g = ns['g']
# limit crossing
ss = np.linspace(2.0, 3.0, 100001)
lim = (math.exp(g) * omega(ss)) ** 2
print('limit s* where (e^g w)^2=0.9 :', ss[np.argmin(abs(lim - 0.9))])
for m in (4, 5, 6, 7, 8, 10, 12, 16, 24):
    N = 10.0 ** m; L = math.log(N)
    sc = None
    for s in np.arange(4.0, 1.999, -0.0025):
        y = math.exp(L / s)
        if 1 - F(N, y) > 0.1:
            sc = s; break
    top = 1 - F(N, math.sqrt(N))
    peak = max((1 - F(N, math.exp(L / s)), s) for s in np.arange(2.0, 2.6, 0.005))
    print('N=1e%d  s_c^F=%s  c=N^(1/s-1/2)=%s  1-F at s=2: %.4f  peak %.4f at s=%.3f' %
          (m, None if sc is None else '%.4f' % sc, None if sc is None else '%.3f' % N ** (1 / sc - .5), top, peak[0], peak[1]))
