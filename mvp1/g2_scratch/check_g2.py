# G2 scratch: verify Prop 2.1 bounds by brute force, compute k_max(N,eps),
# compare with results/D2.json (read-only), and evaluate the Buchstab heuristic.
import json, math, sys
from fractions import Fraction

def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]

PR = primes_upto(200000)

def nu(p, N): return 1 if N % p == 0 else 2

def Tstar(k, N):
    t = 1
    for p in PR[:k]: t *= min(p - nu(p, N), 1 + nu(p, N))
    return t

def Tleg(k, N):
    t = 1
    for p in PR[:k]: t *= 1 + nu(p, N)
    return t

def haar(k, N):
    h = Fraction(1)
    for p in PR[:k]: h *= Fraction(p - nu(p, N), p)
    return h

# 1. brute force: arbitrary windows [x0, x0+L), N even/odd small
def brute():
    worst = 0.0
    for N in list(range(2, 200)) + [210, 2310, 30030, 1024, 999, 4096]:
        for k in range(1, 6):
            Pk = math.prod(PR[:k]); h = haar(k, N)
            G = [1 if all(n % p and (N - n) % p for p in PR[:k]) else 0 for n in range(Pk)]
            ts, tl = Tstar(k, N), Tleg(k, N)
            for x0 in (0, 1, 7, Pk // 3):
                for L in list(range(0, min(Pk, 400))) + [Pk + 5, 2 * Pk + 77]:
                    c = sum(G[(x0 + i) % Pk] for i in range(L))
                    E = abs(c - L * h)
                    assert E <= ts, (N, k, x0, L, E, ts)
                    assert E <= tl - 1 or tl == 1, (N, k, E, tl)
                    if ts: worst = max(worst, float(E) / ts)
    print('brute force OK; max |E|/T* =', worst)

def kmax(N, eps, perN=False):
    # largest k>=2 with beta_k <= eps, beta = T*/(N haar) (uniform worst case if not perN)
    k = 1; best = None
    while True:
        k += 1
        if perN: b = Tstar(k, N) / (N * float(haar(k, N)))
        else:
            b = 6.0 / N
            for p in PR[2:k]: b *= 3 * p / (p - 2)
        if b <= eps: best = (k, PR[k - 1], b)
        else: return best

if __name__ == '__main__':
    brute()
    print('\nk_max table (eps=0.1, uniform worst case):')
    for m in range(4, 13):
        N = 10 ** m
        k, p, b = kmax(N, 0.1)
        lP = sum(math.log(q) for q in PR[:k])
        print(f'N=1e{m}: k_max={k} p={p} beta={b:.4f} x=lnP/lnN={lP/math.log(N):.3f}  '
              f'asym (ln(epsN)-2lnlnlnN)/ln3={(math.log(0.1*N)-2*math.log(math.log(math.log(N))))/math.log(3):.2f}'
              f'  lnlnN/ln3={math.log(math.log(N))/math.log(3):.2f}')
    d = json.load(open(sys.argv[1]))
    print('\nD2 comparison:')
    viol = 0; tot = 0
    for c in d['cases']:
        N = c['N']
        for r in c['layers']:
            tot += 1
            T = Tstar(r['k'], N)
            if T < 10**15 and abs(r['E']) > T + 1e-6: viol += 1
        ku = kmax(N, 0.1); kn = kmax(N, 0.1, perN=True)
        # last k where |E| <= T* is informative (T*/main <=1)
        first_ratio_gt = next((r['k'] for r in c['layers'] if Tstar(r['k'], N) > N * r['Haar']), None)  # int vs float ok for small
        # actual |E|/T* at k = k_max
        lay = {r['k']: r for r in c['layers']}
        rk = lay[kn[0]]
        print(N, 'k_max uni', ku[:2], 'perN', kn[:2], 'k_c', c['critical'] and c['critical']['k'],
              'T*>main from k=', first_ratio_gt,
              'rel@kmax=%.2e |E|/T*@kmax=%.3f' % (rk['relative_error'], abs(rk['E']) / Tstar(rk['k'], N)))
    print('violations of |E|<=T*:', viol, 'of', tot)
