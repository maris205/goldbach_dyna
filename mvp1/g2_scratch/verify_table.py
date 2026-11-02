# Recompute G2 table columns from D2.json + definitions (static check).
import json, math
D2 = json.load(open('results/D2.json'))
def primes(n):
    s = bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2, int(n**.5)+1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]
P = primes(20000)
def beta_sharp(N, k):
    v = 6/N*3**(k-2)
    for p in P[2:k]: v /= (1-2/p)
    return v
def beta_N(N, k):
    v = 6/N
    for p in P[2:k]: v *= (2*p/(p-1) if N % p == 0 else 3*p/(p-2))
    return v
def kmax(N, f, eps=.1):
    k = 1
    while f(N, k+1) <= eps: k += 1
    return k
for m in range(4, 13):
    N = 10**m; k = kmax(N, beta_sharp)
    print('table', m, k, P[k-1], round(sum(math.log(p) for p in P[:k])/math.log(N), 2), round(math.log(math.log(N))/math.log(3), 2))
for c in D2['cases']:
    N = c['N']; L = c['layers']; top = L[-1]
    k1 = kmax(N, beta_sharp); k2 = kmax(N, beta_N)
    mx = max(L, key=lambda r: r['relative_error'])
    print(N, k1, round(sum(math.log(p) for p in P[:k1])/math.log(N), 2), k2, '| top', top['k'], round(top['relative_error'], 4),
          '| max', round(mx['relative_error'], 4), round(math.log(N)/math.log(mx['p']), 2), '| rho at k2', '%.2e' % L[k2-1]['relative_error'])
