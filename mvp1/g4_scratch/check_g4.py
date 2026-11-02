# g4_scratch/check_g4.py -- G4 自检（非 D5 官方结果）
# 1) N_k^* = 2*3^(k-2)/m_k, m_k = 1/2 * prod_{3<=p<=p_k}(1-2/p)   (k>=2)
# 2) 对 k<=12、偶数 4<=N<=ceil(N_k^*)+2：暴力找 Cnt_k(N)=0 的 N（逐个 n 检查 gcd，早停）
# 3) 检查：最大例外 <= N_k^*；并记录最大例外与 p_k^2 的比较
# 4) 顶层 Möbius 恒等式 sum_{0<n<N} mu(n) 1_G(n) = 2*1[N-1 prime] - Cnt   (N<=3000 全量)
from fractions import Fraction
from math import gcd, isqrt
import json, sys

PR = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]

def Nstar(k):
    m = Fraction(1, 2)
    for p in PR[1:k]:
        m *= Fraction(p - 2, p)
    return Fraction(2 * 3 ** (k - 2)) / m, m

def exceptions(k, Nmax):
    P = 1
    for p in PR[:k]:
        P *= p
    cop = [n for n in range(1, Nmax) if gcd(n, P) == 1]
    out = []
    for N in range(4, Nmax + 1, 2):
        hit = False
        for n in cop:
            if n >= N:
                break
            if gcd(N - n, P) == 1:
                hit = True
                break
        if not hit:
            out.append(N)
    return out

def sieve(M):
    s = bytearray([1]) * (M + 1); s[0] = s[1] = 0
    for i in range(2, isqrt(M) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return s

def mobius(M):
    mu = [1] * (M + 1); isp = sieve(M)
    for p in range(2, M + 1):
        if isp[p]:
            for j in range(p, M + 1, p): mu[j] = -mu[j]
            for j in range(p * p, M + 1, p * p): mu[j] = 0
    return mu, isp

res = {}
for k in range(2, 13):
    ns, m = Nstar(k)
    Nmax = int(ns) + 2
    if Nmax % 2: Nmax += 1
    ex = exceptions(k, Nmax)
    res[k] = dict(p_k=PR[k - 1], m_k=float(m), Nstar=float(ns), n_exc=len(ex),
                  max_exc=(max(ex) if ex else None), pk2=PR[k - 1] ** 2,
                  bound_ok=(not ex or max(ex) <= ns), exc_head=ex[:40])
    print(k, res[k], file=sys.stderr, flush=True)

M = 3000
mu, isp = mobius(M)
bad = 0
for N in range(4, M + 1):
    k = sum(1 for p in range(2, isqrt(N) + 1) if isp[p])
    ps = [p for p in range(2, isqrt(N) + 1) if isp[p]]
    G = [n for n in range(1, N) if all(n % p and (N - n) % p for p in ps)]
    lhs = sum(mu[n] for n in G)
    rhs = 2 * isp[N - 1] - len(G)
    bad += (lhs != rhs)
res['mobius_top_identity_N<=3000_failures'] = bad
json.dump(res, open('g4_scratch/check_g4.json', 'w'), indent=1)
print('mobius failures', bad)
