from math import gcd, log, isqrt
from fractions import Fraction as F

def primes(n):
    s = bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2, isqrt(n)+1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]

PR = primes(200000)
# G1.1 set identities mod P_k
for k in range(1, 6):
    ps = PR[:k]; P = 1
    for p in ps: P *= p
    U = {x for x in range(P) if all(x % p for p in ps)}
    for N in range(0, 2*P):
        A = {x for x in U if (N-x) % P in U}
        B = {x for x in U if (x-N) % P in U}
        C = {x for x in U if (x+N) % P in U}
        assert A == B
        assert C == {(x-N) % P for x in B}
        assert C == {(-x) % P for x in B}
        # G1.2 ratio
        R = F(len(A)*P, len(U)**2)
        f = F(1)
        for p in ps:
            f *= F(p-1, p-2) * F(p*(p-2), (p-1)**2) if N % p == 0 and p > 2 else (F(p*(p-2), (p-1)**2) if p > 2 else (2 if N % 2 == 0 else 0))
        assert R == f, (k, N, R, f)
    # G1.3 Haar average
    avg = sum(F(len({x for x in U if (N-x) % P in U})*P, len(U)**2) for N in range(P)) / P
    assert avg == 1
print("G1.1/G1.2/G1.3(Haar) exact checks OK for k<=5")

# finite window bijection
for N in range(4, 3000, 2):
    r = isqrt(N); ps = [p for p in PR if p <= r]; Pset = set(PR)
    mid = sum(1 for n in range(r+1, N-r) if all(n % p and (N-n) % p for p in ps))
    rN = sum(1 for n in range(2, N-1) if n in Pset and N-n in Pset)
    assert mid == sum(1 for n in range(r+1, N-r) if n in Pset and N-n in Pset)
    assert 0 <= rN - mid <= 2*len(ps)
print("finite window OK")

# G1.3 natural density, fixed k: sum_{even N<=X} R_k(N) = X - theta, 0<=theta<=2 prod p/(p-1)
def Rk(N, ps):
    if N % 2: return 0.0
    v = 2.0
    for p in ps[1:]:
        v *= p*(p-2)/(p-1)**2
        if N % p == 0: v *= (p-1)/(p-2)
    return v
for k in (2, 4, 6, 10):
    ps = PR[:k]; bound = 2
    for p in ps[1:]: bound *= p/(p-1)
    for X in (10, 101, 1000, 12345):
        th = X - sum(Rk(N, ps) for N in range(2, X+1, 2))
        assert -1e-9 <= th <= bound + 1e-9, (k, X, th, bound)
print("fixed-k density bound OK")

# singular series sums
C2 = 1.0
for p in PR[1:]: C2 *= 1 - 1/(p-1)**2
print("C2 approx", C2)
def S(N):
    if N % 2: return 0.0
    v = 2*C2; m = N
    while m % 2 == 0: m //= 2
    for p in PR[1:]:
        if p*p > m: break
        if m % p == 0:
            v *= (p-1)/(p-2)
            while m % p == 0: m //= p
    if m > 1: v *= (m-1)/(m-2)
    return v
for X in (10**3, 10**4, 10**5, 10**6):
    s = sum(S(N) for N in range(2, X+1, 2))
    print(X, s - X, (s - X)/log(X))
