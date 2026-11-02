# Check identity Cnt_top - r = 2*[N-1 prime] - 2*#{p<=sqrt N: N-p prime} against D2/D3 fields.
import json, math
def ip(n):
    if n < 2: return False
    for p in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n % p == 0: return n == p
    d, s = n-1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2,3,5,7,11,13,17):
        x = pow(a, d, n)
        if x in (1, n-1): continue
        for _ in range(s-1):
            x = x*x % n
            if x == n-1: break
        else: return False
    return True
D2 = json.load(open('results/D2.json')); D3 = json.load(open('results/D3.json'))
hs = {x['N']: x['r'] for x in D3['high_samples']}
for c in D2['cases']:
    N = c['N']; top = c['layers'][-1]; rt = math.isqrt(N)
    pred = 2*ip(N-1) - 2*sum(1 for p in range(2, rt+1) if ip(p) and ip(N-p))
    print(N, 'top_k', top['k'], 'pi(sqrtN) bound', -2*top['k'], 'identity', pred, 'D2-D3', (top['count']-hs[N]) if N in hs else '-')
