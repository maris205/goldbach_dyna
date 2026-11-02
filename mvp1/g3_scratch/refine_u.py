# 在 1.4855245 附近解 f_u^{210}(0)=0（超稳定），并核对最小周期、柱集结构、反射单调性
from decimal import Decimal as D, getcontext
getcontext().prec=80
def g(u,n=210):
    x=D(0)
    for _ in range(n): x=1-u*x*x
    return x
u0=D('1.4855245'); h=D('1e-9')
# 先扫描找变号
import itertools
grid=[u0+k*h for k in range(-2000,2001)]
vals=[g(u) for u in grid]
roots=[]
for i in range(len(grid)-1):
    if vals[i]==0 or vals[i]*vals[i+1]<0:
        a,b=grid[i],grid[i+1]; fa=vals[i]
        for _ in range(200):
            m=(a+b)/2; fm=g(m)
            if fa*fm<=0: b=m
            else: a,fa=m,fm
        roots.append((a+b)/2)
print("sign changes of f^210(0) in [u0-2e-6,u0+2e-6]:",len(roots))
for r in roots:
    xs=[D(0)]
    for m in range(1,211): xs.append(1-r*xs[-1]**2)
    minper=next(p for p in range(1,211) if abs(xs[p])<D('1e-40'))
    def disjoint(P):
        H=sorted((min(xs[m] for m in range(210) if m%P==a),max(xs[m] for m in range(210) if m%P==a)) for a in range(P))
        return all(H[i][1]<H[i+1][0] for i in range(P-1))
    order=sorted(range(210),key=lambda m: xs[m]); rank={m:i for i,m in enumerate(order)}
    nonmono=all(not(all(x<y for x,y in zip(I,I[1:])) or all(x>y for x,y in zip(I,I[1:])))
                for N in range(210) for I in [[rank[(N-m)%210] for m in order]])
    pm=min(abs(xs[m]+xs[n]) for m in range(210) for n in range(210) if (m,n)!=(0,0))
    if all(disjoint(P) for P in (2,6,30)): print("u=%s minper=%d disjoint(2,6,30)=%s rho_N nonmonotone for all N mod 210=%s min|x_m+x_n|=%.3e"
          %(str(r)[:14],minper,[disjoint(P) for P in (2,6,30)],nonmono,float(pm)))
