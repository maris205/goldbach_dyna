# G3 说明性数值：u_*=1.4855245 处 f(x)=1-u x^2 的吸引周期轨道与反射 rho_N 的单调性
from decimal import Decimal as mpf, getcontext
getcontext().prec=60
class _M:
    @staticmethod
    def nstr(x,k): return "%.*e"%(k,float(x))
mp=_M()

u = mpf('1.4855245')
f = lambda x: 1 - u*x*x
x = mpf(0)
for _ in range(200000): x = f(x)
# 找周期
orb=[x]; y=x
for n in range(1,2000):
    y=f(y); orb.append(y)
per=None
for p in range(1,2000):
    if abs(orb[p]-orb[0])<mpf(10)**-25: per=p;break
print("period of attracting cycle:",per)
# 超稳定近似: 从 0 出发 f^210(0)
z=mpf(0)
for n in range(210): z=f(z)
print("f^210(0) =",mp.nstr(z,8))
# 用 0 的轨道 x_m=f^m(0)（m<210）作为编号 m mod 210 的点
xs=[mpf(0)]
for m in range(1,210): xs.append(f(xs[-1]))
# 检查 mod 2,6,30 的层: 不同剩余类的点落在不相交区间（凸包不交）
def disjoint(P):
    hulls=[]
    for a in range(P):
        pts=[xs[m] for m in range(210) if m%P==a]
        hulls.append((min(pts),max(pts)))
    hulls.sort()
    return all(hulls[i][1]<hulls[i+1][0] for i in range(len(hulls)-1))
for P in (2,6,30,210): print("level P=%d cylinders have disjoint hulls:"%P, disjoint(P))
# rho_N: a -> N-a mod 210; 是否在实数序下单调
order=sorted(range(210),key=lambda m: xs[m])
rank={m:i for i,m in enumerate(order)}
for N in list(range(0,12))+[100,101,209]:
    img=[rank[(N-m)%210] for m in order]
    inc=all(img[i]<img[i+1] for i in range(209)); dec=all(img[i]>img[i+1] for i in range(209))
    fix=[m for m in range(210) if (N-m)%210==m]
    print("N=%3d monotone_inc=%s monotone_dec=%s fixed=%s"%(N,inc,dec,fix))
# sigma(x)=-x: 是否有 x_m = -x_n
print("min |x_m + x_n| over m!=n, (m,n) not both 0:",
      mp.nstr(min(abs(xs[m]+xs[n]) for m in range(210) for n in range(210) if not(m==0 and n==0)),5))
best=min(((abs(xs[m]+xs[n]),m,n) for m in range(210) for n in range(210) if not(m==0 and n==0)))
print("argmin pair:",best[1],best[2],float(xs[best[1]]),float(xs[best[2]]))
for tol in (1e-6,1e-9,1e-12):
    p=next((p for p in range(1,2000) if abs(orb[p]-orb[0])<mpf(tol)),None); print("tol",tol,"period",p)
mult=mpf(1); y=orb[0]
for _ in range(210): mult*= -2*u*y; y=f(y)
print("multiplier of 210-cycle (approx):",float(mult))
