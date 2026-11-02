import json, math
from math import log, pi
d=json.load(open('/Users/maris/goldbach_sd/results/V9.json'))
c2=6-pi**2/6
# third coefficient: sum_{m+n=3} int (-log u)^m (-log(1-u))^n du = 2*6 + 2*int log^2u (-log(1-u))... compute numerically
from math import fsum
def integ(m,n,M=200000):
    # midpoint with substitution u = s^4 near 0? just use symmetric tanh-sinh-lite: u=(1-cos(pi x))/2
    s=0.0
    for i in range(M):
        x=(i+0.5)/M; u=(1-math.cos(pi*x))/2; du=pi*math.sin(pi*x)/2
        if u<=0 or u>=1: continue
        s+=((-log(u))**m)*((-log(1-u))**n)*du
    return s/M
c3=fsum(integ(m,3-m) for m in range(4))
print('c2',c2,'c3',c3)
for r in d['rows']:
    N=r['N'];L=log(N);I=r['I'];S=r['HL_integral']/I
    A=r['Haar_U_squared_logN_squared']
    B=N/(L*L*I)          # = 1/(I L^2/N)
    C=r['HL_integral']/r['r_ordered']   # S I / r
    RS=r['N_Haar_G']/(N*r['Haar_U']**2)/S   # R_k/S
    tot=A*B*C*RS
    print(N, 'IL2/N=%.5f'%(I*L*L/N), '1+2/L=%.5f'%(1+2/L), '3term=%.5f'%(1+2/L+c2/L**2), '4term=%.5f'%(1+2/L+c2/L**2+c3/L**3),
      'A=%.5f B=%.5f RS=%.8f C=%.5f tot=%.5f obs=%.5f'%(A,B,RS,C,tot,r['N_Haar_G_over_r']),
      'naive 4e^-2g/(1+2/L)=%.5f'%(1.2609470067/(1+2/L)), 'A*B=%.5f'%(A*B))
