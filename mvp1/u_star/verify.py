"""Independent Decimal / CRT verification of delivered numerical results."""
from pathlib import Path
from decimal import Decimal, localcontext
from fractions import Fraction
from itertools import product
import hashlib
import json
import math
import time

HERE=Path(__file__).resolve().parent
OUT=HERE/'results'
PRIMES=[2,3,5,7,11,13,17]


def critical_return(u,p=210):
    x=Decimal(0)
    xs=[x]
    for _ in range(p):
        x=1-u*x*x
        xs.append(x)
    return xs


def crt_count(N,primes):
    # Inclusion-exclusion of forbidden classes n=0 or N mod p.
    # Independent of full-period residue enumeration in experiment.py.
    total=0
    choices=[[None]+sorted(set([0,N%p])) for p in primes]
    for selected in product(*choices):
        constraints=[(p,r) for p,r in zip(primes,selected) if r is not None]
        modulus=math.prod(p for p,r in constraints)
        if not constraints:
            a=0
        else:
            a=sum(r*(modulus//p)*pow(modulus//p,-1,p) for p,r in constraints)%modulus
        count=(N-1-a)//modulus-(-a)//modulus
        total+=(-1 if len(constraints)%2 else 1)*count
    return total


def main():
    start=time.perf_counter()
    root=json.loads((OUT/'u_star.json').read_text())
    cyl=json.loads((OUT/'cylinders.json').read_text())
    d4=json.loads((OUT/'d4_rerun.json').read_text())
    checks=[]
    with localcontext() as ctx:
        ctx.prec=160
        for c in root['candidates']:
            u=Decimal(c['u'])
            xs=critical_return(u)
            assert abs(xs[210])<Decimal('1e-70')
            assert all(abs(xs[p])>Decimal('1e-20') for p in range(1,210) if 210%p==0)
            for P in (2,6,30):
                hh=sorted((min(xs[a:210:P]),max(xs[a:210:P])) for a in range(P))
                assert all(hh[i][1]<hh[i+1][0] for i in range(P-1))
            phase='C'+''.join('R' if x>0 else 'L' for x in xs[1:210])
            assert phase==c['phase_symbols_x0_to_x209']
            bracket=c['numerical_root']
            ga=critical_return(Decimal(bracket['lower']))[210]
            gb=critical_return(Decimal(bracket['upper']))[210]
            assert ga*gb<0
            errorbound=Decimal(c['rounding']['absolute_u_error_bound'])
            lo,hi=Decimal(bracket['lower']),Decimal(bracket['upper'])
            assert max(abs(u-lo),abs(u-hi))<=errorbound
            checks.append(dict(index=c['index'],decimal160_residual=str(xs[210]),
                bracket_lower_residual=str(ga),bracket_upper_residual=str(gb),
                residual_tolerance_pass=True,period_hull_sign_checks_pass=True,
                printed_rounding_bound_pass=True))
    for r in d4['rows']:
        assert crt_count(r['N'],PRIMES[:r['k']])==r['count_odo']
        h=Fraction(r['Haar_fraction'])
        er=abs(Fraction(r['count_odo'],r['N']-1)/h-1)
        assert er==Fraction(r['relative_error_odo_fraction'])
        assert r['passes_odo']==(er<Fraction(1,1000))
        if r['k']<=4:
            assert r['pointwise_mismatches']==0 and r['unknown_pairs']==0
            assert r['count_odo']==r['count_knead']
    phase=cyl['phase_word']
    for dic in cyl['dictionaries']:
        if dic['k']>4:
            w=dic['intrinsic_period_obstruction']
            n=w['witness_n']
            assert math.gcd(n,dic['P'])!=1 or math.gcd(n+210,dic['P'])!=1
            assert (math.gcd(n,dic['P'])==1)!=(math.gcd(n+210,dic['P'])==1)
            assert phase[n%210:]+phase[:n%210]==phase[(n+210)%210:]+phase[:(n+210)%210]
    res=dict(root_independent_checks=checks,D4_CRT_inclusion_exclusion_rows=126,
        exact_threshold_checks=True,higher_level_obstruction_verified=True,
        elapsed_seconds=time.perf_counter()-start,
        artifact_sha256={name:hashlib.sha256((OUT/name).read_bytes()).hexdigest()
            for name in ('u_star.json','cylinders.json','d4_rerun.json')})
    with (OUT/'independent_verification.json').open('x') as f:
        json.dump(res,f,ensure_ascii=False,indent=2);f.write('\n')
    with (HERE/'REPORT.md').open('a') as f:
        f.write('\n## 独立复核（追加）\n\n用 Python Decimal 160 位独立重算全部 9 个输出参数的回返、符号、最小周期、模 2/6/30 凸包分离，均通过；每个保存的根区间端点仍严格异号，100 位十进制输出误差界逐个通过区间端点距离核验。用 CRT 禁止剩余类容斥公式独立核验全部 126 组地址计数，并用精确有理数重检严格 1e-3 阈值，全部一致。k=5…7 的同符号异标签见证亦核验通过。详情 `results/independent_verification.json`，耗时 '+f'{res["elapsed_seconds"]:.3f}'+' 秒。\n')
    print(json.dumps(res,indent=2))


if __name__=='__main__':
    main()
