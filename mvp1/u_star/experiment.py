"""Recover finite symbolic dictionaries and rerun D4, without network or parent writes."""
from pathlib import Path
from fractions import Fraction
from collections import defaultdict
import hashlib
import json
import math
import time
import mpmath as mp
import numpy as np

HERE=Path(__file__).resolve().parent
OUT=HERE/'results'
PRIMES=[2,3,5,7,11,13,17]
PERIODS=[math.prod(PRIMES[:k]) for k in range(1,8)]
NS=[100,128,256,512,1000,1024,2048,4096,8192,10000,10002,12000,
    16384,20000,30030,60060,90090,100000]
HORIZON=2*PERIODS[-1]
GUARD=mp.mpf('1e-70')
mp.mp.dps=120


def save(name,data):
    with (OUT/name).open('x',encoding='utf8') as f:
        json.dump(data,f,ensure_ascii=False,indent=2,allow_nan=False)
        f.write('\n')


def stat(decoded,truth,start=210):
    known=decoded>=0
    bad=known & (decoded!=truth)
    unknown=~known
    failed=bad|unknown
    hold=np.arange(len(truth))>=start
    def first(mask):
        ids=np.flatnonzero(mask)
        return int(ids[0]) if len(ids) else None
    return dict(size=len(truth),known=int(known.sum()),unknown=int(unknown.sum()),
        mismatches=int(bad.sum()),matching_rate=float((known & ~bad).sum()/len(truth)),
        known_matching_rate=float((known & ~bad).sum()/max(1,known.sum())),
        first_mismatch=first(bad),first_unknown=first(unknown),first_failure=first(failed),
        heldout_start=start,heldout_size=int(hold.sum()),
        heldout_unknown=int((hold & unknown).sum()),heldout_mismatches=int((hold & bad).sum()),
        heldout_matching_rate=float((hold & known & ~bad).sum()/hold.sum()))


def long_orbit(u,length,dps=120,guard='1e-70',reference=None):
    with mp.workdps(dps):
        u=mp.mpf(u)
        e=mp.mpf(guard)
        x=mp.mpf(0)
        symbols=bytearray()
        maxcritical=mp.mpf(0)
        min_noncritical=mp.inf
        maxdiff=mp.mpf(0)
        maxabs=mp.mpf(0)
        base=[]
        for n in range(length):
            symbols.append(ord('C') if abs(x)<e else ord('R') if x>0 else ord('L'))
            maxabs=max(maxabs,abs(x))
            if n<210:
                base.append(x)
            if n%210==0:
                maxcritical=max(maxcritical,abs(x))
            else:
                min_noncritical=min(min_noncritical,abs(x))
            if n>=210:
                maxdiff=max(maxdiff,abs(x-base[n%210]))
            x=1-u*x*x
        info=dict(length=length,dps=dps,guard=guard,max_abs=str(maxabs),
            max_critical_phase_abs=str(maxcritical),min_noncritical_phase_abs=str(min_noncritical),
            max_difference_from_first_period=str(maxdiff),x_after_length=str(x),
            finite=bool(mp.isfinite(x)))
        if reference is not None:
            compare=np.frombuffer(bytes(symbols),dtype=np.uint8)!=np.frombuffer(reference[:length],dtype=np.uint8)
            ids=np.flatnonzero(compare)
            info.update(symbol_mismatches_to_periodic=int(compare.sum()),
                        first_mismatch_to_periodic=int(ids[0]) if len(ids) else None)
        return bytes(symbols),info


def train(raw,labels,L,exclude_C=False):
    groups=defaultdict(set)
    phases=defaultdict(list)
    for n in range(210):
        key=raw[n:n+L]
        groups[key].add(int(labels[n]))
        phases[key].append(n)
    mapping={key:next(iter(v)) if len(v)==1 and not (exclude_C and b'C' in key) else -1
             for key,v in groups.items()}
    entries=[dict(word=key.decode(),label=mapping[key],training_labels=sorted(v),
                  training_phases=phases[key]) for key,v in groups.items()]
    return mapping,dict(word_length=L,training_positions=[0,209],training_size=210,
        distinct_words=len(groups),conflicting_words=sum(len(v)>1 for v in groups.values()),
        C_words=sum(b'C' in key for key in groups),exclude_C=exclude_C,entries=entries)


def decode(raw,mapping,L,length):
    # Runtime decoder gets only the observed word. No n, gcd, modulo or address argument.
    return np.fromiter((mapping.get(raw[i:i+L],-1) for i in range(length)),dtype=np.int8,count=length)


def direct_sign_test(raw,P):
    s=np.frombuffer(raw,dtype=np.uint8)
    truth=np.gcd(np.arange(len(s),dtype=np.int64),P)==1
    literal=np.where(s==ord('L'),1,0).astype(np.int8)  # sieve L=survive, R=delete; C=delete.
    reverse=np.where(s==ord('R'),1,0).astype(np.int8)
    return dict(convention='Sieve L=survive/R=delete; logistic C treated as delete for this binary diagnostic',
        L_as_survive=stat(literal,truth),R_as_survive=stat(reverse,truth))


def pair_summary(a,n,N,truth):
    left,right=a[n],a[N-n]
    good=(left==1)&(right==1)
    bad=(left==0)|(right==0)
    known=good|bad
    mismatch=known & (good!=truth)
    unknown=~known
    count=int(good.sum())
    return dict(count_knead=count if not unknown.any() else None,
        freq_knead=count/len(n) if not unknown.any() else None,
        decoder_frequency_lower=count/len(n),
        decoder_frequency_upper=(count+int(unknown.sum()))/len(n),
        unknown_pairs=int(unknown.sum()),pointwise_mismatches=int(mismatch.sum()),
        first_mismatch_n=int(n[np.flatnonzero(mismatch)[0]]) if mismatch.any() else None,
        first_unknown_n=int(n[np.flatnonzero(unknown)[0]]) if unknown.any() else None,
        all_match=not unknown.any() and not mismatch.any())


def sign_sensitivity(u,phase):
    length=20001
    ss=np.frombuffer((phase*((length+209)//210))[:length].encode(),dtype=np.uint8)
    noncrit=np.arange(length)%210!=0
    out=[]
    params=[('truncated_7dp','1.4855245',None,None)]
    for d in (8,9,10,11,12,14,20,40):
        for sg in (-1,1):
            params.append((f'u_star{sg:+d}e-{d}',str(mp.mpf(u)+sg*mp.power(10,-d)),d,sg))
    # Adjacent float64 parameter approximation is tested separately below.
    for label,param,d,sg in params:
        with mp.workdps(110):
            x=mp.mpf(0); uu=mp.mpf(param)
            first_raw=None;first_nc=None;first_guard=None
            ncdiff=0;guarddiff=0
            taildiff=0
            for n in range(length):
                raw=ord('C') if x==0 else ord('R') if x>0 else ord('L')
                sym=ord('C') if abs(x)<mp.mpf('1e-70') else raw
                if raw!=ss[n] and first_raw is None:
                    first_raw=n
                if sym!=ss[n]:
                    guarddiff+=1
                    if first_guard is None:first_guard=n
                if noncrit[n] and raw!=ss[n]:
                    ncdiff+=1
                    if n>=10000:taildiff+=1
                    if first_nc is None:first_nc=n
                x=1-uu*x*x
            out.append(dict(label=label,u=param,delta_exponent=d,delta_sign=sg,length=length,
                first_sign_difference_raw=first_raw,
                first_sign_difference_guarded=first_guard,
                guarded_differences=guarddiff,
                first_LR_difference_noncritical=first_nc,
                LR_differences_noncritical=ncdiff,noncritical_tail_differences=taildiff,
                guard='1e-70',dps=110))
    uf=np.float64(u); x=np.float64(0)
    actual=bytearray()
    for _ in range(max(NS)+211):
        actual.append(ord('C') if abs(x)<1e-12 else ord('R') if x>0 else ord('L'))
        x=np.float64(1)-uf*x*x
    observed=np.frombuffer(bytes(actual),dtype=np.uint8)
    expected=np.frombuffer((phase*((len(actual)+209)//210))[:len(actual)].encode(),dtype=np.uint8)
    mm=observed!=expected
    nc=(np.arange(len(actual))%210!=0)&mm
    f64=dict(u_float64=repr(float(uf)),guard='1e-12',length=len(actual),
        first_sign_difference=int(np.flatnonzero(mm)[0]) if mm.any() else None,
        symbol_mismatches=int(mm.sum()),
        first_LR_difference_noncritical=int(np.flatnonzero(nc)[0]) if nc.any() else None,
        LR_differences_noncritical=int(nc.sum()))
    return out,f64,bytes(actual)


def main():
    start=time.perf_counter()
    root=json.loads((OUT/'u_star.json').read_text())
    chosen=root['candidates'][root['selected_index']-1]
    phase=chosen['phase_symbols_x0_to_x209']
    u=root['u']
    length=max(NS)+211
    ref=(phase*((HORIZON+420)//210)).encode()
    actual,orbitinfo=long_orbit(u,length,reference=ref)
    assert orbitinfo['symbol_mismatches_to_periodic']==0
    actual2,precisioninfo=long_orbit(u,5001,dps=170,reference=ref)
    assert actual2==actual[:5001]
    print('High precision orbit:',json.dumps(orbitinfo),flush=True)
    # Beyond the directly iterated 100211 points, use the finite root's symbolic cycle.
    raw=ref[:HORIZON+210]
    indices=np.arange(HORIZON,dtype=np.int64)
    truth=[(np.gcd(indices,P)==1).astype(np.int8) for P in PERIODS]
    dictionaries=[];decoded=[];controls=[]
    for k,(P,t) in enumerate(zip(PERIODS,truth),1):
        best=None
        for L in range(1,211):
            mapping,info=train(raw,t,L)
            if info['conflicting_words']==0:
                best=(L,mapping,info)
                break
        assert best is not None
        L,mapping,info=best
        a=decode(raw,mapping,L,HORIZON)
        direct=decode(actual,mapping,L,max(NS)+1)
        assert np.array_equal(direct,a[:len(direct)])
        diag=stat(a,t)
        info.update(k=k,P=P,diagnostics=diag,
            evaluation_source='First 100001 decoded from directly iterated high-precision orbit and checked against periodic extension; rest from the verified 210-symbol cycle',
            direct_iterated_diagnostics=stat(direct,t[:len(direct)]),
            direct_sign_diagnostic=direct_sign_test(raw[:HORIZON],P),
            valid_as_periodic_sieve_dictionary=bool(k<=4 and diag['mismatches']==0))
        if k>4:
            phase_conflicts=[]
            for r in range(210):
                vals=t[r:P:210]
                if vals.min()!=vals.max():phase_conflicts.append(r)
            witness=next(n for n in range(P-210) if t[n]!=t[n+210])
            q=P//210
            phi_q=sum(math.gcd(i,q)==1 for i in range(q))
            lower=Fraction(48*min(phi_q,q-phi_q),P)
            info.update(intrinsic_period_obstruction=dict(
                same_symbols_at_n_and_n_plus_210=True,witness_n=witness,
                witness_n_plus_210=witness+210,label_n=int(t[witness]),label_next=int(t[witness+210]),
                mixed_phases=phase_conflicts,best_possible_error_rate_fraction=str(lower),
                best_possible_matching_rate=float(1-lower),
                reason='Identical entire symbolic futures at n and n+210 may carry different gcd labels; no longer word can remove this collision.'))
        dictionaries.append(info);decoded.append(a)
        cc=[]
        for LL,exclude in ((32,False),(32,True),(L,True)):
            mm,ii=train(raw,t,LL,exclude)
            aa=decode(actual,mm,LL,max(NS)+1)
            cc.append(dict(word_length=LL,exclude_C=exclude,
                distinct_words=ii['distinct_words'],conflicting_words=ii['conflicting_words'],
                diagnostic=stat(aa,t[:len(aa)])))
        controls.append(dict(k=k,P=P,variants=cc))
        print('dictionary',k,P,'L=',L,json.dumps(diag),flush=True)
    # A phase dictionary has 210 unique words at most. Explicitly verify phase recovery.
    for L in range(1,211):
        words={raw[n:n+L]:n for n in range(210)}
        if len(words)==210:break
    phase_mapping={raw[n:n+L]:n for n in range(210)}
    phase_recovery=np.fromiter((phase_mapping[actual[n:n+L]] for n in range(max(NS)+1)),
                              dtype=np.int16,count=max(NS)+1)
    assert np.array_equal(phase_recovery,np.arange(max(NS)+1)%210)
    cyl=dict(u=u,phase_word=phase,training='Only 210 phase labels; observed word at test time, no test address fallback',
        C_convention='Numerical |x|<1e-70 is C; min noncritical |x| and critical error are separately recorded',
        dictionaries=dictionaries,controls=controls,orbit_validation=orbitinfo,
        precision_validation=precisioninfo,phase_recovery=dict(word_length=L,
            distinct_words=210,checked=max(NS)+1,mismatches=0,
            interpretation='Finite 210-phase decoding with a labeled training cycle; does not establish an infinite odometer conjugacy'),
        horizon=HORIZON,code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    save('cylinders.json',cyl)
    rows=[]
    for N in NS:
        n=np.arange(1,N,dtype=np.int64)
        for k,(P,t,a) in enumerate(zip(PERIODS,truth,decoded),1):
            good=(t[n]==1)&(t[N-n]==1)
            count=int(good.sum())
            rs=np.arange(P,dtype=np.int64)
            residues=rs[(np.gcd(rs,P)==1)&(np.gcd(N-rs,P)==1)]
            exact=int(np.sum((N-1-residues)//P-(-residues)//P))
            h=Fraction(1)
            for p in PRIMES[:k]:h*=Fraction(p-(1 if N%p==0 else 2),p)
            assert Fraction(len(residues),P)==h and exact==count
            ps=pair_summary(a,n,N,good)
            fq=count/(N-1)
            err=abs(Fraction(count,N-1)/h-1)
            fk=ps['freq_knead']
            er_k=None if fk is None else abs(Fraction(ps['count_knead'],N-1)/h-1)
            rows.append(dict(N=N,k=k,P=P,sample_size=N-1,method='full enumeration',
                freq_odo=fq,count_odo=count,exact_count=exact,Haar=float(h),Haar_fraction=str(h),
                absolute_error_odo=abs(fq-float(h)),relative_error_odo=float(err),
                relative_error_odo_fraction=str(err),passes_odo=err<Fraction(1,1000),
                absolute_error_knead=None if fk is None else abs(fk-float(h)),
                relative_error_knead=None if er_k is None else float(er_k),
                relative_error_knead_fraction=None if er_k is None else str(er_k),
                passes_knead=er_k is not None and er_k<Fraction(1,1000),
                interpretation='finite k<=4 dictionary' if k<=4 else 'invalid higher-level extrapolation diagnostic',
                **ps))
    big=[r for r in rows if r['N']>=10000 and r['k']<=4]
    summary=dict(NS_count=len(NS),rows=len(rows),primary_rows=len(NS)*4,
        point_pairs_primary=sum(N-1 for N in NS)*4,
        all_primary_pointwise_match=all(r['all_match'] for r in rows if r['k']<=4),
        all_primary_large_odo_pass=all(r['passes_odo'] for r in big),
        all_primary_large_knead_pass=all(r['passes_knead'] for r in big),
        max_primary_large_relative_error=max(r['relative_error_odo'] for r in big),
        small_primary_failures=[dict(N=r['N'],k=r['k'],relative_error=r['relative_error_odo'])
            for r in rows if r['k']<=4 and r['N']>=1000 and not r['passes_odo']],
        higher_levels_pointwise_match=all(r['all_match'] for r in rows if r['k']>4))
    sensitivity,f64,f64raw=sign_sensitivity(u,phase)
    f64controls=[]
    for info in dictionaries[:4]:
        k,P,LL=info['k'],info['P'],info['word_length']
        mm,_=train(raw,truth[k-1],LL)
        aa=decode(f64raw,mm,LL,max(NS)+1)
        f64controls.append(dict(k=k,P=P,diagnostic=stat(aa,truth[k-1][:len(aa)])))
    old=json.loads((HERE.parent/'results/D4.json').read_text())
    result=dict(u=u,NS=NS,denominator='N-1, n=1..N-1, all enumerated',
        summary=summary,rows=rows,symbol_sensitivity=sensitivity,float64=f64,
        float64_dictionary_controls=f64controls,
        previous_D4_comparison=dict(old_u=old['orbit']['u'],
            old_x210=old['orbit']['x210'],old_summary=old['summary'],
            old_k4_32word_conflicts=old['cylinder_diagnostics'][3]['conflicting_cylinders'],
            note='New main decoder trains one full true period and permits C; old main trained P phases, length 32, excludes C; both parameter and decoder changes matter.'),
        elapsed_seconds=time.perf_counter()-start)
    save('d4_rerun.json',result)
    with (HERE/'REPORT.md').open('a') as f:
        f.write('\n## 柱字典与 D4 结果（追加）\n\n')
        f.write('实际以输出的 100 位参数做 120 位高精度迭代，未在周期端点重置轨道。'+str(length)+' 点的带保护符号与 210 字周期词完全一致；另以 170 位复核前 5001 点。C 定义为 |x|<1e-70；原始临界回返虽非精确零，但与最小非临界距离相隔数十个数量级：`'+json.dumps(orbitinfo,ensure_ascii=False)+'`。这是一组明确有限精度的数值检验。\n\n')
        f.write('直接符号 L/R 不等于里程计筛编码（详见 JSON 各层两个方向匹配率）。字典只用真实周期的 210 个标签训练，测试时只查观察到的前向符号词，不用时间地址补齐。逐层选择首个无训练冲突词长；字典和逐点首错在 `results/cylinders.json`。k≤4 的存活点定义为本次周期轨道中 gcd(相位,P)=1 的点，其符号柱表示由监督标签恢复，不是论文独立给定的禁止柱。\n\n')
        f.write('|k|P|最短标签可辨词长|柱数|扩展测试匹配率|首个失败|后续周期错误|\n|---|---|---|---|---|---|---|\n')
        for d in dictionaries:
            v=d['diagnostics']
            f.write(f'|{d["k"]}|{d["P"]}|{d["word_length"]}|{d["distinct_words"]}|{v["matching_rate"]:.9f}|{v["first_failure"]}|{v["heldout_mismatches"]}|\n')
        f.write('\nk≤4 字典建成，后续周期逐点匹配率 100%；完整 210 相位恢复需 '+str(L)+' 字，实际迭代测试无误。k=5…7 仅能生成一个周期内自洽的拟合字典，不能称地址存活字典建成。周期 210 的相同未来符号不能区分 n 与 n+210，而新的 11/13/17 整除标签会不同。超过直接迭代窗口的测试使用已验证周期词的符号延拓，长度 '+str(HORIZON)+'，不是新模拟的无限塔。\n\n')
        for d in dictionaries[4:]:
            f.write(f'k={d["k"]} 不可辨识见证与最优误差下限：`'+json.dumps(d['intrinsic_period_obstruction'],ensure_ascii=False)+'`。\n\n')
        f.write('D4 共 '+str(len(NS))+' 个偶数、126 组（主目标 72 组 k≤4），全量 n=1…N−1。地址 gcd 计数与独立剩余类进程计数逐组相同，完整周期密度以 Fraction 精确等于 Haar。主目标：`'+json.dumps(summary,ensure_ascii=False)+'`。所有大尺度 k≤4 的符号频率与地址版逐点相等，最大相对 Haar 误差 <1e-3。小尺度不达标项保持披露，不更换分母；高层地址结果也在 JSON 全部保留，但无有限 210 周期符号对照的合法解释。\n\n')
        f.write('|N|k|freq_odo|freq_knead|Haar|绝对误差|相对误差|符号错对数|\n|---|---|---|---|---|---|---|---|\n')
        for r in rows:
            if r['k']<=4:
                f.write(f'|{r["N"]}|{r["k"]}|{r["freq_odo"]:.10f}|{r["freq_knead"]:.10f}|{r["Haar"]:.10f}|{r["absolute_error_odo"]:.4g}|{r["relative_error_odo"]:.9g}|{r["pointwise_mismatches"]}|\n')
        f.write('\n### C 与短词对照\n\n')
        f.write('恢复精确参数本身并不足以让旧译码器全部通过；固定 32 字和排除 C 的对照如下，含 C 导致未知，短词仍可冲突。\n\n')
        f.write('|k|词长|排除C|冲突词|未知点|确定错点|首个失败|\n|---|---|---|---|---|---|---|\n')
        for c in controls[:4]:
            for v in c['variants']:
                s=v['diagnostic']
                f.write(f'|{c["k"]}|{v["word_length"]}|{v["exclude_C"]}|{v["conflicting_words"]}|{s["unknown"]}|{s["mismatches"]}|{s["first_failure"]}|\n')
        f.write('\n### 参数与数值稳定性\n\n')
        f.write('对每组直接迭代 20001 点、110 位精度。原始符号把精确 C 与任何非零值区别开；扰动离开根后在第 210 步丢失 C 是必然现象，因此另报排除临界相位后的首个 L/R 差异。\n\n')
        f.write('|参数|首次原始符号差异|非临界首次L/R差异|非临界差异总数|后10000步非临界差异|\n|---|---|---|---|---|\n')
        for s in sensitivity:
            f.write(f'|{s["label"]}|{s["first_sign_difference_raw"]}|{s["first_LR_difference_noncritical"]}|{s["LR_differences_noncritical"]}|{s["noncritical_tail_differences"]}|\n')
        f.write('\nfloat64 对照：`'+json.dumps(f64,ensure_ascii=False)+'`；使用主字典的 float64 结果见 JSON。高精度且有临界保护时周期一致性显著优于七位截断，但任何非零扰动都不能保留精确 C；稳定性结论依赖 δ 的大小与临界符号约定。\n\n')
        f.write(f'本段计算耗时 {result["elapsed_seconds"]:.3f} 秒（不含求根、绘图、文件完整性审计）。不声称证明无穷里程计共轭或哥德巴赫。\n')
    print(json.dumps(summary,indent=2),flush=True)
    print('Experiment seconds',result['elapsed_seconds'],flush=True)


if __name__=='__main__':
    main()
