from pathlib import Path
import time, math, json, hashlib, platform
from fractions import Fraction as F
import numpy as np
from scipy.fft import rfft, irfft, next_fast_len
from scipy.integrate import quad
ROOT=Path(__file__).resolve().parent
REPORT=ROOT/'REPORT.md'
marker=b'<!-- RESULTS APPEND BELOW; CRITERIA ABOVE ARE IMMUTABLE -->\n'
baseline=REPORT.read_bytes().split(marker,1)[0]+marker
digest=hashlib.sha256(baseline).hexdigest()
C2=0.6601618158468695739278121100145557784326233602847334133194484233
PS=[2,3,5,7,11,13,17,19,23,29,31,37]
start=time.perf_counter()
(ROOT/'results').mkdir(exist_ok=True)
summary=[]
def save(name,data,t,status='通过'):
    extra={
        'V2':('k=1..6; every N mod P_k','Fraction(P_k*direct_G_count,phi(P_k)^2) 与逐素数有理数乘积比较'),
        'V3':('k=1..6; complete periods','以 V1 的直接集合基数作全周期整数求和，再用 Fraction 精确计算均值'),
        'V4':('k=1..6 finite layers; X=10^3..10^8 integer trend','直接枚举基数的平方和算有限层精确矩；V6 完整奇异级数数组计算整数均值'),
        'V8':('cutoffs 10^1..10^6; 10 sample N at each cutoff','有限素数乘积与先去除2的奇异级数因数分解结果比较，检查适用条件及容差界'),
        'V9':('N=10^4,10^5,10^6,10^7,10^8','完整素数布尔筛直接计算有序对；局部乘积及对称数值积分作比较')}
    if name in extra:
        data.update(range=extra[name][0],method=extra[name][1])
    if 'failures' not in data:
        data['failures']=[] if status!='失败' else ['independent S crosscheck failed']
    data.update(status=status,elapsed_seconds=time.perf_counter()-t,criteria_sha256=digest)
    (ROOT/'results'/f'{name}.json').write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    summary.append((name,status,data['elapsed_seconds']))
    print(name,status,round(data['elapsed_seconds'],3),flush=True)
    return data
def rat(x): return {'numerator':x.numerator,'denominator':x.denominator,'float':float(x)}
def singular(n):
    while n%2==0: n//=2
    s=2*C2; p=3
    while p*p<=n:
        if n%p==0:
            s*=(p-1)/(p-2)
            while n%p==0:n//=p
        p+=2
    if n>1:s*=(n-1)/(n-2)
    return s
def sieve(n):
    a=np.ones(n+1,dtype=bool);a[:2]=False
    for p in range(2,math.isqrt(n)+1):
        if a[p]:a[p*p::p]=False
    return a
def fsum(a):
    return math.fsum(float(np.sum(a[i:i+8192],dtype=np.float64)) for i in range(0,len(a),8192))
def prodformula(n,ps):
    v=F(1)
    for p in ps: v*=F(p,p-1) if n%p==0 else F(p*(p-2),(p-1)**2)
    return v
# V1: direct complete-period arrays, every shift N, no probabilistic sampling.
t=time.perf_counter();enum=[]; failures=[];counter=None
for k in range(1,7):
    ps=PS[:k];P=math.prod(ps);x=np.arange(P);u=np.gcd(x,P)==1; neg=(-x)%P
    counts=np.empty(P,dtype=np.int64); failcounts={a:0 for a in ['G_plus','Pi_translate','Pi_reflect','center_symmetry']}
    for n in range(P):
        g=u & u[(n-x)%P];gp=u & np.roll(u,n);pi=u & np.roll(u,-n)
        checks={'G_plus':np.array_equal(g,gp),'Pi_translate':np.array_equal(pi,np.roll(g,-n)),
                'Pi_reflect':np.array_equal(pi,g[neg])}
        if n%2==0:
            center=np.roll(g,-n//2); checks['center_symmetry']=np.array_equal(center,center[neg])
        for label,ok in checks.items():
            if not ok:
                failcounts[label]+=1
                if len(failures)<20:failures.append({'k':k,'N':n,'identity':label})
        counts[n]=np.count_nonzero(g)
        if k==2 and n==2: counter={'k':2,'N':2,'P':6,'U':np.flatnonzero(u).tolist(),'G':np.flatnonzero(g).tolist(),'Pi':np.flatnonzero(pi).tolist(),'G_equals_Pi':bool(np.array_equal(g,pi))}
    enum.append((k,P,int(u.sum()),counts));print('V1 k',k,'done',flush=True)
    if any(failcounts.values()): failures.append({'k':k,'failure_counts':failcounts})
ok=not failures and counter['G']==[1] and counter['Pi']==[5] and not counter['G_equals_Pi']
v1=save('V1',{'method':'直接逐 N、逐 x 布尔集合枚举完整周期；np.roll(a,s)[x]=a[x−s]','ranges':[{'k':k,'P':P,'N_count':P,'ordered_N_x_checks':P*P} for k,P,_,_ in enum],'counterexample':counter,'failures':failures},t,'通过' if ok else '失败')
t=time.perf_counter();rows=[];fails=[]
for k,P,phi,g in enum:
    for n,c in enumerate(g):
        actual=F(P*int(c),phi**2); expect=prodformula(n,PS[:k])
        if actual!=expect:fails.append({'k':k,'N':n,'actual':rat(actual),'expected':rat(expect)})
    rows.append({'k':k,'P':P,'phi':phi,'tested_N':P,'odd_zero':bool(np.all(g[1::2]==0)),'N_zero_ratio':rat(F(P*int(g[0]),phi**2)),'N_zero_expected':rat(F(P,phi))})
v2=save('V2',{'rows':rows,'failures':fails},t,'通过' if not fails else '失败')
t=time.perf_counter();rows=[];fails=[]
for k,P,phi,g in enum:
    full=F(sum(map(int,g)),phi**2);even=2*full
    rows.append({'k':k,'full_average':rat(full),'even_average':rat(even)})
    if full!=1 or even!=2:fails.append(k)
v3=save('V3',{'rows':rows,'failures':fails},t,'通过' if not fails else '失败')
# D1 independent local residue enumeration and integer CRT cardinality recurrence.
t=time.perf_counter();ns=np.arange(2,10001,2,dtype=np.int64);ss=np.array([singular(int(n)) for n in ns]);rows=[];allrat=[];fails=[]
counts=[1]*5000;P=phi=1
for k,p in enumerate(PS[:9],1):
    local=[sum(a!=0 and (r-a)%p!=0 for a in range(p)) for r in range(p)]
    P*=p;phi*=p-1; values=[]
    for i,n0 in enumerate(ns):
        n=int(n0);counts[i]*=local[n%p];actual=F(P*counts[i],phi**2);expect=prodformula(n,PS[:k])
        if actual!=expect:fails.append({'k':k,'N':n})
        values.append(float(actual)/ss[i])
    a=np.array(values);e=np.abs(a-1)
    rows.append({'k':k,'p_k':p,'P_k':P,'tested_N':5000,'min':float(a.min()),'q01_q05_q25_q50_q75_q95_q99':np.quantile(a,[.01,.05,.25,.5,.75,.95,.99]).tolist(),'max':float(a.max()),'mean':float(a.mean()),'std_population':float(a.std()),'mean_abs_error_from_1':float(e.mean()),'max_abs_error_from_1':float(e.max()),'within_1_percent':float(np.mean(e<=.01)),'within_5_percent':float(np.mean(e<=.05))})
    allrat.append({'k':k,'R_over_S':values})
d1=save('D1',{'method':'显式局部 a=0..p−1 枚举；CRT 基数从 c 递推至 c×local_count；精确 Fraction 比较独立乘积','N_values':ns.tolist(),'S_values':ss.tolist(),'S_factorization':'先完全去除因子 2，再分解奇数部分；未运行 starter.py','rows':rows,'distributions':allrat,'failures':fails,'largest_array_length':5000},t,'通过' if not fails and P==223092870 else '失败')
# Shared prime sieve, includes X=1e8 for V9 and X/2 odd prime factors for V6.
t=time.perf_counter();isp=sieve(100_000_000);primes=np.flatnonzero(isp);sieve_time=time.perf_counter()-t
print('shared prime sieve',sieve_time,flush=True)
# V5: exact finite-window certificates plus integer convolution and boundary counts.
t=time.perf_counter();M=1_000_000
spf=np.arange(M+1,dtype=np.int32)
for p0 in primes[primes<=1000]:
    p=int(p0);view=spf[p*p::p];mask=view==np.arange(p*p,M+1,p,dtype=np.int32);view[mask]=p
m=np.arange(2,M+1,dtype=np.int64);prime=isp[2:M+1];composite_certificate=bool(np.all(spf[2:][~prime].astype(np.int64)**2<=m[~prime]));prime_certificate=bool(np.all(spf[2:][prime]==m[prime]))
L=next_fast_len(2*M+1);f=rfft(isp[:M+1].astype(float),n=L);conv=irfft(f*f,n=L);ev=np.arange(4,M+1,2);raw=conv[ev];residual=float(np.max(np.abs(raw-np.rint(raw))));r=np.rint(raw).astype(np.int64)
boundary=np.zeros(M//2+1,dtype=np.int64)
for p0 in primes[primes<=1000]:
    p=int(p0);lo=max(4,p*p);lo+=(lo%2)
    n=np.arange(lo,M+1,2,dtype=np.int64)
    boundary[n//2]+=2*isp[n-p].astype(np.int64)
boundary[2]-=1 # N=4: (2,2) belongs to both small-coordinate sets.
b=boundary[ev//2];middle=r-b;root=np.sqrt(ev);pk=np.searchsorted(primes,root,side='right');delta=r-middle
checks={'composites_spf_squared_le_m':composite_certificate,'primes_spf_equal_m':prime_certificate,'fft_rounding_residual_below_1e-6':residual<1e-6,'middle_nonnegative':bool(np.all(middle>=0)),'boundary_bounds_all':bool(np.all((delta>=0)&(delta<=2*pk)))}
samples=np.unique(np.r_[np.arange(4,2001,2),np.linspace(4,M,100,dtype=int)//2*2]);samplefails=[]
for n0 in samples:
    n=int(n0);x=np.arange(2,n-1);pair=isp[x]&isp[n-x];direct=int(pair.sum());lo=math.isqrt(n)+1;hi=n-lo
    xx=np.arange(lo,hi+1);survive=np.ones(len(xx),dtype=bool)
    for p0 in primes[primes<=math.isqrt(n)]:
        p=int(p0);survive &= (xx%p!=0)&((n-xx)%p!=0)
    c=int(survive.sum());idx=n//2-2
    if direct!=int(r[idx]) or c!=int(middle[idx]):samplefails.append({'N':n,'direct_r':direct,'fft_r':int(r[idx]),'direct_sieve_middle':c,'certified_middle':int(middle[idx])})
checks['direct_crosschecks']=not samplefails
v5=save('V5',{'N_first':4,'N_last':M,'even_N_count':len(ev),'method':'全部整数 m=2..10^6 的最小素因子证书 + 线性 FFT 素数对计数 + 逐小素数累加边界有序对；直接窗口筛交叉核对样本','certificate_quantification':'每个窗口坐标 m 满足 sqrt(N)<m<N−sqrt(N)，故 m<N 且 N<m²。若 m 为合数，spf(m)²≤m<N，所以被 p≤sqrt(N) 筛掉；若 m 为素数，m>sqrt(N)，故不被筛掉。证书穷尽每个 m，逻辑覆盖所有偶 N 和两个坐标，无周期平均替代。','boundary_method':'每个 p≤sqrt(N) 若 N−p 为素数，加 2；唯一双重覆盖 N=4,(2,2) 减 1。中段=完整素数对−边界，集合恒等式由证书逐点保证。','checks':checks,'fft_length':L,'fft_max_rounding_residual':residual,'direct_sample_N':samples.tolist(),'failures':samplefails,'min_middle':int(middle.min()),'max_boundary':int(b.max()),'max_boundary_over_2pi':float(np.max(b/(2*pk))),'examples':[{'N':int(ev[i]),'r':int(r[i]),'middle':int(middle[i]),'boundary':int(b[i]),'bound_2pi':int(2*pk[i])} for i in [0,1,48,4998,len(ev)-1]],'shared_sieve_seconds':sieve_time},t,'通过' if all(checks.values()) else '失败')
del f,conv,raw,spf,m,prime,root
# V7 sampled exact rational error, vectorized integer count prefixes.
t=time.perf_counter();rows=[];fails=[];P=phi=1;gc=np.ones(500000,dtype=np.int64);ne=np.arange(2,1000001,2,dtype=np.int64)
base=set(range(2,1001))|set(map(int,np.linspace(1001,1000000,2000)))
for k,p in enumerate(PS,1):
    P*=p;phi*=p-1;gc*=np.where(ne%p==0,p-1,p-2);pref=np.cumsum(gc,dtype=np.int64)
    samples=set(base)
    for j in [1,2,3,10]:
        for d in [-1,0,1]:
            z=j*P+d
            if 2<=z<=1000000:samples.add(z)
    maxnum=None;minnum=None;argmax=None
    for X in sorted(samples):
        num=X*phi*phi-P*int(pref[X//2-1]);den=phi*phi
        if num<0 or num*phi>=P*den:fails.append({'k':k,'X':X,'theta':rat(F(num,den))})
        if maxnum is None or num>maxnum:maxnum=num;argmax=X
        if minnum is None or num<minnum:minnum=num
    rows.append({'k':k,'p_k':p,'P_k':P,'samples':len(samples),'X_values':sorted(samples),'min_theta':rat(F(minnum,phi*phi)),'max_sample_theta':rat(F(maxnum,phi*phi)),'argmax_sample_X':argmax,'upper_bound':rat(F(P,phi)),'max_over_bound':float(F(maxnum,P*phi))})
v7=save('V7',{'rows':rows,'method':'所有 N≤10^6 的局部允许类数精确 int64 乘积和前缀和；样本 X 用 Python 大整数交叉相乘比较严格界','failures':fails},t,'通过' if not fails else '失败')
# V8 truncation ratios and admissibility.
t=time.perf_counter();rows=[];fails=[]
for limit in [10,100,1000,10000,100000,1000000]:
    ps=primes[primes<=limit];pk=int(ps[-1]);c=float(np.prod(1-1/(ps[1:].astype(float)-1)**2));bound=1/(1-1/(pk-1))
    for n in [2,4,6,30,210,2310,30030,10000,9998,510510]:
        rem=n
        while rem%2==0:rem//=2
        factors=[];p=3
        while p*p<=rem:
            if rem%p==0:
                factors.append(p)
                while rem%p==0:rem//=p
            p+=2
        if rem>1:factors.append(rem)
        ratio=2*c*math.prod((p-1)/(p-2) for p in factors if p<=pk)/singular(n);valid=pk>=max(factors,default=2)
        good=not valid or 1-1e-12<=ratio<=bound+1e-12
        row={'N':n,'cutoff':limit,'p_k':pk,'max_odd_prime_factor':max(factors,default=1),'applicable':valid,'R_over_S':ratio,'upper_bound':bound,'passed_if_applicable':good if valid else None};rows.append(row)
        if not good:fails.append(row)
v8=save('V8',{'rows':rows,'failures':fails,'float_tolerance':1e-12},t,'通过' if not fails else '失败')
# V9 direct count to 1e8 (boolean slicing, no FFT at 1e8).
t=time.perf_counter();rows=[];target=4*math.exp(-2*float(np.euler_gamma))
for n in [10000,100000,1000000,10000000,100000000]:
    ts=time.perf_counter();rr=int(np.count_nonzero(isp[2:n-1]&isp[n-2:1:-1]));ps=primes[primes<=math.isqrt(n)];mu=float(np.prod(1-1/ps.astype(float)));S=singular(n);R=2*float(np.prod(1-1/(ps[1:].astype(float)-1)**2))
    for p0 in ps[1:]:
        p=int(p0)
        if n%p==0:R*=(p-1)/(p-2)
    I,err=quad(lambda x:1/(math.log(x)*math.log(n-x)),2,n/2,epsabs=1e-7,epsrel=1e-11,limit=300);I*=2;err*=2
    main=S*n/math.log(n)**2;integ=S*I;naive=n*mu*mu*R
    rows.append({'N':n,'r_ordered':rr,'Haar_U':mu,'Haar_U_squared_logN_squared':mu*mu*math.log(n)**2,'target_4exp_minus2gamma':target,'HL_main':main,'HL_integral':integ,'I':I,'I_error_estimate':err,'N_Haar_G':naive,'r_over_HL_main':rr/main,'r_over_HL_integral':rr/integ,'N_Haar_G_over_r':naive/rr,'overestimate_percent':100*(naive/rr-1),'row_elapsed_seconds':time.perf_counter()-ts})
v9=save('V9',{'rows':rows,'integral_convention':'有效端点为 [2,N−2]，通过对称性计算 2∫[2,N/2]；原文 ∫[2,N] 在 N−t=1 处奇异，不能当作普通收敛积分','interpretation':'只展示有限样本与常数的差异；不声称实际 r 的渐近成立，也不强求有限样本高估恰为26%。','shared_sieve_seconds':sieve_time},t,'只给趋势')
# V6 full singular series for even inputs <=1e8, no cutoff of odd factors.
t=time.perf_counter();s=np.full(50_000_000,2*C2,dtype=np.float64)
for p0 in primes[1:]:
    p=int(p0)
    if p>50_000_000:break
    s[p-1::p]*=(p-1)/(p-2)
print('V6 multiplicative sieve done',time.perf_counter()-t,flush=True)
err=float(np.max(np.abs(s[:5000]-ss)));v6rows=[];v4trend=[]
for X in [1000,10000,100000,1000000,10000000,100000000]:
    ts=time.perf_counter();a=s[:X//2];total=fsum(a);E=total-X
    # alternate block length for a cheap independent summation-order sensitivity check
    alt=math.fsum(float(np.sum(a[i:i+65536],dtype=np.float64)) for i in range(0,len(a),65536))
    second=math.fsum(float(np.dot(a[i:i+8192],a[i:i+8192])) for i in range(0,len(a),8192))/(X//2)
    v6rows.append({'X':X,'count_even':X//2,'sum_S':total,'E':E,'E_over_logX':E/math.log(X),'abs_E_over_logX':abs(E)/math.log(X),'alternate_block_sum_difference':alt-total,'row_reduction_seconds':time.perf_counter()-ts})
    v4trend.append({'X':X,'even_S_squared_mean':second})
# spot exact decimal high-precision factorization to gauge aggregate float error.
import decimal
with decimal.localcontext() as ctx:
    ctx.prec=60;dc=decimal.Decimal('0.6601618158468695739278121100145557784326233602847334133194484233');decerr=[]
    for n in [2,6,30,210,2310,30030,510510,9699690,100000000]:
        rem=n
        while rem%2==0:rem//=2
        val=2*dc;p=3
        while p*p<=rem:
            if rem%p==0:
                val*=decimal.Decimal(p-1)/decimal.Decimal(p-2)
                while rem%p==0:rem//=p
            p+=2
        if rem>1:val*=decimal.Decimal(rem-1)/decimal.Decimal(rem-2)
        decerr.append({'N':n,'decimal_S':str(val),'float_minus_decimal':str(decimal.Decimal(float(s[n//2-1]))-val)})
v6=save('V6',{'rows':v6rows,'max_X_completed':100000000,'method':'50000000 个偶数上的 float64 乘法筛；每个奇素数 p≤50000000 乘到索引 p−1,p−1+p,...，覆盖全部奇素因子；8192 项分块和后 math.fsum','C2_decimal':'0.6601618158468695739278121100145557784326233602847334133194484233','S_factorization_first_strip_2':True,'crosscheck_all_5000_even_N_max_abs_error':err,'crosscheck_passed':err<=1e-12,'decimal_spot_checks':decerr,'numerical_note':'arm64 numpy.longdouble 等同 float64。报告不同分块归约差及60位十进制样本误差；不是严谨的全和误差界。E 的数量级远大于浮点归约敏感度。','shared_sieve_seconds':sieve_time},t,'只给趋势' if err<=1e-12 else '失败')
# V4 exact finite moments + trend only infinite integer average.
t=time.perf_counter();rows=[];fails=[]
for k,P,phi,g in enum:
    actual=F(P*sum(int(c)**2 for c in g),phi**4);expect=F(2)
    for p in PS[1:k]:expect*=1+F(1,(p-1)**3)
    rows.append({'k':k,'actual_second_moment':rat(actual),'expected':rat(expect),'exact_equal':actual==expect})
    if actual!=expect:fails.append(k)
ps=primes[(primes>=3)&(primes<=1000000)];constant=4*float(np.prod(1+1/(ps.astype(float)-1)**3));last=int(ps[-1]);log_tail_bound=1/(2*(last-1)**2)
for row in v4trend:row.update(target_truncated_at_1e6=constant,mean_minus_target=row['even_S_squared_mean']-constant)
v4=save('V4',{'finite_layer_rows':rows,'finite_layer_status':'通过' if not fails else '失败','integer_average_trend':v4trend,'target_product':constant,'product_max_prime':last,'omitted_log_tail_upper_bound':log_tail_bound,'tail_bound_method':'log∏(1+a)≤Σa；以所有整数尾部替代素数，用积分 1/(2(p_k−1)²) 作上界','failures':fails,'trend_shared_V6_seconds':v6['elapsed_seconds'],'interpretation':'有限二阶矩严格核对；整数 S² 平均仅为数值趋势，不判断猜测。'},t,'只给趋势' if not fails else '失败')
# Append only: criteria bytes above must be unchanged.
assert REPORT.read_bytes().startswith(baseline)
out=['\n\n## 本次运行结果\n',f'判据 SHA-256：`{digest}`。Python {platform.python_version()}，NumPy {np.__version__}；C₂ 使用预定常数。共享 10⁸ 素数筛耗时 {sieve_time:.3f} 秒。各项耗时不重复计共享筛。\n','|项目|状态|耗时（秒）|','|---|---|---|']
for name,status,dt in sorted(summary,key=lambda z:(z[0]=='D1',z[0])):out.append(f'|{name}|{status}|{dt:.3f}|')
out+=['\nV1：全部 '+str(sum(P for _,P,_,_ in enum))+' 个 N 层位、'+str(sum(P*P for _,P,_,_ in enum))+' 个 (N,x) 位置，方向反例：'+json.dumps(counter,ensure_ascii=False)+'。失败明细见 JSON。',f'\nV2、V3：全部周期精确有理数比较；V2 失败 {len(v2["failures"])}，V3 失败 {len(v3["failures"])}。',f'\nV4：有限层状态 {v4["finite_layer_status"]}；整数二阶矩只给趋势，对照截断常数 {constant:.12f}，对数尾部上界 {log_tail_bound:.3g}。趋势计算复用 V6 的完整奇异级数数组。','\n|X|S² 偶数均值|与截断常数之差|','|---|---|---|']
for row in v4trend:out.append(f'|{row["X"]}|{row["even_S_squared_mean"]:.12f}|{row["mean_minus_target"]:.8g}|')
out += [f'\nV5：499999 个偶数全部通过={all(checks.values())}；FFT 最大舍入残差 {residual:.3g}，直接样本 {len(samples) if False else len(v5["direct_sample_N"])} 个，最大边界差 {v5["max_boundary"]}。所有 m≤10⁶ 的最小素因子证书与窗口不等式逐点推出中段集合恒等式；计算有序素数对和小坐标边界后相减得到中段计数，并在样本上直接筛中段核对。没有逐 N 展开长度 N 的全数组。完整证书逻辑见 V5.json。',f'\nV6：完成到 X=10⁸；完整奇素因子乘法筛与先去除 2 的独立分解在全部 5000 个偶数上最大差 {err:.3g}；不对大 O 或 −1/2 极限作结论。','\n|X|E(X)|E/log X|不同分块和之差|','|---|---|---|---|']
for row in v6rows:out.append(f'|{row["X"]}|{row["E"]:.12f}|{row["E_over_logX"]:.9f}|{row["alternate_block_sum_difference"]:.3g}|')
out+=['\nV7：采样点均采用精确整数比较；以下最大值仅是所列样本的最大值。','\n|k|样本数|最大 θ|最大 θ/上界|','|---|---|---|---|']
for row in v7['rows']:out.append(f'|{row["k"]}|{row["samples"]}|{row["max_sample_theta"]["float"]:.9f}|{row["max_over_bound"]:.9f}|')
out+=[f'\nV8：{sum(x["applicable"] for x in v8["rows"])} 个满足最大素因子条件的样本，失败 {len(v8["failures"])}；其余条件不满足的样本明确标为不适用。最大截断素数 {last}。','\nV9：积分按 [2,N−2] 计算。原文的上限 N 会跨越 log(N−t)=0，不是普通收敛积分；因此没有直接数值求该字面表达。Mertens 常数为 '+f'{target:.12f}。有限样本偏差如下，不把“约26%”当作已经达到。','\n|N|r(N)|Haar(U)² log²N|r/HL主项|r/HL积分|N Haar(G)/r|高估%|','|---|---|---|---|---|---|---|']
for row in v9['rows']:out.append(f'|{row["N"]}|{row["r_ordered"]}|{row["Haar_U_squared_logN_squared"]:.9f}|{row["r_over_HL_main"]:.9f}|{row["r_over_HL_integral"]:.9f}|{row["N_Haar_G_over_r"]:.9f}|{row["overestimate_percent"]:.4f}|')
out+=['\nD1：45000 组精确 CRT 局部枚举核验，失败 '+str(len(d1['failures']))+'；P₉=223092870，最大分布数组仅 5000 项。所有 5000 项的 R/S 分布随 k 的数据见 D1.json。分位数使用 NumPy 线性插值，总体标准差 ddof=0。','\n|k|均值|标准差|最小|中位|最大|平均绝对偏差|最大绝对偏差|1%内比例|','|---|---|---|---|---|---|---|---|---|']
for row in d1['rows']:out.append(f'|{row["k"]}|{row["mean"]:.9f}|{row["std_population"]:.9f}|{row["min"]:.9f}|{row["q01_q05_q25_q50_q75_q95_q99"][3]:.9f}|{row["max"]:.9f}|{row["mean_abs_error_from_1"]:.9f}|{row["max_abs_error_from_1"]:.9f}|{row["within_1_percent"]:.4f}|')
out += [f'\n运行命令：`.venv/bin/python validate_stage.py`。总耗时 {time.perf_counter()-start:.3f} 秒。脚本、报告和 results/V1.json … V9.json、D1.json 为本阶段输出；未运行 starter.py，未修改 tex，未访问 g1_scratch 内容。', '\n浮点限制：V6 的 E 是两个约 10⁸ 的量相减，归约变化和十进制抽样检查保留在 JSON 中；这些检查不能替代严格的全局舍入误差证明。V4、V6、V9 的有限数值不构成渐近证明。']
with REPORT.open('wb') as h:h.write(baseline+('\n'.join(out)+'\n').encode())
assert REPORT.read_bytes().startswith(baseline)
print('DONE',time.perf_counter()-start,flush=True)
