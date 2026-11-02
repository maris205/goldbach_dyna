#!/usr/bin/env python3
"""Offline, create-only D2/D3 experiments. Run with .venv/bin/python stage2.py."""
from pathlib import Path
import math, time, json, hashlib, platform, subprocess, tempfile
import numpy as np
import scipy
from scipy.fft import rfft, irfft, next_fast_len
from scipy.integrate import quad
from scipy.interpolate import CubicSpline
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results'
REPORT=ROOT/'REPORT_stage2.md'
MARK=b'<!-- STAGE2 RESULTS BELOW; CRITERIA ABOVE ARE IMMUTABLE -->\n'
C2=0.6601618158468695739278121100145557784326233602847334
START=time.perf_counter()
TIMES={}
def log(s): print(s,flush=True)
def dump(name,obj):
    with (OUT/name).open('x') as f: json.dump(obj,f,ensure_ascii=False,indent=2,allow_nan=False)
def sieve(n):
    a=np.ones(n+1,dtype=bool);a[:2]=False
    for p in range(2,math.isqrt(n)+1):
        if a[p]:a[p*p::p]=False
    return a
def factor(n):
    result={}
    while n%2==0: result[2]=result.get(2,0)+1;n//=2
    p=3
    while p*p<=n:
        if n%p==0:
            result[p]=0
            while n%p==0:result[p]+=1;n//=p
        p+=2
    if n>1:result[n]=1
    return result
def singular(n):
    s=2*C2
    for p in factor(n):
        if p>2:s*=(p-1)/(p-2)
    return s
def direct_r(n,prime,ps):
    p=ps[:np.searchsorted(ps,n//2,side='right')]
    half=int(np.count_nonzero(prime[n-p]))
    return 2*half-int(n%2==0 and prime[n//2])
def integral(n):
    # Symmetry, followed by t=exp(u), controls endpoint scales.
    val,err=quad(lambda u: math.exp(u)/(u*math.log(n-math.exp(u))),math.log(2),math.log(n/2),epsabs=1e-9,epsrel=2e-12,limit=200)
    return 2*val,2*err
def independent_integral(n):
    v,e=quad(lambda t: 1/(math.log(t)*math.log(n-t)),2,n/2,epsabs=1e-9,epsrel=2e-11,limit=300)
    return 2*v,2*e

def d2(prime):
    start=time.perf_counter(); cases={};labels={}
    for m in range(4,9):
        target=10**m
        candidates=[('target',target),('power2',2**round(math.log2(target))),('multiple30030',max(1,round(target/30030))*30030),('multiple510510',max(1,round(target/510510))*510510),('ordinary',target+2468)]
        for label,n in candidates:labels.setdefault(n,[]).append(f'10^{m}:{label}')
    ps=np.flatnonzero(prime[:math.isqrt(max(labels))+1])
    checks=[]
    for n in sorted(labels):
        t=time.perf_counter();alive=np.ones(n,dtype=bool);alive[0]=False;count=n-1;h=1.;lnp=0.;rows=[];critical=None
        for k,p0 in enumerate(ps[:np.searchsorted(ps,math.isqrt(n),side='right')],1):
            p=int(p0);removed=0
            v=alive[::p];removed+=int(np.count_nonzero(v));v[:]=False
            if n%p:
                v=alive[n%p::p];removed+=int(np.count_nonzero(v));v[:]=False
            count-=removed;h*=1-(1 if n%p==0 else 2)/p;lnp+=math.log(p)
            e=count-n*h;rel=abs(e)/(n*h);x=lnp/math.log(n)
            row=dict(k=k,p=p,count=count,newly_removed=removed,Haar=h,E=e,relative_error=rel,x=x)
            rows.append(row)
            if critical is None and rel>.1:critical=dict(k=k,p=p,x=x,relative_error=rel)
            if n<=40000:
                # Full integer gcd for small layers, then independent per-prime gcd,
                # avoiding fixed-width overflow of the full primorial.
                if k<=8:
                    product=math.prod(int(q) for q in ps[:k])
                    brute=sum(math.gcd(a,product)==1 and math.gcd(n-a,product)==1 for a in range(1,n))
                else:
                    xx=np.arange(1,n,dtype=np.int64);ok=np.ones(n-1,dtype=bool)
                    for q in ps[:k]:ok &= (np.gcd(xx,q)==1)&(np.gcd(n-xx,q)==1)
                    brute=int(ok.sum())
                assert brute==count,(n,k,brute,count)
                checks.append(dict(N=n,k=k,count=count))
        assert count==int(np.count_nonzero(alive))
        cases[n]=dict(N=n,labels=labels[n],S=singular(n),critical=critical,layers=rows,seconds=time.perf_counter()-t)
        log(f'D2 N={n}: kc={critical}, {time.perf_counter()-t:.2f}s')
    result=dict(status='passed',criteria_sha256=CRITERIA,method='direct sequential boolean sieve; count only newly removed positions; 0<n<N',cases=list(cases.values()),gcd_checks=checks,elapsed_seconds=time.perf_counter()-start)
    dump('D2.json',result);TIMES['D2']=result['elapsed_seconds'];return result

def statistics(ns,rr,ss,ii,coverage):
    ratio=rr/(ss*ii);bins=[]
    for m in range(0,9):
        sel=(ns>=10**m)&(ns<10**(m+1))
        a=ratio[sel]
        if len(a):bins.append(dict(lower=10**m,upper_exclusive=10**(m+1),count=len(a),coverage=coverage,mean=float(a.mean()),variance_population=float(a.var()),std_population=float(a.std()),min=float(a.min()),max=float(a.max())))
    def row(i):
        n=int(ns[i]);return dict(N=n,r=int(rr[i]),I=float(ii[i]),S=float(ss[i]),ratio=float(ratio[i]),factorization={str(p):e for p,e in factor(n).items()})
    lo=np.argpartition(ratio,19)[:20];hi=np.argpartition(ratio,len(ratio)-20)[-20:]
    low=[row(int(i)) for i in lo[np.argsort(ratio[lo])]];high=[row(int(i)) for i in hi[np.argsort(ratio[hi])[::-1]]]
    # Separate ten million from its preceding half-open bin.
    return dict(bins=bins,minimum20=low,maximum20=high),ratio

def d3(prime):
    start=time.perf_counter();limit=10_000_000;t=time.perf_counter()
    odd=prime[1:limit:2].astype(np.float64);size=next_fast_len(2*len(odd)-1)
    spectrum=rfft(odd,n=size,workers=2);del odd
    spectrum*=spectrum
    raw=irfft(spectrum,n=size,workers=2)[:limit//2];del spectrum
    nearest=np.rint(raw);residual=float(np.max(np.abs(raw-nearest)))
    assert residual<1e-6,residual
    rr=nearest.astype(np.int64);del raw,nearest
    rr[1]=1 # N=4: (2,2); all other even N use odd+odd.
    ns=np.arange(6,limit+1,2,dtype=np.int64);rr=rr[2:]
    assert np.all(rr>=0)
    assert np.all(rr%2==prime[ns//2])
    TIMES['D3_fft']=time.perf_counter()-t;log(f'D3 FFT {TIMES["D3_fft"]:.2f}s residual={residual:g}')
    t=time.perf_counter();ps=np.flatnonzero(prime).astype(np.int32)
    test=np.unique(np.r_[np.arange(6,2001,2),np.geomspace(2002,limit,150).astype(np.int64)//2*2,2**np.arange(3,24)])
    test=test[test<=limit];checks=[]
    for n0 in test:
        n=int(n0);actual=direct_r(n,prime,ps);fft=int(rr[(n-6)//2]);assert fft==actual,(n,fft,actual)
        checks.append(dict(N=n,FFT=fft,direct=actual))
    TIMES['D3_count_validation']=time.perf_counter()-t
    t=time.perf_counter()
    knots=np.unique(np.r_[np.arange(6,1002,2),np.exp(np.linspace(math.log(1000),math.log(1e8),2200)),1e8])
    logknots=np.unique(np.log(knots));knots=np.exp(logknots)
    vals=np.array([integral(float(n))[0] for n in knots]);interp=CubicSpline(logknots,np.log(vals))
    rng=np.random.default_rng(20261102)
    testI=np.unique(np.r_[np.arange(6,2002,2),np.geomspace(2002,1e8,350).astype(np.int64)//2*2,rng.integers(1001,50_000_001,500)*2])
    integral_checks=[]
    for n0 in testI:
        n=int(n0);ref,err=independent_integral(n);estimate=float(np.exp(interp(math.log(n))));rel=abs(estimate/ref-1)
        assert rel<1e-6,(n,rel)
        integral_checks.append(dict(N=n,reference=ref,interpolated=estimate,relative_error=rel,quad_relative_error_estimate=err/ref))
    ii=np.exp(interp(np.log(ns)))
    TIMES['D3_integral']=time.perf_counter()-t;log(f'D3 integral validation max={max(x["relative_error"] for x in integral_checks):g}')
    t=time.perf_counter();ss=np.full(limit//2,2*C2)
    for p0 in ps[1:np.searchsorted(ps,limit//2,side='right')]:
        p=int(p0);ss[p-1::p]*=(p-1)/(p-2)
    ss=ss[2:]
    serr=0.
    for n0 in np.unique(np.r_[np.arange(6,10002,2),test]):
        n=int(n0);err=abs(ss[(n-6)//2]-singular(n));serr=max(serr,err);assert err<1e-12
    TIMES['D3_singular']=time.perf_counter()-t
    full,rat=statistics(ns,rr,ss,ii,'all even N')
    # Reproducible random sample and a separately identified targeted stratum.
    t=time.perf_counter();randoms=[]
    for lower in range(10_000_000,100_000_000,10_000_000):
        randoms.extend((2*rng.choice(5_000_000,300,replace=False)+(lower+2)).tolist())
    randomset=set(randoms);special=set(range(30030*((10_000_000//30030)+1),100_000_001,30030))
    special.update(2**k for k in range(24,27));special.add(100_000_000)
    special.update(c['N'] for c in D2['cases'] if 10_000_000<c['N']<=100_000_000)
    hn=np.array(sorted(randomset|special),dtype=np.int64);hr=np.empty(len(hn),dtype=np.int64)
    for j,n0 in enumerate(hn):
        hr[j]=direct_r(int(n0),prime,ps)
        if j%500==0:log(f'D3 high samples {j}/{len(hn)} elapsed {time.perf_counter()-t:.1f}s')
    hs=np.array([singular(int(n)) for n in hn]);hi=np.exp(interp(np.log(hn)))
    high,hrat=statistics(hn,hr,hs,hi,'random + targeted sample, biased toward 30030 multiples')
    israndom=np.array([int(n) in randomset for n in hn]);randomstats,_=statistics(hn[israndom],hr[israndom],hs[israndom],hi[israndom],'stratified uniform random, 300 per subinterval')
    sample_rows=[dict(N=int(n),r=int(r),S=float(s),I=float(i),ratio=float(a),random=int(n) in randomset,targeted=int(n) in special) for n,r,s,i,a in zip(hn,hr,hs,hi,hrat)]
    allmin=sorted(full['minimum20']+high['minimum20'],key=lambda x:x['ratio'])[:20]
    allmax=sorted(full['maximum20']+high['maximum20'],key=lambda x:x['ratio'],reverse=True)[:20]
    powers=[]
    for k in range(3,27):
        n=2**k
        if n<=limit:j=(n-6)//2;powers.append(dict(N=n,ratio=float(rat[j]),S=float(ss[j])))
        elif n in hn:j=int(np.searchsorted(hn,n));powers.append(dict(N=n,ratio=float(hrat[j]),S=float(hs[j])))
    # A small-N effect is separated from the asymptotic range.
    tail=ns>=10000;ids=np.flatnonzero(tail);tailrat=rat[tail];idx=ids[np.argpartition(tailrat,19)[:20]]
    tailmin=[]
    for j in idx[np.argsort(rat[idx])]:
        n=int(ns[j]);tailmin.append(dict(N=n,r=int(rr[j]),S=float(ss[j]),I=float(ii[j]),ratio=float(rat[j]),factorization={str(p):e for p,e in factor(n).items()}))
    TIMES['D3_high_sampling']=time.perf_counter()-t
    hist_edges=np.linspace(.90,1.10,241);plotbins=[]
    for m in range(3,7):
        v=rat[(ns>=10**m)&(ns<10**(m+1))];counts,edges=np.histogram(v,hist_edges)
        plotbins.append(dict(label=f'10^{m} <= N < 10^{m+1}',x=((edges[1:]+edges[:-1])/2).tolist(),y=(counts/(len(v)*np.diff(edges))).tolist(),outside=int(np.count_nonzero((v<hist_edges[0])|(v>hist_edges[-1])))))
    counts,edges=np.histogram(hrat[israndom],hist_edges);plotbins.append(dict(label='10^7 < N <= 10^8 random',x=((edges[1:]+edges[:-1])/2).tolist(),y=(counts/(israndom.sum()*np.diff(edges))).tolist(),outside=int(np.count_nonzero((hrat[israndom]<.9)|(hrat[israndom]>1.1)))))
    concentration=dict(S_threshold=1.5,full_mean_S=float(ss.mean()),full_fraction_S_le_1_5=float(np.mean(ss<=1.5)),tail_minimum20_fraction_S_le_1_5=sum(r['S']<=1.5 for r in tailmin)/20,tail_minimum20_mean_S=float(np.mean([r['S'] for r in tailmin])),tail_minimum20_power2_count=sum((r['N']&(r['N']-1))==0 for r in tailmin))
    result=dict(status='passed',criteria_sha256=CRITERIA,C2=C2,concentration=concentration,range_full=dict(first=6,last=limit,count=len(ns)),excluded=[dict(N=2,reason='empty integration domain'),dict(N=4,reason='zero integral, ratio undefined')],FFT=dict(length=size,max_rounding_residual=residual,parity_check=True,nonnegative_check=True,direct_checks=checks,method='odd prime indicator linear convolution; no cyclic wrap; float64 rounding checked'),integral=dict(method='symmetry + log substitution quad, log-log cubic spline',knots=len(knots),max_validation_relative_error=max(x['relative_error'] for x in integral_checks),checks=integral_checks),singular_max_absolute_crosscheck_error=serr,full=full,high=high,high_random=randomstats,high_samples=sample_rows,minimum20=allmin,maximum20=allmax,minimum20_N_at_least_10000=tailmin,powers_of_two=powers,distribution=dict(range=[.9,1.1],curves=plotbins),timings=TIMES.copy(),elapsed_seconds=time.perf_counter()-start)
    dump('D3.json',result);TIMES['D3_total']=result['elapsed_seconds'];return result

# Native offline plotting keeps the Python environment untouched.
SWIFT=r'''
import AppKit
import Foundation
let args=CommandLine.arguments
let data=try Data(contentsOf: URL(fileURLWithPath:args[1]))
let spec=try JSONSerialization.jsonObject(with:data) as! [String:Any]
let W=1800.0,H=1100.0
let rep=NSBitmapImageRep(bitmapDataPlanes:nil,pixelsWide:Int(W),pixelsHigh:Int(H),bitsPerSample:8,samplesPerPixel:4,hasAlpha:true,isPlanar:false,colorSpaceName:.deviceRGB,bytesPerRow:0,bitsPerPixel:0)!
NSGraphicsContext.saveGraphicsState();NSGraphicsContext.current=NSGraphicsContext(bitmapImageRep:rep)
NSColor.white.setFill();NSBezierPath(rect:NSRect(x:0,y:0,width:W,height:H)).fill()
func text(_ s:String,_ x:Double,_ y:Double,_ size:Double=19){(s as NSString).draw(at:NSPoint(x:x,y:y),withAttributes:[.font:NSFont.systemFont(ofSize:size),.foregroundColor:NSColor.black])}
let xmin=spec["xmin"] as! Double,xmax=spec["xmax"] as! Double,ymin=spec["ymin"] as! Double,ymax=spec["ymax"] as! Double
let logy=spec["logy"] as! Bool
func px(_ x:Double)->Double{return 120+(x-xmin)/(xmax-xmin)*1180}
func py(_ y:Double)->Double{let v=logy ? log10(max(y,pow(10,ymin))) : y;return 130+(v-ymin)/(ymax-ymin)*840}
for i in 0...6 {let v=xmin+(xmax-xmin)*Double(i)/6;let x=px(v);NSColor.lightGray.setStroke();let p=NSBezierPath();p.move(to:NSPoint(x:x,y:130));p.line(to:NSPoint(x:x,y:970));p.stroke();text(String(format:"%.3g",v),x-20,95)}
for i in 0...6 {let v=ymin+(ymax-ymin)*Double(i)/6;let y=130+840*Double(i)/6;NSColor.lightGray.setStroke();let p=NSBezierPath();p.move(to:NSPoint(x:120,y:y));p.line(to:NSPoint(x:1300,y:y));p.stroke();text(logy ? String(format:"1e%.1f",v) : String(format:"%.2g",v),18,y-8)}
let curves=spec["curves"] as! [[String:Any]]
for (j,c) in curves.enumerated(){let col=NSColor(calibratedHue:Double(j%25)/Double(max(5,curves.count)),saturation:0.75,brightness:0.72,alpha:1);col.setStroke();let xs=c["x"] as! [Double],ys=c["y"] as! [Double];let path=NSBezierPath();path.lineWidth=1.6;var started=false
for (x,y) in zip(xs,ys){if x<xmin || x>xmax || y<0 {continue};let yy=py(y);if !started {path.move(to:NSPoint(x:px(x),y:yy));started=true}else{path.line(to:NSPoint(x:px(x),y:yy))}}
path.stroke();let ly=960-Double(j)*30;col.setFill();NSBezierPath(rect:NSRect(x:1340,y:ly+3,width:25,height:5)).fill();text(c["label"] as! String,1378,ly-5,17)
if let cx=c["criticalX"] as? Double,let cy=c["criticalY"] as? Double{col.setFill();NSBezierPath(ovalIn:NSRect(x:px(cx)-5,y:py(cy)-5,width:10,height:10)).fill()}}
if logy {NSColor.red.setStroke();let p=NSBezierPath();p.lineWidth=2;p.move(to:NSPoint(x:120,y:py(0.1)));p.line(to:NSPoint(x:1300,y:py(0.1)) );p.stroke()}
text(spec["title"] as! String,120,1030,29);text(spec["xlabel"] as! String,450,48,23);text(spec["ylabel"] as! String,120,990,20)
text(spec["note"] as! String,120,8,16)
NSGraphicsContext.restoreGraphicsState();let png=rep.representation(using:.png,properties:[:])!;try png.write(to:URL(fileURLWithPath:args[2]))
'''
def plots(d2,d3):
    t=time.perf_counter();curves=[]
    for c in d2['cases']:
        rows=c['layers'];v=dict(label=f'N={c["N"]}',x=[r['x'] for r in rows],y=[max(r['relative_error'],1e-8) for r in rows])
        if c['critical']:v.update(criticalX=c['critical']['x'],criticalY=c['critical']['relative_error'])
        curves.append(v)
    specs=[('D2_error.png',dict(title='D2: window error along sieve layers',xlabel='x = ln(P_k) / ln(N)',ylabel='|count - N Haar| / (N Haar), log scale',note='Red line: 10%; dots: first strict crossing; zero errors clipped to 1e-8.',xmin=0.,xmax=max(c['x'][-1] for c in curves),ymin=-8.,ymax=0.,logy=True,curves=curves)),('D3_distribution.png',dict(title='D3: normalized Goldbach count distributions',xlabel='ratio = r(N) / (S(N) I(N))',ylabel='Probability density',note='Displayed ratio range 0.9 to 1.1; omitted tails counted in D3.json. High range: random sample only.',xmin=.9,xmax=1.1,ymin=0.,ymax=1.12*max(max(c['y']) for c in d3['distribution']['curves']),logy=False,curves=d3['distribution']['curves']))]
    with tempfile.TemporaryDirectory(prefix='goldbach_stage2_') as temp:
        path=Path(temp);(path/'plot.swift').write_text(SWIFT)
        for name,spec in specs:
            (path/'spec.json').write_text(json.dumps(spec));subprocess.run(['/usr/bin/swift','-module-cache-path',str(path/'cache'),str(path/'plot.swift'),str(path/'spec.json'),str(OUT/name)],check=True)
    TIMES['plots']=time.perf_counter()-t

def table(headers,rows):return '\n'.join(['|'+'|'.join(headers)+'|','|'+'|'.join(['---']*len(headers))+'|']+['|'+'|'.join(map(str,row))+'|' for row in rows])
def report(d2,d3,protected):
    content=['\n## 实验结果\n',f'成功标准全部达成；总耗时 {time.perf_counter()-START:.2f} 秒。运行前判据 SHA-256：`{CRITERIA}`。Python {platform.python_version()}，NumPy {np.__version__}，SciPy {scipy.__version__}。',f'已有文件 {len(protected)} 个全部 SHA-256 未改变。仅本地计算，无联网、commit 或 push。', '\n### D2 临界层\n',table(['N','𝔖(N)','标签','k_c','p_c','x_c'],[(c['N'],f'{c["S"]:.8f}',','.join(c['labels']),c['critical']['k'] if c['critical'] else '未超过',c['critical']['p'] if c['critical'] else '—',f'{c["critical"]["x"]:.6f}' if c['critical'] else '—') for c in d2['cases']]),f'共 {len(d2["cases"])} 个不同 N，逐层 gcd 独立复核 {len(d2["gcd_checks"])} 组；所有层最终存活数另以 count_nonzero 核验。小尺度最近正 510510 倍数为 510510，最近正 30030 倍数为 30030，因此两者超出对应目标量级，未冒充量级内样本。2 的幂按对数距离选最近。','D2 过筛会删除小素数本身，故末层计数不是完整 r(N)。P_k 仅累计对数，不构造巨大整数（小值 gcd 验证除外）。有限窗口误差首次过阈值不要求后续单调。\n','### D3 分段统计\n']
    bins=d3['full']['bins']+d3['high_random']['bins']
    content.append(table(['区间','个数','覆盖','均值','总体方差'],[(f'[{b["lower"]},{b["upper_exclusive"]})',b['count'],b['coverage'],f'{b["mean"]:.10f}',f'{b["variance_population"]:.10g}') for b in bins]))
    content.extend([f'全量偶数 {d3["range_full"]["count"]} 个（6..10^7）；高段去重样本 {len(d3["high_samples"])} 个，随机标记 {sum(r["random"] for r in d3["high_samples"])} 个，另加 30030 倍数、2 的幂和 D2 样本。随机与定向标记可重叠。10^7 作为全量单点落入 [10^7,10^8)，不与随机段混作同一总体；10^8 单点另列。完整定向混合分段见 D3.json high。',table(['端点 N','ratio','样本性质'],[(r['N'],f'{r["ratio"]:.10f}','定向端点') for r in d3['high_samples'] if r['N']==100000000]),'\n### 合并已计算集合的最小 20\n'])
    def extreme_table(rows):return table(['N','r(N)','𝔖(N)','ratio','分解'],[(r['N'],r['r'],f'{r["S"]:.10f}',f'{r["ratio"]:.10f}',' × '.join(f'{p}^{e}' if e>1 else p for p,e in r['factorization'].items())) for r in rows])
    content.extend([extreme_table(d3['minimum20']),'\n### 合并已计算集合的最大 20\n',extreme_table(d3['maximum20']),'\n### 去除小值效应：全量 N≥10000 的最小 20\n',extreme_table(d3['minimum20_N_at_least_10000'])])
    lo=d3['minimum20'];tail=d3['minimum20_N_at_least_10000']
    content.extend([f'最小 20 中 2 的幂有 {sum((r["N"]&(r["N"]-1))==0 for r in lo)} 个，𝔖 最小值 2C₂ 的项有 {sum(abs(r["S"]-2*C2)<1e-12 for r in lo)} 个；平均 𝔖={np.mean([r["S"] for r in lo]):.6f}。N≥10000 的最小 20 中 2 的幂有 {sum((r["N"]&(r["N"]-1))==0 for r in tail)} 个，平均 𝔖={np.mean([r["S"] for r in tail]):.6f}。这些统计不足以说明低比值集中在 2 的幂；最小比值亦有明显小 N 波动，不能只看未分段的极值。全部 2 的幂结果见 JSON powers_of_two。', '\n### 数值验证和耗时\n',f'FFT 长度 {d3["FFT"]["length"]}，最大舍入残差 {d3["FFT"]["max_rounding_residual"]:.3g}，非负与对角项奇偶校验全部通过，直接计数抽查 {len(d3["FFT"]["direct_checks"])} 个无差异。浮点 FFT 结合全量舍入残差、奇偶性和直接抽查提供数值精确性证据，并非形式化逐项整数算术证明。',f'积分 {d3["integral"]["knots"]} 个节点，独立原变量 quad 验证 {len(d3["integral"]["checks"])} 点，最大相对误差 {d3["integral"]["max_validation_relative_error"]:.3g} < 10^-6。这是抽查精度证据，不是对所有插值点的解析误差界。𝔖 交叉核验最大绝对差 {d3["singular_max_absolute_crosscheck_error"]:.3g}。',table(['部分','秒'],[(k,f'{v:.3f}') for k,v in TIMES.items()]),'D3_total 包含各 D3 子项，不重复相加；共享素数筛与画图另列。','\n### 产物与局限\n','- `results/D2.json`：全部逐层表、临界层、gcd 检验。\n- `results/D3.json`：分段统计、分范围与合并的双向极值、积分/计数核验、高段全部样本及标签。\n- `results/D2_error.png`：所有 N 一条线，10% 红线及首次跨越标记；零误差以 10^-8 显示。\n- `results/D3_distribution.png`：分量级密度，图外尾部计数见 JSON。绘图由本机 AppKit 生成，无需安装依赖。','全量原始 4999998 个比值不逐项灌入 JSON，只存统计、极值和核验信息，可用 stage2.py 重现。10^7 以上不是全量，所报极值只对已计算集合有效；含全部 30030 倍数的混合样本有定向偏差，分段主表使用随机子集。小 N 的积分模型偏差较大；N=2、4 无定义比值。C₂ 使用给定高精度常数的 float64 值，不以有限素数乘积代替。','任何结果均为有限实验，不构成哥德巴赫猜想或渐近误差界的证明。'])
    with REPORT.open('a') as f:f.write('\n\n'.join(content)+'\n')

if __name__=='__main__':
    for filename in ['D3.json','D2_error.png','D3_distribution.png','stage2_manifest.json']:
        if (OUT/filename).exists():raise SystemExit(f'Create-only: already exists {filename}')
    base=REPORT.read_bytes();assert base.endswith(MARK), 'report already has results'
    CRITERIA=hashlib.sha256(base).hexdigest()
    protected={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.rglob('*') if p.is_file() and '.venv' not in p.parts and p.name not in ['stage2.py','REPORT_stage2.md'] and p.suffix not in ['.log','.err'] and p.name not in ['G2.json']}
    t=time.perf_counter();prime=sieve(100_000_000);TIMES['shared_prime_sieve']=time.perf_counter()-t
    D2=json.loads((OUT/'D2.json').read_text()) if (OUT/'D2.json').exists() else d2(prime)
    assert D2['criteria_sha256']==CRITERIA
    TIMES['D2']=D2['elapsed_seconds']
    D3=d3(prime);plots(D2,D3)
    for name,h in protected.items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h,name
    assert hashlib.sha256(REPORT.read_bytes()).hexdigest()==CRITERIA
    report(D2,D3,protected)
    dump('stage2_manifest.json',dict(criteria_sha256=CRITERIA,protected_existing_files=protected,all_unchanged=True,timings=TIMES,total_seconds=time.perf_counter()-START))
    log(f'COMPLETE {time.perf_counter()-START:.2f}s')
