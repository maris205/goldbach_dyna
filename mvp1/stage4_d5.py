#!/usr/bin/env python3
"""Offline D5. Only exclusive-create outputs and append the prewritten report.
Run: PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -B stage4_d5.py
"""
from pathlib import Path
from fractions import Fraction
import hashlib
import json
import math
import platform
import time
import subprocess
import tempfile
import numpy as np
import scipy
from scipy.fft import rfft, irfft, next_fast_len

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'results'
REPORT = ROOT / 'REPORT_stage4.md'
PRIMES = [2,3,5,7,11,13,17,19,23,29,31,37]
X = 2_000_000
YS = [100,1000,10000,100000,1000000]
MOD, GENERATOR = 998244353, 3


def save_new(path, obj):
    with path.open('x', encoding='utf-8') as f:
        json.dump(obj, f, indent=2, ensure_ascii=False, allow_nan=False)
        f.write('\n')


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def make_ntt(size):
    assert (MOD-1) % size == 0
    rev = np.zeros(size, dtype=np.int64)
    idx = np.arange(size, dtype=np.int64)
    for _ in range(size.bit_length()-1):
        rev = (rev << 1) | (idx & 1)
        idx >>= 1
    stages = []
    length = 2
    while length <= size:
        w = pow(GENERATOR, (MOD-1)//length, MOD)
        tw = np.ones(1, dtype=np.int64)
        step = w
        while len(tw) < length//2:
            tw = np.concatenate((tw, tw * step % MOD))
            step = step * step % MOD
        inv = np.concatenate((tw[:1], MOD-tw[:0:-1])) if len(tw)>1 else tw.copy()
        stages.append((length,tw,inv))
        length *= 2
    return rev, stages


def ntt(a, rev, stages, inverse=False):
    a = a[rev].copy()
    for length, tw, inv in stages:
        blocks = a.reshape(-1,length)
        half = length//2
        u = blocks[:,:half].copy()
        v = blocks[:,half:] * (inv if inverse else tw) % MOD
        blocks[:,:half] = (u+v) % MOD
        blocks[:,half:] = (u-v) % MOD
    if inverse:
        a = a * pow(len(a), MOD-2, MOD) % MOD
    return a


def direct_counts(N):
    # Independent full-length pair sieve: kill n=0 mod p and n=N mod p.
    a = np.ones(N, dtype=bool)
    counts = []
    for p in PRIMES:
        a[::p] = False
        a[N % p::p] = False
        counts.append(int(np.count_nonzero(a)))
    return counts


def haar_num(N, primes):
    out = 1
    for p in primes:
        out *= p-1 if N % p == 0 else p-2
    return out


def exact_ratio_min(ns, counts, numerators, period):
    # Common P cancels: compare m/(N*hnum) using Python unbounded integers.
    bm, bd, bn, ties = None, None, None, 0
    for n,m,h in zip(ns.tolist(), counts.tolist(), numerators.tolist()):
        d = n*h
        if bm is None or m*bd < bm*d:
            bm,bd,bn,ties = m,d,n,1
        elif m*bd == bm*d:
            ties += 1
    f = Fraction(bm*period,bd)
    return dict(value=float(f), exact=str(f), N=bn, m=bm, ties=ties)


SWIFT = r'''
import AppKit
import Foundation
let args=CommandLine.arguments
let s=try JSONSerialization.jsonObject(with:Data(contentsOf:URL(fileURLWithPath:args[1]))) as! [String:Any]
let W=1300.0,H=760.0
let rep=NSBitmapImageRep(bitmapDataPlanes:nil,pixelsWide:Int(W),pixelsHigh:Int(H),bitsPerSample:8,samplesPerPixel:4,hasAlpha:true,isPlanar:false,colorSpaceName:.deviceRGB,bytesPerRow:0,bitsPerPixel:0)!
NSGraphicsContext.saveGraphicsState();NSGraphicsContext.current=NSGraphicsContext(bitmapImageRep:rep)
NSColor.white.setFill();NSBezierPath(rect:NSRect(x:0,y:0,width:W,height:H)).fill()
func txt(_ t:String,_ x:Double,_ y:Double,_ z:Double=18){(t as NSString).draw(at:NSPoint(x:x,y:y),withAttributes:[.font:NSFont.systemFont(ofSize:z),.foregroundColor:NSColor.black])}
let low=s["low"] as! Double,high=s["high"] as! Double
let ticks=s["ticks"] as! [String]
let series=s["series"] as! [[String:Any]]
func px(_ j:Int)->Double{return 100+Double(j)/Double(max(1,ticks.count-1))*900}
func py(_ v:Double)->Double{return 110+(v-low)/(high-low)*520}
for i in 0...5{let v=low+(high-low)*Double(i)/5;let y=py(v);NSColor.lightGray.setStroke();let p=NSBezierPath();p.move(to:NSPoint(x:100,y:y));p.line(to:NSPoint(x:1000,y:y));p.stroke();txt(String(format:"%.3g",v),20,y-8)}
for (j,t) in ticks.enumerated(){txt(t,px(j)-15,75)}
for (i,c) in series.enumerated(){let col=NSColor(calibratedHue:CGFloat(i)/CGFloat(max(1,series.count)),saturation:0.85,brightness:0.72,alpha:1);col.setStroke();col.setFill();let p=NSBezierPath();p.lineWidth=2
let ys=c["y"] as! [Double]
for (j,v) in ys.enumerated(){let point=NSPoint(x:px(j),y:py(v));if j==0{p.move(to:point)}else{p.line(to:point)};NSBezierPath(ovalIn:NSRect(x:point.x-3,y:point.y-3,width:6,height:6)).fill()}
p.stroke();let ly=625-Double(i)*36;NSBezierPath(rect:NSRect(x:1040,y:ly,width:22,height:4)).fill();txt(c["name"] as! String,1070,ly-8,16)}
txt(s["title"] as! String,100,695,24);txt(s["title"] as! String,100,695,24);txt(s["xlabel"] as! String,390,28,19)
NSGraphicsContext.current?.flushGraphics();NSGraphicsContext.restoreGraphicsState();try rep.representation(using:.png,properties:[:])!.write(to:URL(fileURLWithPath:args[2]),options:.withoutOverwriting)
'''


def plot_new(path, series, title, xticks, ylim):
    assert not path.exists()
    spec={'series':[{'name':n,'y':v} for n,v in series], 'title':title,
          'ticks':list(map(str,xticks)), 'low':float(ylim[0]), 'high':float(ylim[1]),
          'xlabel':'k' if len(xticks)==12 else 'Y (logarithmic axis; full windows [Y, 2Y])'}
    with tempfile.TemporaryDirectory(prefix='D5_plot_') as temp:
        temp=Path(temp)
        (temp/'plot.swift').write_text(SWIFT)
        (temp/'spec.json').write_text(json.dumps(spec))
        subprocess.run(['/usr/bin/swift','-module-cache-path',str(temp/'cache'),
                        str(temp/'plot.swift'),str(temp/'spec.json'),str(path)],check=True)


def main():
    start=time.perf_counter()
    assert REPORT.exists() and '运行前固定的成功标准' in REPORT.read_text()
    for name in ['D5.json','stage4_sha256_after.json','D5_max_exception.png','D5_min_ratio.png']:
        assert not (OUT/name).exists(), name
    baseline=json.loads((OUT/'stage4_sha256_before.json').read_text())
    ns=np.arange(2,X+1,2,dtype=np.int64)
    odd=np.arange(1,X,2,dtype=np.int64)
    alive=np.ones(len(odd),dtype=bool)
    size=1 << (2*len(odd)-2).bit_length()
    t=time.perf_counter();rev,stages=make_ntt(size)
    times={'ntt_setup':time.perf_counter()-t}
    print('NTT size',size,'X',X,flush=True)
    # Small unit check tests transform roots and linear convolution before full run.
    sr,ss=make_ntt(16)
    test=np.array([1,0,1,1,0,1,0,0]+[0]*8,dtype=np.int64)
    assert np.array_equal(ntt(ntt(test,sr,ss),sr,ss,True),test)
    assert np.array_equal(ntt(ntt(test,sr,ss)**2 % MOD,sr,ss,True)[:15],np.convolve(test[:8],test[:8]))
    exceptions={};growth=[];checks=[];allcounts=[]
    hnums=np.ones(len(ns),dtype=np.int64)
    period=1
    for k,p in enumerate(PRIMES,1):
        kt=time.perf_counter();part={}
        t=time.perf_counter()
        # A direct prime-by-prime Boolean sieve, compressed to odd positions.
        if p != 2:alive[odd % p == 0]=False
        part['boolean_sieve']=time.perf_counter()-t
        period*=p
        hnums*=np.where(ns % p == 0,p-1,p-2)
        t=time.perf_counter()
        a=np.zeros(size,dtype=np.int64);a[:len(alive)]=alive
        transform=ntt(a,rev,stages)
        counts=ntt(transform*transform % MOD,rev,stages,True)[:len(ns)].copy()
        # coefficient j counts N=2(j+1); max <= N/2 < MOD.
        assert np.all((counts>=0)&(counts<=ns//2))
        part['exact_integer_convolution']=time.perf_counter()-t
        t=time.perf_counter()
        fs=next_fast_len(2*len(alive)-1)
        fa=rfft(alive.astype(np.float64),n=fs)
        fc=irfft(fa*fa,n=fs)[:len(ns)]
        residual=float(np.max(np.abs(fc-np.rint(fc))))
        assert np.array_equal(np.rint(fc).astype(np.int64),counts)
        part['independent_float_convolution']=time.perf_counter()-t
        ex=ns[counts==0].tolist()
        exceptions[str(k)]={'p_k':p,'P_k':period,'exceptions':ex,'count':len(ex),'max_N':max(ex) if ex else None}
        t=time.perf_counter()
        for Y in YS:
            sel=(ns>=Y)&(ns<=2*Y)
            nn,mm,hh=ns[sel],counts[sel],hnums[sel]
            mi=int(mm.min());at=int(nn[np.flatnonzero(mm==mi)[0]])
            ratio=exact_ratio_min(nn,mm,hh,period)
            ratio['Haar_exact']=str(Fraction(haar_num(ratio['N'],PRIMES[:k]),period))
            growth.append({'k':k,'Y':Y,'interval':[Y,2*Y],'even_N_count':len(nn),
                'min_m':mi,'min_m_N':at,'min_m_ties':int(np.count_nonzero(mm==mi)),
                'min_ratio':ratio})
        part['exact_window_minima']=time.perf_counter()-t
        allcounts.append(counts)
        checks.append({'k':k,'fft_coefficients_checked':len(ns),'fft_matches_integer':True,
                       'fft_max_integer_residual':residual})
        part['total']=time.perf_counter()-kt;times[str(k)]=part
        print('k',k,'exceptions',ex,'seconds',round(part['total'],3),flush=True)
    t=time.perf_counter()
    rng=np.random.default_rng(20261102)
    samples=set(range(2,2001,2))
    for lo,hi in [(2002,10000),(10002,100000),(100002,1000000),(1000002,X)]:
        samples.update((2*rng.integers(lo//2,hi//2+1,size=64)).tolist())
    for item in exceptions.values():samples.update(item['exceptions'])
    for g in growth:samples.update([g['min_m_N'],g['min_ratio']['N']])
    for p in PRIMES:
        samples.update([2*p,2*p*p,2*p*(X//(2*p))])
    samples.update([X-2,X])
    direct_records=[]
    for N in sorted(samples):
        got=direct_counts(N)
        expected=[int(c[N//2-1]) for c in allcounts]
        assert got==expected,(N,got,expected)
        direct_records.append({'N':N,'counts_by_k':got})
    times['direct_pair_sieve_crosscheck']=time.perf_counter()-t
    # Independently enumerate small genuine prime-pair facts; these are not m_k.
    a=np.ones(201,dtype=bool);a[:2]=False
    for p in range(2,15):
        if a[p]:a[p*p::p]=False
    small=[]
    for N in range(2,101,2):
        pairs=[[n,N-n] for n in range(2,N//2+1) if a[n] and a[N-n]]
        small.append({'N':N,'unordered_prime_pairs':pairs})
    t=time.perf_counter()
    plot_new(OUT/'D5_max_exception.png',[('max N (0 = none)',[exceptions[str(k)]['max_N'] or 0 for k in range(1,13)])],
             'D5: maximum empty-set exception vs k',list(range(1,13)),(0,110))
    plot_new(OUT/'D5_min_ratio.png',[(f'k={k}',[g['min_ratio']['value'] for g in growth if g['k']==k]) for k in range(1,13)],
             'D5: normalized minimum by window',['1e2','1e3','1e4','1e5','1e6'],(0,1.05))
    times['plots']=time.perf_counter()-t
    versions={'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
              'platform':platform.platform(),'plot_backend':'Swift system AppKit'}
    cross={'all_consistent':True,'integer_modulus':MOD,'primitive_root':GENERATOR,'transform_size':size,
           'integer_coefficient_bound':X//2,'no_modular_wrap':X//2<MOD,'ntt_unit_tests_passed':True,
           'FFT_by_k':checks,'direct_pair_sieve_seed':20261102,'direct_N_count':len(samples),
           'direct_k_N_checks':12*len(samples),'direct_records':direct_records,
           'crt_used':False,'crt_note':'CRT product is not used to compute interval counts.'}
    times['experiment_before_output']=time.perf_counter()-start
    result={'convention':'positive even N=2..X; 0<n<N; ordered n count; endpoints excluded',
            'X':X,'primes':PRIMES,'k_max':12,'scanning':'all even N, no sampling for exceptions or minima',
            'exceptions_by_k':exceptions,'min_growth':growth,'crosscheck':cross,
            'small_actual_Goldbach_checks':small,'versions':versions,'timing_seconds':times,
            'plots':['results/D5_max_exception.png','results/D5_min_ratio.png'],
            'limits':['Fixed k model; 1 is admissible and primes <=p_k are removed.',
                      'No analytic exceptional-set bound or Goldbach proof is claimed.']}
    save_new(OUT/'D5.json',result)
    with REPORT.open('a',encoding='utf-8') as f:
        f.write('\n### 完成情况与方法\n\n')
        f.write(f'标准 2–6 全部达成：k=1,…,12，全部 {X//2:,} 个正偶数 N≤{X:,}；五个窗口全部枚举，未截短 [10^6,2×10^6]。文件完整性标准 1 见末尾审计。\n\n')
        f.write('主算法先逐素数在布尔数组 A(n) 中划掉倍数；A(1)=1，A(0)=0。只压缩存储奇数 A(2i+1)，其线性自卷积第 j 项恰为 m_k(2j+2)。这是两个布尔存活指示函数的离散求和，没有使用 Haar 预测计数。用模 998244353、原根 3、长度 2^21 的 NTT 整数运算计算，补零避免循环混叠；每个系数≤10^6<模数，因此恢复的是精确整数。int64 模乘最大值小于 10^18，不溢出。\n\n')
        f.write(f'独立 scipy 实数 FFT 比对全部 12,000,000 个系数，全部一致；最大距整数误差 {max(c["fft_max_integer_residual"] for c in checks):.6g}，最终计数不依赖浮点舍入。另对 {len(samples)} 个 N、{12*len(samples)} 个 (k,N) 逐 N 重新建立长度 N 的布尔数组，依次划掉 n≡0 (mod p) 与 n≡N (mod p) 两类，计数全部一致。含 N≤2000 全量、四层各64个固定种子抽样、全部例外及窗口见证、各 p|N 的小/大样本。抽样仅用于交叉检查，主扫描全量。直接复核点及各层计数全部存于 JSON。未采用 CRT 公式交叉检查。\n\n')
        f.write('Haar 以精确有理数计算：共同分母 P_k，分子逐 p 乘 p−1（p|N）或 p−2（p∤N）。每个窗口的归一化最小值由所有偶数的 Python 任意精度整数交叉乘法比较，包含并列计数，故不是浮点候选估计。min m 与 min ratio 分别优化，见证可能不同。\n\n')
        f.write('### 全部空集例外\n\n|k|p_k|例外个数|最大例外 N|全部例外 N|\n|---|---|---|---|---|\n')
        for k,e in exceptions.items():
            f.write(f'|{k}|{e["p_k"]}|{e["count"]}|{e["max_N"] if e["max_N"] is not None else "无"}|{", ".join(map(str,e["exceptions"])) or "无"}|\n')
        f.write('\n### 最小窗口增长（闭区间，两种最小值各自取最小 N 见证）\n\n|k|Y|min m|见证 N|min m/(N Haar)|比值见证 N|比值并列数|\n|---|---|---|---|---|---|---|\n')
        for g in growth:
            r=g['min_ratio']
            f.write(f'|{g["k"]}|{g["Y"]}|{g["min_m"]}|{g["min_m_N"]}|{r["value"]:.12g}|{r["N"]}|{r["ties"]}|\n')
        f.write('\nJSON 同时保存每个比值的精确分数、见证 Haar 分数、min m 的并列数及扫描个数。图：`results/D5_max_exception.png`（无例外用 0 作图标记）、`results/D5_min_ratio.png`（Y 等距对数刻度）。\n\n')
        f.write('### 数学解释与局限\n\n固定 k 时，对每个 N 的剩余类集合以 P_k 为周期，其一个完整周期存活数恰为 P_k·Haar(G_k^N)。对偶数 N，每个奇素数最多禁两类，2 只禁一类，故该数至少为 1；m_k(N)=N·Haar+O(P_k)，且 N≥P_k 时必含完整周期，因而固定层的空集例外严格有限。这里的 O 常数可取 1，并且这一界不对增长的 k 给出有效一致相对误差。实验与固定层无大例外一致，超出 X 的结论只能依据上述周期论证，不能借有限扫描外推具体最大例外。\n\n')
        f.write('“真正哥德巴赫例外来自 k(N)→∞ 时的筛误差与素数定理级误差”须理解为研究方向：固定层放过大量合数和 1；把有限层候选数转成真实素数对，需要随 N 增长的筛深度、处理小素数等边界项，并控制筛余项和素数对分布误差。只有素数定理的一点分布信息并不足以解决二素数相关性/筛的奇偶障碍，不能把这句话当成已证明的误差分解。p_k≥sqrt(N) 时存活的正数除 1 外只能是大于 p_k 的素数，但被筛掉的小素数对仍需单独补回。\n\n')
        f.write('本模型的空集与经典哥德巴赫例外不同：例如小素数本身会被本层划掉，而 1 被保留。另以直接埃氏筛独立枚举 N=2,…,100 的全部无序素数对，结果存于 JSON；此有限范围仅 N=2 无素数对（4=2+2，6=3+3），绝不据此断言任意范围的经典例外。没有声称改进任务提到的 Pintz X^0.72 或任何解析例外集上界，也没有核验/引用该文献。\n\n')
        f.write('### 环境与耗时\n\nPython/包版本：`'+json.dumps(versions,ensure_ascii=False)+'`。命令：`PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -B stage4_d5.py`。无需安装包；系统 Swift 调用 AppKit 只用于画图。\n\n|项|秒|\n|---|---|\n')
        for name,val in times.items():
            f.write(f'|{name}|{val["total"] if isinstance(val,dict) else val:.6f}|\n')
        f.write('\n每层布尔筛、整数卷积、独立 FFT、精确最小值比较的分项时间均见 D5.json。计时不含编写代码与最终输出/文件审计。未试 10^7，也未对 X 外抽样。\n')
    # Audit after experiment and report output, including all original venv files.
    t=time.perf_counter();changed=[];missing=[];external=[];after={}
    for name,old in baseline['files'].items():
        p=ROOT/name
        if not p.is_file():missing.append(name);continue
        current=sha(p);after[name]=current
        if current!=old:
            entry={'path':name,'before':old,'after':current}
            (external if name=='run_stage4.log' else changed).append(entry)
    audit={'checked':len(baseline['files']),'files':after,'changed':changed,'missing':missing,
           'external_log_changes':external,'all_protected_unchanged':not changed and not missing,
           'audit_seconds':time.perf_counter()-t,
           'new_files_at_audit':sorted(str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file() and str(p.relative_to(ROOT)) not in baseline['files'])}
    save_new(OUT/'stage4_sha256_after.json',audit)
    with REPORT.open('a',encoding='utf-8') as f:
        f.write(f'\n### 运行后 SHA-256 审计\n\n核验运行前全部 {audit["checked"]} 个文件；受保护文件全部不变：{audit["all_protected_unchanged"]}。缺失：{missing}；受保护文件变化：{changed}。外部 run_stage4.log 的变化单列：{json.dumps(external,ensure_ascii=False)}。审计耗时 {audit["audit_seconds"]:.6f} 秒。前后逐文件摘要为 results/stage4_sha256_before.json 与 results/stage4_sha256_after.json；新文件列表见后者（自身在创建时尚不在列表）。\n')
    assert audit['all_protected_unchanged'],audit
    print('DONE',json.dumps({'elapsed':time.perf_counter()-start,'direct_N':len(samples),'protected_unchanged':True}),flush=True)


if __name__=='__main__':
    main()
