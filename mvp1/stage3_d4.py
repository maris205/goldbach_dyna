#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""D4: exact odometer counts and an explicitly provisional symbol-cylinder decoder.

No network, no writes to pre-existing artifacts except appending this stage's report.
Run with PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -B stage3_d4.py.
"""
from pathlib import Path
from fractions import Fraction
from collections import defaultdict
import hashlib
import json
import math
import platform
import subprocess
import tempfile
import time
import numpy as np

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'results'
REPORT = ROOT / 'REPORT_stage3.md'
NS = sorted(set([100,128,256,512,1000,1024,2048,4096,8192,10000,
                 10002,12000,16384,20000,30030,60060,90090,100000]))
PRIMES = [2,3,5,7]
PERIODS = [2,6,30,210]
U = np.float64(1.4855245)
EPS = 1e-12


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(baseline):
    missing, changed = [], []
    for name, old in baseline['files'].items():
        p = ROOT / name
        if not p.is_file():
            missing.append(name)
        elif digest(p) != old:
            changed.append(dict(path=name, before=old, after=digest(p)))
    return dict(checked=len(baseline['files']), missing=missing, changed=changed,
                all_unchanged=not missing and not changed)


def save_new(name, data):
    with (OUT/name).open('x', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=2, allow_nan=False)
        f.write('\n')


def orbit(u, length):
    x = np.empty(length, dtype=np.float64)
    x[0] = 0.
    with np.errstate(over='raise', invalid='raise'):
        for i in range(length-1):
            x[i+1] = np.float64(1.) - u*x[i]*x[i]
    assert np.isfinite(x).all()
    return x


def cylinders(x, length, eps=EPS):
    # Forward L/R itinerary is a genuine finite partition cylinder. C is unknown.
    symbols = np.where(x > eps, 2, np.where(x < -eps, 1, 0)).astype(np.uint8)
    windows = np.lib.stride_tricks.sliding_window_view(symbols, length)
    keys = [row.tobytes() for row in windows[:max(NS)+1]]
    decoded, info = [], []
    for p in PERIODS:
        labels = defaultdict(set)
        for n in range(p):
            if b'\x00' not in keys[n]:
                labels[keys[n]].add(int(math.gcd(n,p) == 1))
        a = np.full(len(keys), -1, dtype=np.int8)
        for n, word in enumerate(keys):
            candidates = labels.get(word, set())
            if b'\x00' not in word and len(candidates) == 1:
                a[n] = next(iter(candidates))
        truth = np.gcd(np.arange(len(keys)),p) == 1
        known = a >= 0
        mismatch = known & (a != truth)
        heldout = np.arange(len(keys)) >= p+length
        info.append(dict(P=p, word_length=length, guard=eps,
            training_positions=[0,p-1], training_cylinders=len(labels),
            conflicting_cylinders=sum(len(v)>1 for v in labels.values()),
            unknown=int((~known).sum()), known_mismatch=int(mismatch.sum()),
            first_mismatch_indices=np.flatnonzero(mismatch)[:12].tolist(),
            heldout_start=p+length, heldout_known=int((known & heldout).sum()),
            heldout_mismatch=int((mismatch & heldout).sum()),
            invalid_C_words=sum(b'\x00' in w for w in keys)))
        decoded.append(a)
    return decoded, info


def pair_summary(a, n, N, truth, heldout_start):
    left, right = a[n], a[N-n]
    # Three-valued AND: a known false operand makes the pair false.
    good = (left == 1) & (right == 1)
    bad = (left == 0) | (right == 0)
    known = good | bad
    mismatch = known & (good != truth)
    unknown = int((~known).sum())
    low = int(good.sum())/len(n)
    high = (int(good.sum())+unknown)/len(n)
    hold = (n >= heldout_start) & (N-n >= heldout_start)
    return dict(freq_knead=low if not unknown else None,
        freq_knead_lower=low, freq_knead_upper=high,
        confirmed_survivors=int(good.sum()), unknown_pairs=unknown,
        known_pairs=int(known.sum()), pointwise_mismatches=int(mismatch.sum()),
        first_mismatches=n[mismatch][:12].tolist(),
        first_unknowns=n[~known][:12].tolist(),
        heldout_pairs=int(hold.sum()), heldout_known=int((hold & known).sum()),
        heldout_mismatches=int((hold & mismatch).sum()),
        all_match=not unknown and not mismatch.any())


SWIFT = r'''
import AppKit
import Foundation
let args=CommandLine.arguments
let s=try JSONSerialization.jsonObject(with:Data(contentsOf:URL(fileURLWithPath:args[1]))) as! [String:Any]
let W=1600.0,H=1000.0
let rep=NSBitmapImageRep(bitmapDataPlanes:nil,pixelsWide:Int(W),pixelsHigh:Int(H),bitsPerSample:8,samplesPerPixel:4,hasAlpha:true,isPlanar:false,colorSpaceName:.deviceRGB,bytesPerRow:0,bitsPerPixel:0)!
NSGraphicsContext.saveGraphicsState();NSGraphicsContext.current=NSGraphicsContext(bitmapImageRep:rep)
NSColor.white.setFill();NSBezierPath(rect:NSRect(x:0,y:0,width:W,height:H)).fill()
func txt(_ t:String,_ x:Double,_ y:Double,_ z:Double=18){(t as NSString).draw(at:NSPoint(x:x,y:y),withAttributes:[.font:NSFont.systemFont(ofSize:z),.foregroundColor:NSColor.black])}
let ymin=s["ymin"] as! Double,ymax=s["ymax"] as! Double
func px(_ x:Double)->Double{return 125+(log10(x)-2)/3*1120}
func py(_ y:Double)->Double{return 130+(y-ymin)/(ymax-ymin)*720}
for i in 0...6{let x=125+Double(i)/6*1120;NSColor.lightGray.setStroke();let p=NSBezierPath();p.move(to:NSPoint(x:x,y:130));p.line(to:NSPoint(x:x,y:850));p.stroke();txt(String(format:"%.0f",pow(10,2+Double(i)/2)),x-25,95)}
for i in 0...5{let v=ymin+(ymax-ymin)*Double(i)/5;let y=py(v);NSColor.lightGray.setStroke();let p=NSBezierPath();p.move(to:NSPoint(x:125,y:y));p.line(to:NSPoint(x:1245,y:y));p.stroke();txt(String(format:"%.4g",v),20,y-8)}
let cs=s["curves"] as! [[String:Any]]
let colors=[NSColor.systemBlue,NSColor.systemOrange,NSColor.systemGreen,NSColor.systemPurple]
for (j,c) in cs.enumerated(){let col=colors[c["color"] as! Int];col.setStroke();col.setFill();let xs=c["x"] as! [Double],ys=c["y"] as! [Double];let p=NSBezierPath();p.lineWidth=2
if c["dashed"] as! Bool{p.setLineDash([8,5],count:2,phase:0)}
for i in 0..<xs.count{let point=NSPoint(x:px(xs[i]),y:py(ys[i]));if i==0{p.move(to:point)}else{p.line(to:point)};if !(c["dashed"] as! Bool){NSBezierPath(ovalIn:NSRect(x:point.x-3,y:point.y-3,width:6,height:6)).fill()}}
p.stroke();let ly=840-Double(j)*35;NSBezierPath(rect:NSRect(x:1280,y:ly,width:25,height:4)).fill();txt(c["label"] as! String,1320,ly-8,16)}
if let target=s["target"] as? Double{NSColor.red.setStroke();let p=NSBezierPath();p.lineWidth=2;p.move(to:NSPoint(x:125,y:py(target)));p.line(to:NSPoint(x:1245,y:py(target)));p.stroke();txt("target 0.001",1260,py(target)+8,16)}
txt(s["title"] as! String,125,930,27);txt(s["ylabel"] as! String,125,875,19);txt("N (logarithmic axis)",550,45,21);txt(s["note"] as! String,125,10,16)
NSGraphicsContext.restoreGraphicsState();try rep.representation(using:.png,properties:[:])!.write(to:URL(fileURLWithPath:args[2]))
'''


def plots(rows):
    curves, errors = [], []
    for k in range(1,5):
        rr = [r for r in rows if r['k']==k]
        for key, label, dash in [('freq_odo','freq',False),('Haar','Haar',True)]:
            curves.append(dict(label=f'k={k} {label}',color=k-1,
                dashed=dash,x=NS,y=[r[key] for r in rr]))
        errors.append(dict(label=f'k={k}',color=k-1,dashed=False,x=NS,
                           y=[r['relative_error_odo'] for r in rr]))
    specs=[('D4_frequency.png',dict(title='D4: full enumeration vs Haar',
        ymin=0.,ymax=.55,ylabel='count / (N - 1)',curves=curves,
        note='Solid + dots: odometer; dashed: Haar. Logistic cylinder frequencies unresolved.')),
        ('D4_relative_error.png',dict(title='D4: odometer relative error',
        ymin=0.,ymax=1.08*max(r['relative_error_odo'] for r in rows),
        ylabel='abs(freq / Haar - 1)',curves=errors,target=.001,
        note='All n=1,...,N-1. Red line: strict relative-error target 0.001.'))]
    with tempfile.TemporaryDirectory(prefix='D4_plot_') as temp:
        temp=Path(temp);(temp/'plot.swift').write_text(SWIFT)
        for name,spec in specs:
            assert not (OUT/name).exists()
            (temp/'spec.json').write_text(json.dumps(spec))
            subprocess.run(['/usr/bin/swift','-module-cache-path',str(temp/'cache'),
                str(temp/'plot.swift'),str(temp/'spec.json'),str(OUT/name)],check=True)


def main():
    start=time.perf_counter()
    baseline=json.loads((OUT/'stage3_sha256_before.json').read_text())
    criteria=REPORT.read_bytes()
    assert b'## ' in criteria and not (ROOT/'stage3_d4.py').name in baseline['files']
    for name in ['D4.json','stage3_sha256_after.json','D4_frequency.png','D4_relative_error.png']:
        assert not (OUT/name).exists(), f'Refusing overwrite: {name}'
    before=audit(baseline)
    print('Pre-run integrity:',json.dumps(before),flush=True)
    x=orbit(U,max(NS)+421)
    decoded, cylinfo=cylinders(x,32)
    rows=[]
    for N in NS:
        n=np.arange(1,N,dtype=np.int64)
        for k,P in enumerate(PERIODS,1):
            truth=(np.gcd(n,P)==1) & (np.gcd(N-n,P)==1)
            count=int(truth.sum())
            residues=[r for r in range(P) if math.gcd(r,P)==1 and math.gcd(N-r,P)==1]
            # Count each progression in [1,N-1], independent of vector gcd.
            exact=sum((N-1-r)//P - (-r)//P for r in residues)
            haar=Fraction(1)
            for p in PRIMES[:k]:
                haar *= Fraction(p-(1 if N%p==0 else 2),p)
            assert Fraction(len(residues),P)==haar
            assert count==exact
            freq=count/(N-1); h=float(haar)
            sym=pair_summary(decoded[k-1],n,N,truth,P+32)
            sym['absolute_error_knead']=None if sym['freq_knead'] is None else abs(sym['freq_knead']-h)
            sym['relative_error_knead']=None if sym['freq_knead'] is None else abs(sym['freq_knead']/h-1)
            rows.append(dict(N=N,k=k,P=P,sample_size=N-1,method='full enumeration',
                count=count,exact_count=exact,period_survivors=len(residues),
                Haar_fraction=str(haar),Haar=h,freq_odo=freq,count_over_N=count/N,
                absolute_error_odo=abs(freq-h),relative_error_odo=abs(freq/h-1),
                count_over_N_absolute_error=abs(count/N-h),
                theoretical_absolute_bound=P/(N-1), **sym))
            assert abs(freq-h) <= P/(N-1)
    sensitivities=[]
    for length,eps in [(8,EPS),(16,EPS),(64,EPS),(420,EPS),(32,0.),(32,1e-10)]:
        dd, info=cylinders(x,length,eps)
        sensitivities.append(dict(length=length,guard=eps,layers=info,
            max_N_pairs=[pair_summary(a,np.arange(1,max(NS)),max(NS),
                (np.gcd(np.arange(1,max(NS)),P)==1)&(np.gcd(max(NS)-np.arange(1,max(NS)),P)==1),P+length)
                for a,P in zip(dd,PERIODS)]))
    # Adjacent decimal rounding endpoints: sensitivity only, not parameter recovery.
    param=[]
    signs=np.sign(x)
    for du in [-5e-8,5e-8]:
        y=orbit(np.float64(float(U)+du),len(x))
        differing=np.flatnonzero(np.sign(y)!=signs)
        param.append(dict(u=float(U)+du,max_abs=float(abs(y).max()),
            first_sign_difference=int(differing[0]) if len(differing) else None,
            sign_differences=int(len(differing))))
    orbit_info=dict(u=float(U),dtype=str(x.dtype),length=len(x),
        max_abs=float(abs(x).max()),overflow=False,all_finite=bool(np.isfinite(x).all()),
        guard=EPS,near_critical_indices=np.flatnonzero(abs(x)<=EPS)[:30].tolist(),
        near_critical_count=int((abs(x)<=EPS).sum()),min_nonzero_abs=float(abs(x[1:]).min()),
        x210=float(x[210]),
        recurrence=[dict(P=P,max_abs_return_difference=float(abs(x[P:]-x[:-P]).max()),
            sign_disagreements=int((signs[P:]!=signs[:-P]).sum()),
            compared=len(x)-P) for P in PERIODS],parameter_sensitivity=param)
    plots(rows)
    failures1000=[dict(N=r['N'],k=r['k'],error=r['relative_error_odo']) for r in rows if r['N']>=1000 and r['relative_error_odo']>=.001]
    high=[r for r in rows if r['N']>=10000]
    summary=dict(N_count=len(NS),row_count=len(rows),pair_evaluations=sum(N-1 for N in NS)*4,
        max_relative_error_N_ge_10000=max(r['relative_error_odo'] for r in high),
        all_N_ge_10000_pass=all(r['relative_error_odo']<.001 for r in high),
        failures_N_ge_1000=failures1000,
        knead_all_match=all(r['all_match'] for r in rows),
        knead_rows_resolved=sum(r['freq_knead'] is not None for r in rows))
    after=audit(baseline)
    result=dict(task='stage3 D4 only',criteria_sha256=hashlib.sha256(criteria).hexdigest(),
        code_sha256=digest(Path(__file__)),python=platform.python_version(),numpy=np.__version__,
        parameter_source='User task and local spec_stage3.txt / prompt_G3.txt; no higher precision or star-product implementation found.',
        search=dict(roots=['/Users/maris/goldbach_sd','/Users/maris/Documents','/Users/maris/Desktop','/Users/maris/Downloads','/Users/maris/app_dev'],
            extensions=['tex','md','py','txt','ipynb'],max_file_size='2M for external search',
            patterns=['1.485524','W2...W3...W5','kneading star/product','Chinese star product / superstable'],
            finding='Only task restatements locally; external hits were unrelated SymPy Wild variables. No primary paper located.'),
        symbolic_method='Provisional supervised finite-cylinder decoder from first P phases; not an independently established forbidden-cylinder family or conjugacy. Unknown means no justified classification. Pair bounds retain full denominator.',
        source_convention='G1: S_p[n]=R iff p divides n, starting n=0; R is sieve deletion. Logistic R means x>0, a distinct encoding requiring a dictionary.',
        NS=NS,orbit=orbit_info,cylinder_diagnostics=cylinfo,sensitivity=sensitivities,
        rows=rows,summary=summary,integrity_before=before,integrity_after=after,
        elapsed_seconds=time.perf_counter()-start,
        plots=['results/D4_frequency.png','results/D4_relative_error.png'])
    save_new('D4.json',result)
    lines=['\n## 实验结果（仅追加）\n',
        f'地址版主目标：{"达成" if summary["all_N_ge_10000_pass"] else "未达成"}。共 {len(NS)} 个 N、{len(rows)} 组，全量枚举；N≥10000 最大相对误差 {summary["max_relative_error_N_ge_10000"]:.9g}。符号版完整一致目标：{"达成" if summary["knead_all_match"] else "未达成"}。',
        '\n参数与来源：按用户任务及 spec_stage3.txt / prompt_G3.txt 转述，采用 u=1.4855245。重新搜索本目录 tex/md/py/txt，以及 Documents/Desktop/Downloads/app_dev 的 tex/md/py/txt/ipynb（外部文件上限 2 MB）；没有找到更高精度来源、星积实现或论文1全文。外部命中仅为无关 SymPy 变量。未联网。此数值只给到小数点后七位；±5e-8 仅作假设按四舍五入给值的敏感性范围，不能视为真实误差界。',
        '\n全部 N：'+', '.join(map(str,NS))+'。含三个 30030 倍数，无须用替代数。每组分母 N-1，无抽样。',
        f'\n临界轨道：float64，共 {len(x)} 点（含柱词前瞻），max |x_n|={orbit_info["max_abs"]:.17g}，全部有限、无溢出；|x|≤1e-12 共 {orbit_info["near_critical_count"]} 点；最小非零 |x|={orbit_info["min_nonzero_abs"]:.9g}；x_210={x[210]:.17g}。',
        '\n实际符号采用 L=x<0、R=x>0，|x|≤1e-12 记 C/未知。G1 的 S_p=RL^(p-1) 则是从 n=0 起、R 表示 p|n 的筛编码，两者不能直接等同。实现用真实前向 32 字柱词、初始 P_k 个相位的 gcd 标签尝试学习存活柱；训练后仅按柱词分类，不按 n 取模补齐。此方法是对有限柱译码的诊断尝试，不是已有理论给出的禁止柱构造。若柱词冲突或未见则未知；两条轨道以三值逻辑相与。',
        '\n频率区间为 [确定双存活数/(N-1), (确定双存活数+未知数)/(N-1)]。有未知时 freq_knead=null，报告表以 NA 表示；不能把已知子样本频率冒充全量频率。即便区间含 Haar 也不能称吻合。训练标签含地址信息，故匹配也不能独立证明共轭。\n',
        '|k|P|训练柱数|冲突柱数|单轨未知数|已知点错误|首个错误位置|\n|---|---|---|---|---|---|---|']
    for d in cylinfo:
        lines.append(f'|{PERIODS.index(d["P"])+1}|{d["P"]}|{d["training_cylinders"]}|{d["conflicting_cylinders"]}|{d["unknown"]}|{d["known_mismatch"]}|{d["first_mismatch_indices"][:3]}|')
    lines += ['\n周期检验（针对整个计算窗口，非理论否定）：',
        '|P|符号不等位置数|最大轨道回返差|\n|---|---|---|']
    for d in orbit_info['recurrence']:
        lines.append(f'|{d["P"]}|{d["sign_disagreements"]}|{d["max_abs_return_difference"]:.9g}|')
    lines += ['\n参数敏感性：'+json.dumps(param,ensure_ascii=False),
        '\n长度 8/16/64/420 与临界保护 0、1e-10 的敏感性结果见 D4.json。不能仅将失败归咎于 float64：参数截断、未验证的星积/相位模型以及有限柱译码本身均可能导致失败；没有建立真实物流映射与指定里程计的共轭。一个有限 210 周期超稳定轨道也不等于无限素数塔里程计。',
        '\n计数验证：每组向量 gcd 计数等于独立剩余类进程计数，完整周期密度以 Fraction 精确等于 Haar 乘积。count/N 另存 JSON；主表采用 count/(N-1)。对每个完整 P 周期精确、尾段误差有界，因此固定 k 时误差为 O(1/N)，并逐组核验粗界 P/(N-1)。',
        '\n### 全部结果\n',
        '|N|k|freq_odo|freq_knead|Haar|绝对误差 odo|相对误差 odo|knead 区间|未知对|确定错对|\n|---|---|---|---|---|---|---|---|---|---|']
    for r in rows:
        f='NA' if r['freq_knead'] is None else f'{r["freq_knead"]:.9f}'
        lines.append(f'|{r["N"]}|{r["k"]}|{r["freq_odo"]:.9f}|{f}|{r["Haar"]:.9f}|{r["absolute_error_odo"]:.3g}|{r["relative_error_odo"]:.3g}|[{r["freq_knead_lower"]:.6f},{r["freq_knead_upper"]:.6f}]|{r["unknown_pairs"]}|{r["pointwise_mismatches"]}|')
    lines += ['\n符号版误差：有未知的组不能给出唯一 |freq-Haar| 或 |freq/Haar-1|，JSON 为 null；频率界、逐点错误及未知位置前12例保留。',
        '\nN≥1000 未达到严格 10^-3 的组：'+json.dumps(failures1000,ensure_ascii=False),
        '\n图：[频率与 Haar](results/D4_frequency.png)；[相对误差 vs N](results/D4_relative_error.png)。图采用对数横轴；符号频率未完全可辨识时不画成实测点。',
        f'\n运行 `.venv/bin/python -B stage3_d4.py`，实验及绘图耗时 {result["elapsed_seconds"]:.3f} 秒；Python {platform.python_version()}，NumPy {np.__version__}。成功标准前缀 SHA-256：`{result["criteria_sha256"]}`。',
        f'\n已有文件核验：{after["checked"]} 个，全部未变={after["all_unchanged"]}。变化明细：'+json.dumps(after['changed'],ensure_ascii=False)+'；缺失：'+json.dumps(after['missing'],ensure_ascii=False)+ '。见 results/stage3_sha256_before.json 与 results/stage3_sha256_after.json。',
        '\n结论：地址实验核验有限周期局部密度，不能证明哥德巴赫。符号测试没有可靠恢复全部存活柱，不能宣称真实 kneading 双存活等同地址版。']
    assert REPORT.read_bytes()==criteria
    with REPORT.open('a',encoding='utf-8') as f:
        f.write('\n'.join(lines)+'\n')
    assert REPORT.read_bytes().startswith(criteria)
    final=audit(baseline)
    final.update(timestamp_unix=time.time(),criteria_sha256=result['criteria_sha256'],
                 report_prefix_unchanged=REPORT.read_bytes().startswith(criteria),
                 stage_outputs_sha256={name:digest(ROOT/name) for name in
                    ['stage3_d4.py','REPORT_stage3.md','results/D4.json','results/D4_frequency.png','results/D4_relative_error.png']})
    save_new('stage3_sha256_after.json',final)
    print(json.dumps(summary,indent=2),flush=True)
    print('Final integrity:',json.dumps(final['changed']),flush=True)


if __name__=='__main__':
    main()
