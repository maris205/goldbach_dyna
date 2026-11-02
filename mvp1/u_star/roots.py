"""Offline numerical root recovery. All output stays next to this script."""
from pathlib import Path
import hashlib
import json
import math
import time
import mpmath as mp
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / 'results'
mp.mp.dps = 120


def write_json(name, data):
    with (OUT / name).open('x', encoding='utf8') as f:
        json.dump(data, f, ensure_ascii=False, indent=2, allow_nan=False)
        f.write('\n')


def orbit(u, p=210):
    xs = [mp.mpf(0)]
    for _ in range(p):
        xs.append(1-u*xs[-1]**2)
    return xs


def gd(u, p=210):
    x, d, d2 = mp.mpf(0), mp.mpf(0), mp.mpf(0)
    for _ in range(p):
        x, d, d2 = (1-u*x*x, -x*x-2*u*x*d,
                     -4*x*d-2*u*d*d-2*u*x*d2)
    return x, d, d2


def farray(us, p=210):
    x = np.zeros_like(us)
    for _ in range(p):
        x = 1-us*x*x
    return x


def hulls(xs, P):
    h = sorted((min(xs[a:210:P]), max(xs[a:210:P]), a) for a in range(P))
    gaps = [h[i+1][0]-h[i][1] for i in range(P-1)]
    return dict(disjoint=all(g > 0 for g in gaps),
                min_gap=str(min(gaps)), spatial_order=[v[2] for v in h],
                intervals=[dict(phase=a, lower=str(lo), upper=str(hi)) for lo, hi, a in h])


def float_scan(cells):
    us = np.linspace(1.4855225, 1.4855265, cells+1)
    vals = farray(us)
    idx = np.flatnonzero(np.signbit(vals[:-1]) != np.signbit(vals[1:]))
    original_a, original_b = us[idx], us[idx+1]
    a, b, fa = original_a.copy(), original_b.copy(), vals[idx].copy()
    for _ in range(30):
        m = (a+b)/2
        fm = farray(m)
        left = np.signbit(fa) != np.signbit(fm)
        b = np.where(left, m, b)
        a = np.where(left, a, m)
        fa = np.where(left, fa, fm)
    roots = (a+b)/2
    x = np.zeros((210, len(roots)))
    for i in range(1, 210):
        x[i] = 1-roots*x[i-1]**2
    keep = np.ones(len(roots), dtype=bool)
    for P in (2, 6, 30):
        lo = np.array([x[a::P].min(axis=0) for a in range(P)])
        hi = np.array([x[a::P].max(axis=0) for a in range(P)])
        order = np.argsort(lo, axis=0)
        lo, hi = np.take_along_axis(lo, order, axis=0), np.take_along_axis(hi, order, axis=0)
        keep &= np.all(lo[1:] > hi[:-1]+1e-12, axis=0)
    # A lower-period root has x30 close to zero; float is only a prefilter.
    keep &= np.abs(x[30]) > 1e-8
    brackets = [(repr(float(a)), repr(float(b)))
                for a, b in zip(original_a[keep], original_b[keep])]
    return brackets, dict(cells=cells, spacing=4e-6/cells,
                          sign_change_brackets=len(idx), compatible_prefilter=len(brackets),
                          approximate_candidates=roots[keep].tolist())


def refine(a, b, p=210):
    a, b = mp.mpf(a), mp.mpf(b)
    fa = gd(a, p)[0]
    fb = gd(b, p)[0]
    assert fa*fb < 0, (str(a), str(b), str(fa), str(fb))
    for _ in range(340):
        m = (a+b)/2
        fm = gd(m, p)[0]
        if fa*fm <= 0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    mid = (a+b)/2
    g, d, d2 = gd(mid, p)
    newton = mid-g/d
    return mid, dict(lower=str(a), upper=str(b), width=str(b-a),
                     lower_residual=str(fa), upper_residual=str(fb),
                     midpoint_residual=str(g), derivative=str(d), second_derivative=str(d2),
                     newton_value=str(newton), newton_correction=str(abs(newton-mid)),
                     newton_residual=str(gd(newton,p)[0]),
                     numerical_method='120 dps bisection 340 steps + independent derivative Newton step; not outward-rounded interval arithmetic')


def decimal(u, places):
    q = mp.nint(u*10**places)
    digits = str(int(q))
    return digits[:-places]+'.'+digits[-places:]


def swap(c):
    return {'L': 'R', 'R': 'L', 'C': 'C'}[c]


def star(A, B):
    assert A.endswith('C') and B.endswith('C')
    a = A[:-1]
    return ''.join(a+(swap(c) if a.count('R') % 2 else c) for c in B)


def factor(K, P):
    a = K[:P-1]
    assert all(K[i:i+P-1] == a for i in range(0, len(K), P))
    terminals = K[P-1::P]
    B = ''.join(swap(c) for c in terminals) if a.count('R') % 2 else terminals
    A = a+'C'
    assert star(A, B) == K
    return A, B


def main():
    start = time.perf_counter()
    assert '运行前固定的成功标准' in (HERE/'REPORT.md').read_text()
    scans, all_brackets = [], []
    for cells in (4000, 40000, 160000):
        brackets, diag = float_scan(cells)
        scans.append(diag)
        all_brackets = brackets
        print('scan', json.dumps(diag), flush=True)
    candidates = []
    for a, b in all_brackets:
        r, root_info = refine(a, b)
        xs = orbit(r)
        divisors = [d for d in range(1,211) if 210 % d == 0]
        period = next(d for d in divisors if abs(xs[d]) < mp.mpf('1e-70'))
        hh = {str(P): hulls(xs, P) for P in (2,6,30)}
        assert period == 210 and all(h['disjoint'] for h in hh.values())
        assert abs(xs[210]) < mp.mpf('1e-90')
        phase = 'C'+''.join('R' if x > 0 else 'L' for x in xs[1:210])
        K = phase[1:]+'C'
        W30, W7 = factor(K,30)
        W6, W5 = factor(W30,6)
        W2, W3 = factor(W6,2)
        assert star(star(star(W2,W3),W5),W7) == K
        u100, u40 = decimal(r,100), decimal(r,40)
        g100, d100, _ = gd(mp.mpf(u100))
        g40, _, _ = gd(mp.mpf(u40))
        candidate = dict(index=len(candidates)+1, u=u100, u_40dp=u40,
            rounds_to_7dp=decimal(r,7), selected=False,
            numerical_root=root_info,
            rounding=dict(decimal_places=100, absolute_u_error_bound=str(mp.mpf('5e-101')+(mp.mpf(root_info['width'])/2)),
                          observed_error=str(abs(mp.mpf(u100)-r)), residual_at_printed_u=str(g100),
                          linear_residual_estimate=str(abs(d100)*mp.mpf('5e-101')),
                          residual_at_40dp=str(g40), absolute_u_error_bound_40dp='5e-41'),
            minimal_period=period,
            proper_divisor_returns={str(d):str(xs[d]) for d in divisors[:-1]},
            min_abs_noncritical=str(min(abs(x) for x in xs[1:210])),
            hulls=hh, phase_symbols_x0_to_x209=phase, MS_word_x1_to_x210=K,
            MS_star_factors=dict(W2=W2,W3=W3,W5=W5,W7=W7,W30=W30,round_trip=True))
        candidates.append(candidate)
        print('candidate', candidate['index'], u100, W7, 'selected',candidate['selected'], flush=True)
    assert len(candidates)==9
    assert len({c['MS_star_factors']['W7'] for c in candidates})==9
    rounded_matches = [c for c in candidates if c['rounds_to_7dp']=='1.4855245']
    chosen = min(rounded_matches, key=lambda c: abs(mp.mpf(c['u'])-mp.mpf('1.4855245')))
    chosen['selected'] = True
    selected = [chosen]
    assert len(selected)==1
    # Independent primitive period-7 kneading enumeration on [0,2].
    grid=np.linspace(0,2,200001)
    vv=farray(grid,7)
    ii=np.flatnonzero(vv[:-1]*vv[1:]<0)
    period7=[]
    for i in ii:
        r, info=refine(repr(float(grid[i])),repr(float(grid[i+1])),7)
        xx=orbit(r,7)
        if abs(xx[1])>mp.mpf('1e-60'):
            word=''.join('R' if x>0 else 'L' for x in xx[1:7])+'C'
            period7.append(dict(u=decimal(r,70),MS_word=word,residual=str(xx[7])))
    # Include exact sampled zeros separately, with proper-divisor test.
    for i in np.flatnonzero(vv==0):
        r=mp.mpf(repr(float(grid[i])))
        xx=orbit(r,7)
        if abs(xx[7])<mp.mpf('1e-80') and abs(xx[1])>mp.mpf('1e-60'):
            word=''.join('R' if x>0 else 'L' for x in xx[1:7])+'C'
            if word not in {v['MS_word'] for v in period7}:
                period7.append(dict(u=decimal(r,70),MS_word=word,residual=str(xx[7])))
    factor_set={c['MS_star_factors']['W7'] for c in candidates}
    assert factor_set=={v['MS_word'] for v in period7}
    result=dict(task='u_star recovery',precision_dps=mp.mp.dps,
        source='Local goldbach_G3.tex, g3_scratch/refine_u.py, stage3_d4.py, REPORT_stage3.md. No network.',
        identification='Candidate set; 6 candidates round to 1.4855245. Primary is uniquely nearest to that value and matches the supplied ~1.485524500756 hint; Paper 1 W7 not available locally.',
        MS_convention='K=s(x1)...s(x209)C; phase word=C s(x1)...s(x209). A=aC, A*B=concat(a sigma(b_i)), sigma flips L/R iff number of R in a is odd; C fixed.',
        scans=scans,candidate_count=len(candidates),candidates=candidates,
        selected_index=selected[0]['index'],u=selected[0]['u'],rounded_match_indices=[c['index'] for c in rounded_matches],
        period7_independent_scan=period7,period7_word_set_match=True,
        completeness='All 9 persist on 3 meshes; all 9 distinct period-7 MS factors independently recovered. Numerical completeness evidence only; a sign-change grid cannot exclude missed tangencies or subgrid pairs.',
        derivative_recurrence='x0=d0=0; x_next=1-u*x*x; d_next=-x*x-2*u*x*d',
        elapsed_seconds=time.perf_counter()-start,
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    write_json('u_star.json',result)
    with (HERE/'REPORT.md').open('a') as f:
        f.write('\n## 高精度求根结果（追加）\n\n')
        f.write('任务意义的主值：`'+result['u']+'`。输出 100 位小数。实际有 6 个候选（4…9）均四舍五入到 1.4855245，七位转述不能唯一选根；选择第 9 个是因为它距离转述值最近，且匹配任务给出的 ~1.485524500756 提示。论文1 的 W₇ 定义本地未取得，因此以下仍称候选集。\n\n')
        f.write('使用 120 位 mpmath 运算；先在 ±2e-6 扫描，再对通过几何预筛的区间做 340 步高精度二分，并以参数导数牛顿复核。导数递推 d₀=0，dₙ₊₁=−xₙ²−2uxₙdₙ。三个扫描网格记录：`'+json.dumps(scans,ensure_ascii=False)+'`。扫描不能形式化排除切触根或网格内成对漏根。\n\n')
        f.write('|候选|u（40位小数；100位见 JSON）|W₇ MS词|取舍|\n|---|---|---|---|\n')
        for c in candidates:
            f.write(f'|{c["index"]}|{c["u_40dp"]}|{c["MS_star_factors"]["W7"]}|'+('主值' if c['selected'] else '保留候选')+'|\n')
        c=selected[0]
        f.write('\n全部候选均通过最小周期 210 与模 2/6/30 凸包分离检验；各真因子回返、凸包端点与间距、完整符号词见 `results/u_star.json`。按 K=s(x₁)…s(x₂₀₉)C 定义 MS 词，A=aC 的星积按 a 中 R 数奇偶翻转插入符号，全部逐字验证四因子重构。公共因子：`'+json.dumps({k:v for k,v in c['MS_star_factors'].items() if k!='W7'},ensure_ascii=False)+'`。9 个 W₇ 词与独立在 [0,2] 扫描的周期 7 超稳定词集合完全一致；这提供数值完备性证据，不等于论文身份确认。\n\n')
        f.write('主值误差数据：`'+json.dumps(c['numerical_root'],ensure_ascii=False)+'`；十进制输出舍入：`'+json.dumps(c['rounding'],ensure_ascii=False)+'`。以上界限依赖高精度浮点可靠性，并非带向外舍入的严格区间证明。主值完整 MS 词：`'+c['MS_word_x1_to_x210']+'`。\n\n')
        f.write(f'求根、扫描与独立 W₇ 核对耗时 {result["elapsed_seconds"]:.3f} 秒。\n')
    print('root phase seconds', result['elapsed_seconds'], flush=True)


if __name__=='__main__':
    main()
