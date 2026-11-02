"""Reproduce exactly the five requested Goldbach experiments."""
from pathlib import Path
import hashlib
import json
import math
import platform
import time

import numpy as np
import scipy
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parent
C2 = 0.6601618
PRIMES7 = [2, 3, 5, 7, 11, 13, 17]
ENUM_NS = [30, 32, 98, 210, 1000, 1024, 2310, 30030, 100000, 510510]
DOC = ROOT / 'STARTER.md'
MARKER = '以下结果由运行脚本追加。'
original = DOC.read_text()
criteria = original.split(MARKER)[0] + MARKER + '\n'
criteria_hash = hashlib.sha256(criteria.encode()).hexdigest()
started = time.time()

def sieve(limit):
    flags = np.ones(limit + 1, dtype=bool)
    flags[:2] = False
    for p in range(2, math.isqrt(limit) + 1):
        if flags[p]:
            flags[p*p::p] = False
    return flags

isprime = sieve(10_000_000)
primes = np.flatnonzero(isprime)

def singular(n):
    remaining = n
    s = 2 * C2
    for p0 in primes:
        p = int(p0)
        if p * p > remaining:
            break
        if remaining % p == 0:
            if p > 2:
                s *= (p - 1) / (p - 2)
            while remaining % p == 0:
                remaining //= p
    if remaining > 2:
        s *= (remaining - 1) / (remaining - 2)
    return s

def count(n):
    return int(np.count_nonzero(isprime[2:n-1] & isprime[n-2:1:-1]))

def integral(n):
    # Symmetry improves endpoint resolution; QUADPACK errors are estimates.
    value, error = quad(lambda t: 1 / (math.log(t) * math.log(n-t)),
                        2, n / 2, epsabs=1e-8, epsrel=1e-11, limit=200)
    return 2 * value, 2 * error

enumeration = []
for k in range(1, 8):
    modulus = math.prod(PRIMES7[:k])
    residues = np.arange(modulus, dtype=np.int64)
    unit = np.gcd(residues, modulus) == 1
    unit_count = int(unit.sum())
    for n in ENUM_NS:
        g = unit & unit[(n-residues) % modulus]
        alt = unit & unit[(residues-n) % modulus]
        g_count = int(g.sum())
        ratio = (g_count / modulus) / (unit_count / modulus)**2
        formula = 2.0
        for p in PRIMES7[1:k]:
            formula *= p/(p-1) if n % p == 0 else (1-2/p)/(1-1/p)**2
        error = abs(ratio-formula)
        s = singular(n)
        enumeration.append(dict(N=n, k=k, p_k=PRIMES7[k-1], P_k=modulus,
            U_count=unit_count, G_count=g_count, R=ratio, formula=formula,
            absolute_error=error, formula_assertion=bool(error < 1e-12),
            set_assertion=bool(np.array_equal(g, alt)), singular_series=s,
            R_minus_S=ratio-s, R_over_S=ratio/s))

true_counts = []
for n in [10000, 20000, 50000, 100000, 200000, 500000, 1000000]:
    i, err = integral(n)
    r, s = count(n), singular(n)
    true_counts.append(dict(N=n, r=r, I=i, I_error_estimate=err,
                            singular_series=s, HL=s*i, r_over_HL=r/(s*i)))

# Linear convolution of the prime indicator counts ordered pairs exactly,
# after rounding FFT values; independently verify every interval count below.
from scipy.fft import rfft, irfft, next_fast_len
fft_length = next_fast_len(2 * 110000 + 1)
f = rfft(isprime[:110001].astype(float), n=fft_length)
conv = irfft(f*f, n=fft_length)
interval_rows = []
for n in range(100000, 110000, 2):
    r = int(round(float(conv[n])))
    i, err = integral(n)
    s = singular(n)
    interval_rows.append(dict(N=n, r=r, I=i, I_error_estimate=err,
        singular_series=s, r_over_I=r/i, r_over_SI=r/(s*i)))
raw = np.array([x['r_over_I'] for x in interval_rows])
adjusted = np.array([x['r_over_SI'] for x in interval_rows])
ss = np.array([x['singular_series'] for x in interval_rows])
direct_counts_agree = all(x['r'] == count(x['N']) for x in interval_rows)
def stats(a):
    return dict(mean=float(a.mean()), std_population=float(a.std(ddof=0)),
                correlation_with_S=float(np.corrcoef(a, ss)[0,1]))
interval = dict(start_inclusive=100000, stop_exclusive=110000, count=5000,
    std_ddof=0, raw=stats(raw), normalized=stats(adjusted),
    fft_max_rounding_residual=float(np.max(np.abs(conv[100000:110000:2]-np.rint(conv[100000:110000:2])))),
    all_fft_counts_verified_directly=direct_counts_agree,
    references=dict(raw_std=0.77, normalized_std=0.0185, raw_correlation=0.9987),
    rows=interval_rows)
interval['reference_differences'] = dict(
    raw_std_minus_reference=float(raw.std()-0.77),
    normalized_std_minus_reference=float(adjusted.std()-0.0185),
    raw_correlation_minus_reference=interval['raw']['correlation_with_S']-0.9987)
interval['reference_analysis'] = (
    'Raw std agrees at the reference precision; normalized std and correlation do not. '
    'All 5000 ordered counts are independently verified. The stated full singular series '
    'uses every odd prime divisor, not just primes <= sqrt(N) or the first seven primes. '
    'Changing ddof from 0 to 1 multiplies std by sqrt(5000/4999), which cannot explain '
    'the discrepancy. The supplied reference computation is unavailable, so its exact '
    'cause cannot be determined; no definition or reference has been changed.')

gamma = float(np.euler_gamma)
target = 4 * math.exp(-2 * gamma)
legendre = []
for n in [10000, 100000, 1000000, 10000000]:
    ps = primes[primes <= math.isqrt(n)]
    densities = [((int(p)-1) if n % int(p) == 0 else (int(p)-2))/int(p) for p in ps]
    density = math.prod(densities)
    r = count(n)
    ratio = n*density/r
    legendre.append(dict(N=n, sqrt_floor=math.isqrt(n), k=len(ps), r=r,
        density_product=density, predicted_count=n*density, ratio=ratio,
        ratio_minus_asymptotic_constant=ratio-target,
        ratio_minus_1_2609=ratio-1.2609,
        reference_1_0988_difference=ratio-1.0988 if n == 10000000 else None))

# Multiplicative sieve for every positive even N, independently cross-checked
# with trial-factorization singular() at every requested mean-series value.
even_s = np.full(100000, 2*C2, dtype=float)
for p0 in primes[1:]:
    p = int(p0)
    if p > 200000:
        break
    even_s[p-1::p] *= (p-1)/(p-2)
series_agree = all(abs(float(even_s[n//2-1])-singular(n)) < 1e-12
                   for n in range(2, 200001, 2))
mean = float(even_s.mean())
mean_result = dict(first_N=2, last_N=200000, count=100000,
    mean=mean, reference=1.9868, theoretical_limit=2,
    mean_minus_reference=mean-1.9868, mean_minus_theory=mean-2,
    all_values_verified_by_factorization=series_agree)
mean_result['reference_analysis'] = (
    'The reference 1.9868 is not reproduced with all odd prime divisors and the stated '
    '100000 positive even inputs. Multiplicative sieve and independent factorization '
    'agree at every input. Finite cutoff of prime divisors would reduce the mean but '
    'would define a different function; the reference provenance is unknown. The '
    'rounding of C2 cannot explain a difference of 0.0131358.')

assertions = dict(formula_all=all(x['formula_assertion'] for x in enumeration),
    sets_all=all(x['set_assertion'] for x in enumeration),
    formula_max_absolute_error=max(x['absolute_error'] for x in enumeration),
    interval_counts_crosscheck=direct_counts_agree, mean_series_crosscheck=series_agree)
result = dict(metadata=dict(python=platform.python_version(), numpy=np.__version__,
    scipy=scipy.__version__, C2=C2, gamma=gamma, legendre_constant=target,
    success_criteria_sha256=criteria_hash, elapsed_seconds=time.time()-started,
    conventions='positive even N; ordered prime pairs; natural logarithms; population standard deviation'),
    assertions=assertions, step1=enumeration, step2=true_counts, step3=interval,
    step4=legendre, step5=mean_result)
(ROOT/'results').mkdir(exist_ok=True)
(ROOT/'results/starter.json').write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+'\n')

def table(headers, rows):
    return '\n'.join(['| '+' | '.join(headers)+' |', '| '+' | '.join(['---']*len(headers))+' |'] +
        ['| '+' | '.join(str(v) for v in row)+' |' for row in rows])

out = ['\n\n## 1. 剩余类枚举\n',
    'C₂ 固定为题设的 0.6601618；有限 k 的 R 不要求已等于极限 𝔖。\n',
    table(['N','k','P_k','|U_k|','|G_k^N|','R','𝔖','R−𝔖','公式误差','集合相同'],
        [[x['N'],x['k'],x['P_k'],x['U_count'],x['G_count'],f"{x['R']:.10f}",
          f"{x['singular_series']:.10f}",f"{x['R_minus_S']:.10f}",f"{x['absolute_error']:.3e}",
          x['set_assertion']] for x in sorted(enumeration,key=lambda x:(x['N'],x['k']))]),
    f"\n结论：70 组公式断言{'全部成立' if assertions['formula_all'] else '存在失败'}，最大误差 {assertions['formula_max_absolute_error']:.3e}，集合断言{'全部成立' if assertions['sets_all'] else '存在失败'}；有限素数截断解释了 R 与 𝔖 的差异。\n",
    '## 2. 有序素数对与 Hardy–Littlewood 预测\n',
    table(['N','r','I','积分误差估计','𝔖','HL','r/HL'],
        [[x['N'],x['r'],f"{x['I']:.6f}",f"{x['I_error_estimate']:.2e}",f"{x['singular_series']:.8f}",
          f"{x['HL']:.6f}",f"{x['r_over_HL']:.8f}"] for x in true_counts]),
    '\n结论：表中 r/HL 衡量每个 N 的预测偏差，有限样本偏离 1 不构成渐近猜想的反例。\n',
    '## 3. 连续偶数区间统计\n',
    table(['量','均值','总体标准差','与 𝔖 的 Pearson 相关'],
        [[name,f"{v['mean']:.9f}",f"{v['std_population']:.9f}",f"{v['correlation_with_S']:.9f}"]
          for name,v in [('r/I',interval['raw']),('r/(𝔖I)',interval['normalized'])]]),
    '\n参考值保留为：标准差 0.77 → 0.0185，原始比值与 𝔖 的相关系数 0.9987；两种归一化量与 𝔖 的相关系数均已列出。',
    f"\n5000 个 FFT 计数逐一与直接布尔筛计数核对：{direct_counts_agree}；最大舍入残差 {interval['fft_max_rounding_residual']:.3e}。\n",
    f"结论：乘入奇异级数后标准差从 {raw.std():.6f} 降至 {adjusted.std():.6f}，原始相关为 {interval['raw']['correlation_with_S']:.7f}，可直接按题设精度与参考值核对。\n",
    f"参考差异：原始标准差按两位小数吻合 0.77，但归一化标准差与 0.0185 相差 {adjusted.std()-0.0185:.9f}，原始相关与 0.9987 相差 {interval['raw']['correlation_with_S']-0.9987:.9f}，后二者未复现。",
    '\n原因核查：全部 5000 个有序计数已独立核对；𝔖 使用所有奇素因子，不能仅取前七个或 √N 以下的素因子；ddof=1 仅使标准差乘 √(5000/4999)，不足以解释差异，参考计算来源未给出，无法确定其具体原因。\n',
    '## 4. Legendre 层局部密度乘积\n',
    '每个素数的密度为 p|N 时 (p−1)/p，否则为 (p−2)/p；p=2 的密度为 1/2，包含在乘积内。\n',
    table(['N','k=π(⌊√N⌋)','r','N×密度乘积','比值','比值−4e^(−2γ)'],
        [[x['N'],x['k'],x['r'],f"{x['predicted_count']:.6f}",f"{x['ratio']:.9f}",
          f"{x['ratio_minus_asymptotic_constant']:.9f}"] for x in legendre]),
    f"\n按 γ 计算的常数为 {target:.10f}，题设参考保留为 1.2609；N=10⁷ 的比值与参考 1.0988 之差为 {legendre[-1]['reference_1_0988_difference']:.9f}。\n",
    f"结论：N=10⁷ 得 {legendre[-1]['ratio']:.7f}，局部密度独立乘积在该有限规模下仍高估真实计数，尚未接近常数 {target:.7f}。\n",
    '## 5. 奇异级数均值\n',
    table(['偶数范围','样本数','均值','参考','理论极限','均值−参考','均值−2'],
        [['2..200000',100000,f'{mean:.10f}',1.9868,2,f'{mean-1.9868:.10f}',f'{mean-2:.10f}']]),
    f"\n结论：严格按指定范围得到 {mean:.7f}，与参考 1.9868 的差为 {mean-1.9868:.7f}，与理论极限 2 的差为 {mean-2:.7f}。\n",
    '\n参考 1.9868 未复现：全部 100000 项均经乘法筛与独立因数分解交叉核对；截断素因子乘积会压低均值，但不符合题设的完整 𝔖 定义，C₂ 的小数截断也不足以解释约 0.0131 的差异，故保留原参考并将原因标为未确定。\n',
    '## 复现与差异说明\n',
    '运行命令：`.venv/bin/python starter.py`；原始 70 组枚举、5000 行区间数据、积分误差及环境版本见 `results/starter.json`。',
    '\n两种集合相同是因为 U_k 在取负下封闭；局部计数来自 CRT：p|N 排除一个剩余类，否则排除两个。',
    '\nLegendre 层的乘积是完整周期中的精确密度，但把它乘以 N 来估计区间计数需要额外的均匀性近似，且筛去的小素数本身可能出现在素数对中；固定 N 的有限筛积与大 N 渐近常数不能直接等同。',
    '\nC₂ 是题设截断小数，所有计算统一使用该值；参考数字为给定的有限精度对照，不用于拟合参数。',
    f'\n运行前判据 SHA-256：`{criteria_hash}`。',
    f'\n额外交叉核对：全部区间真实计数一致={direct_counts_agree}；全部均值序列与独立因数分解结果一致={series_agree}。\n']
if not assertions['formula_all'] or not assertions['sets_all']:
    out.append('\n断言失败明细（保留原断言）：\n' + json.dumps(
        [x for x in enumeration if not x['formula_assertion'] or not x['set_assertion']], ensure_ascii=False, indent=2))
DOC.write_text(criteria + '\n'.join(out))
assert hashlib.sha256((DOC.read_text().split(MARKER)[0]+MARKER+'\n').encode()).hexdigest() == criteria_hash
print(json.dumps(dict(assertions=assertions, interval={k:v for k,v in interval.items() if k!='rows'},
    legendre=legendre, mean=mean_result, elapsed_seconds=time.time()-started),ensure_ascii=False,indent=2))
