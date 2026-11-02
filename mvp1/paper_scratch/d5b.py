import json
from fractions import Fraction as F
d=json.load(open('results/D5.json'))
for r in d['min_growth']:
    if r['k'] in (12,) or r['Y']==1000000:
        mr=r['min_ratio']; print(r['k'],r['Y'],r['min_m'],r['min_m_N'],mr['exact'],float(F(mr['exact'])),mr['N'],mr['ties'])
print(d['exceptions_by_k']['12']['exceptions'])
print(all(n< d['primes'][int(k)-1]**2 for k,v in d['exceptions_by_k'].items() for n in v['exceptions']))
cc=d['crosscheck']; print({k:(v if not isinstance(v,list) else len(v)) for k,v in cc.items()})
# D4
d4=json.load(open('results/D4.json')); print(d4['summary']['max_relative_error_N_ge_10000'], d4['summary']['all_N_ge_10000_pass'], d4['summary']['knead_rows_resolved'], d4['orbit']['parameter_sensitivity'])
w=max((r for r in d4['rows'] if r['N']>=10000), key=lambda r:r['relative_error_odo']); print(w['N'],w['k'],w['count'],w['Haar_fraction'])
print([ (c['P'],c['training_cylinders'],c['conflicting_cylinders'],c['known_mismatch']) for c in d4['cylinder_diagnostics']])
d1=json.load(open('results/D1.json')); print(d1['status'], len(d1['failures']), [ (r['k'],round(r['min'],6),round(r['max'],6),r['within_1_percent']) for r in d1['rows']])
