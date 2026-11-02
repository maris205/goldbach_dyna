import json
d=json.load(open('results/D5.json'))
print([k for k in d.keys()])
for k,v in d.items():
    if k in ('exceptions_by_k',): continue
    s=repr(v); print(k, s[:300])
