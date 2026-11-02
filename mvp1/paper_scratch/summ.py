import json,sys
def sk(o,d=0,maxd=2):
    if d>maxd: return
    if isinstance(o,dict):
        for k,v in list(o.items())[:40]:
            t=type(v).__name__
            s=repr(v)[:90] if not isinstance(v,(dict,list)) else (f"len={len(v)}")
            print("  "*d+f"{k}: {t} {s}")
            sk(v,d+1,maxd)
    elif isinstance(o,list) and o:
        sk(o[0],d+1,maxd)
for f in sys.argv[1:]:
    print("=====",f); sk(json.load(open(f)))
