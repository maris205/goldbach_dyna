import re,sys
s=open(sys.argv[1],encoding='utf-8').read()
t=re.sub(r'(?<!\\)%.*','',s)
t2=t.replace('\\{','').replace('\\}','')
print('braces',t2.count('{')-t2.count('}'))
b=re.findall(r'\\(begin|end)\{([^}]*)\}',t); st=[]
ok=True
for k,e in b:
    if k=='begin': st.append(e)
    else:
        if not st or st[-1]!=e: print('mismatch',e,st[-1:]); ok=False
        else: st.pop()
print('env ok',ok and not st)
labs=re.findall(r'\\label\{([^}]*)\}',t); refs=set(re.findall(r'\\ref\{([^}]*)\}',t))
import collections
print('dup labels',[l for l,c in collections.Counter(labs).items() if c>1])
print('undef refs',refs-set(labs))
print('nonprefixed',[l for l in labs if not l.startswith('g3:')])
cites=set(c for g in re.findall(r'\\cite\{([^}]*)\}',t) for c in g.split(','))
bib=set(re.findall(r'\\bibitem\{([^}]*)\}',t)); print('undef cites',cites-bib)
# $ parity per paragraph-ish
d=re.sub(r'\\\$','',t).count('$'); print('dollar parity',d%2)
print('\\[ vs \\]',t.count('\\['),t.count('\\]'))
