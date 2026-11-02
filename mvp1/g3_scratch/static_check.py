import re,sys
s=open(sys.argv[1],encoding='utf-8').read()
body=re.sub(r'(?<!\\)%.*','',s)
st=[]
for m in re.finditer(r'\\(begin|end)\{([^}]*)\}',body):
    if m.group(1)=='begin': st.append(m.group(2))
    else:
        if not st or st[-1]!=m.group(2): print("ENV MISMATCH at",m.start(),m.group(2),st[-1:]); 
        else: st.pop()
print("unclosed envs:",st)
b=0
for i,c in enumerate(body):
    if c=='{' and body[i-1]!='\\': b+=1
    elif c=='}' and body[i-1]!='\\': b-=1
    if b<0: print("brace underflow at",i); break
print("brace balance:",b)
labels=re.findall(r'\\label\{([^}]*)\}',body); refs=re.findall(r'\\(?:eq)?ref\{([^}]*)\}',body)
import collections
print("dup labels:",[k for k,v in collections.Counter(labels).items() if v>1])
print("undefined refs:",sorted(set(refs)-set(labels)))
print("labels w/o g3: prefix:",[l for l in labels if not l.startswith('g3:')])
cites=set(x.strip() for c in re.findall(r'\\cite(?:\[[^]]*\])?\{([^}]*)\}',body) for x in c.split(','))
bibs=set(re.findall(r'\\bibitem\{([^}]*)\}',body))
print("undefined cites:",cites-bibs,"unused bib:",bibs-cites)
# 未转义的 $ 配对（粗略：按段落计数 $ 个数奇偶）
for k,para in enumerate(re.split(r'\n\s*\n',body)):
    n=len(re.findall(r'(?<!\\)\$',para))
    if n%2: print("odd $ in paragraph",k,para[:80].replace('\n',' '))
g2=open('goldbach_G2.tex',encoding='utf-8').read(); g1=open('goldbach_props.tex',encoding='utf-8').read()
old=set(re.findall(r'\\label\{([^}]*)\}',g1+g2))
print("labels clashing with G1/G2:",set(labels)&old)
