import re,sys
s=open(sys.argv[1],encoding='utf-8').read()
body=re.sub(r'(?<!\\)%.*','',s)
# braces
d=0;bad=[]
for i,ch in enumerate(body):
    if ch=='{' and body[i-1]!='\\': d+=1
    elif ch=='}' and body[i-1]!='\\': d-=1
    if d<0: bad.append(i); d=0
print('brace final depth',d,'neg',len(bad))
# brackets \[ \]
print('\\[',body.count('\\['),'\\]',body.count('\\]'))
# environments
st=[];err=[]
for m in re.finditer(r'\\(begin|end)\{([^}]*)\}',body):
    if m.group(1)=='begin': st.append(m.group(2))
    else:
        if not st or st[-1]!=m.group(2): err.append((m.group(2),st[-1] if st else None,body[:m.start()].count('\n')+1))
        else: st.pop()
print('env errors',err,'unclosed',st)
# dollars per paragraph
t=body.replace('\\$','')
t=re.sub(r'\\\[.*?\\\]','',t,flags=re.S)
paras=t.split('\n\n'); odd=[]
ln=1
for p in paras:
    if p.count('$')%2: odd.append(ln)
    ln+=p.count('\n')+2
print('odd $ paragraphs at lines',odd)
labels=re.findall(r'\\label\{([^}]*)\}',body)
refs=set(re.findall(r'\\(?:ref|eqref)\{([^}]*)\}',body))
import collections
print('dup labels',[k for k,v in collections.Counter(labels).items() if v>1])
print('missing refs',sorted(refs-set(labels)))
print('non-paper labels',[l for l in labels if not l.startswith('paper:')])
print('unused labels',sorted(set(labels)-refs))
cites=set()
for c in re.findall(r'\\cite(?:\[[^]]*\])?\{([^}]*)\}',body):
    cites|={x.strip() for x in c.split(',')}
bib=re.findall(r'\\bibitem\{([^}]*)\}',body)
print('missing bib',sorted(cites-set(bib)),'uncited bib',sorted(set(bib)-cites),'dup bib',[k for k,v in collections.Counter(bib).items() if v>1])
print('lines',s.count('\n'))
