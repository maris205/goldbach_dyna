import re,sys
s=open(sys.argv[1],encoding='utf-8').read()
# strip comments
t=re.sub(r'(?<!\\)%.*','',s)
# braces
depth=0;line=1;bad=[]
for i,ch in enumerate(t):
    if ch=='\n': line+=1
    if ch in '{}' and (i==0 or t[i-1]!='\\'):
        depth+= 1 if ch=='{' else -1
        if depth<0: bad.append(line); depth=0
print('brace final depth',depth,'negative at',bad[:5])
st=[]
for m in re.finditer(r'\\(begin|end)\{([^}]*)\}',t):
    ln=t[:m.start()].count('\n')+1
    if m.group(1)=='begin': st.append((m.group(2),ln))
    else:
        if not st or st[-1][0]!=m.group(2): print('env mismatch',m.group(2),'line',ln,'stack top',st[-1] if st else None)
        else: st.pop()
print('unclosed envs',st)
# math $ parity per paragraph (rough): count unescaped $ overall
d=len(re.findall(r'(?<!\\)\$',t)); print('dollar count',d,'even' if d%2==0 else 'ODD')
print('\\[ vs \\]',len(re.findall(r'\\\[',t)),len(re.findall(r'\\\]',t)))
labels=re.findall(r'\\label\{([^}]*)\}',t); refs=set(re.findall(r'\\(?:ref|eqref)\{([^}]*)\}',t))
print('undefined refs',refs-set(labels),'dup labels',{l for l in labels if labels.count(l)>1})
cites=set(c.strip() for g in re.findall(r'\\cite(?:\[[^\]]*\])?\{([^}]*)\}',t) for c in g.split(','))
bib=set(re.findall(r'\\bibitem\{([^}]*)\}',t)); print('undefined cites',cites-bib,'unused bib',bib-cites)
print('\\left/\\right',t.count('\\left'),t.count('\\right'),' \\bigl/\\bigr',len(re.findall(r'\\[Bb]igg?l',t)),len(re.findall(r'\\[Bb]igg?r',t)))
