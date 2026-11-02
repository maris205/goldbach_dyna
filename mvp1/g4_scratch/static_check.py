import re, sys
src = open(sys.argv[1], encoding='utf-8').read()
body = re.sub(r'(?<!\\)%.*', '', src)
# braces
depth = 0; bad = []
for i, ch in enumerate(body):
    if ch == '{' and body[i-1] != '\\': depth += 1
    elif ch == '}' and body[i-1] != '\\':
        depth -= 1
        if depth < 0: bad.append(i); depth = 0
print('brace final depth', depth, 'negatives', len(bad))
# environments
st = []; errs = []
for m in re.finditer(r'\\(begin|end)\{([^}]+)\}', body):
    if m.group(1) == 'begin': st.append(m.group(2))
    else:
        if not st or st[-1] != m.group(2): errs.append((m.group(2), st[-1] if st else None, body[:m.start()].count('\n')+1))
        else: st.pop()
print('env errors', errs, 'unclosed', st)
# math $ parity per paragraph-ish: count unescaped $
d = len(re.findall(r'(?<!\\)\$', body)); print('dollar count', d, 'even' if d % 2 == 0 else 'ODD')
print('\\[ vs \\]', body.count('\\['), body.count('\\]'))
labels = re.findall(r'\\label\{([^}]+)\}', body)
refs = re.findall(r'\\(?:ref|eqref)\{([^}]+)\}', body)
print('dup labels', [l for l in set(labels) if labels.count(l) > 1])
print('non-g4 labels', [l for l in labels if not l.startswith('g4:')])
print('undefined refs', sorted(set(refs) - set(labels)))
cites = set(c.strip() for g in re.findall(r'\\cite(?:\[[^]]*\])?\{([^}]+)\}', body) for c in g.split(','))
bibs = re.findall(r'\\bibitem\{([^}]+)\}', body)
print('undefined cites', sorted(cites - set(bibs)), 'uncited bibs', sorted(set(bibs) - cites))
defined = set(re.findall(r'\\newcommand\{\\(\w+)\}', body))
used = set(re.findall(r'\\([A-Za-z]+)', body))
known = set('''documentclass usepackage newcommand widehat mathbb mu mathfrak ensuremath mathsf par noindent fbox textbf smallskip Omega mathrm lambda mathcal
theoremstyle newtheorem renewcommand arabic proofname begin end title large author date maketitle section subsection label ref cite emph le ge ne to in
cap cup setminus bigcup prod sum int frac dfrac log sqrt sup ll gg gamma alpha beta delta varepsilon varphi xi eta nu rho omega chi pi iff Rightarrow
mid nmid times cdot cdots dots ldots infty approx sim subseteq subset supseteq bigl bigr Bigl Bigr quad qquad text item hline url bibitem appendix
equiv bmod pm leftrightarrow mapsto circ max min lim hat sharp exp neq prec langle rangle overline tilde ast Pi Delta theconj thequest theprop
leftmargin label kappa sigma tau iota zeta Gamma lfloor rfloor ell partial cdot pmod operatorname mathbf left right vert lvert rvert colon ee
the thebibliography'''.split())
print('possibly undefined macros', sorted(used - defined - known))
