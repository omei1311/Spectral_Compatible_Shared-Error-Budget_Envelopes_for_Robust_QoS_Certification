# -*- coding: utf-8 -*-
"""V6.5.2: three targeted fixes on paper_v6_5.md."""
import re

p = 'paper/paper_v6_5.md'
s = open(p, encoding='utf-8').read()
n = [0]


def sp(old, new):
    global s
    pat = r'\s+'.join(re.escape(w) for w in old.split())
    new_s, k = re.subn(pat, new.replace('\\', '\\\\'), s, count=1,
                       flags=re.S)
    assert k == 1, (k, old[:60])
    s = new_s
    n[0] += 1


# 1. Contribution 3 downgrade
sp("including bounded-nominal-norm experiments that quantify what "
   "additional physical information would buy",
   "including exploratory experiments on norm-restricted information "
   "classes")

# 2. Sec 7.7 boundary-family sentence downgrade
sp("by large-norm geometry, and a norm bound that excludes it can "
   "tighten the class optimum back toward the boundary value",
   "by large-norm geometry. A norm bound below the family norm "
   "excludes this exhibited large-norm near-worst-case construction; "
   "the resulting bracket is consistent with a potentially much "
   "tighter restricted-class optimum, though only this particular "
   "construction is excluded, not all bounded-norm constructions "
   "approaching $Q_G$")

# 3. spectral cache wording: dimension-dependent, not fixed-size
sp("or a fixed-size spectral cache ($λ,w,c$) sufficient for "
   "optimal-value queries — supporting cheap bounds",
   "or a dimension-dependent spectral cache ($λ,w,c$, of size growing "
   "with the ambient dimension) sufficient for optimal-value queries "
   "— supporting cheap bounds")
sp("or the spectral cache ($λ,w,c$) for optimal-value queries",
   "or the dimension-dependent spectral cache ($λ,w,c$) for "
   "optimal-value queries")

sp("Version 6.5.1 • English revision • 8 October 2026",
   "Version 6.5.2 • English revision • 8 October 2026")

open(p, 'w', encoding='utf-8', newline='').write(s)
print(f'paper_v6_5.md updated ({n[0]} fixes)')
