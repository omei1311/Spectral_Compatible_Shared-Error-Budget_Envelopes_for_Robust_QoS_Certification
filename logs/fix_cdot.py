# -*- coding: utf-8 -*-
p = 'paper/paper_v6_2.md'
s = open(p, encoding='utf-8').read()
bad = "acts as $ℓ(u^H\\cdot)+ℓ_2(u'^H\\cdot)$ on $\\{v\\}^⊥$"
if bad not in s:
    # locate actual form
    i = s.find('acts as')
    print('form:', repr(s[i:i+70]))
    raise SystemExit(1)
good = "acts as $ℓ(u^Hz)+ℓ_2(u'^Hz)$ on inputs $z⊥v$"
s = s.replace(bad, good)
open(p, 'w', encoding='utf-8', newline='').write(s)
print('cdot fixed')
