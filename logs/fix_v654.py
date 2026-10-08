# -*- coding: utf-8 -*-
"""Renumber refs [10]-[14] -> [9]-[13] via placeholders (safe)."""
p = 'paper/paper_v6_5.md'
s = open(p, encoding='utf-8').read()
for old, ph in [("[10]", "@A@"), ("[11]", "@B@"), ("[12]", "@C@"),
                ("[13]", "@D@"), ("[14]", "@E@")]:
    s = s.replace(old, ph)
for ph, new in [("@A@", "[9]"), ("@B@", "[10]"), ("@C@", "[11]"),
                ("@D@", "[12]"), ("@E@", "[13]")]:
    s = s.replace(ph, new)
assert "@" not in s
open(p, 'w', encoding='utf-8', newline='').write(s)
import re
refs = re.findall(r'^\[(\d+)\]', s, flags=re.M)
print('reference numbers now:', refs)
