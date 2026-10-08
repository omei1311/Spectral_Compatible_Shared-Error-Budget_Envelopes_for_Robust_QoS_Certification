# -*- coding: utf-8 -*-
import re

s = open('paper/paper_v6_5.md', encoding='utf-8').read()
for m in re.finditer(r'"', s):
    i = m.start()
    print(repr(s[max(0, i-40):i+40]).replace('\\n', '|'))
