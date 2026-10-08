# -*- coding: utf-8 -*-
"""Remove script-introduced quotation artifacts from paper_v6_5.md, then
apply the V6.5.1 claim-downgrade patches on the cleaned text."""
import re

p = 'paper/paper_v6_5.md'
s = open(p, encoding='utf-8').read()
n_art = len(re.findall(r'"\n"', s))
s = re.sub(r'"\n"', '\n', s)          # join-wrapped artifact pairs
# stray lone boundary quotes left after pair removal
s = re.sub(r'(?<=[^\s])"[ \t]*\n', '\n', s)   # closing quote at EOL
s = re.sub(r'\n[ \t]*"(?=[^\n]*\n)', '\n', s)  # opening quote at BOL
left = len(re.findall(r'"', s))
print(f'artifacts removed: {n_art}; total quotes left: {left} '
      '(expected 4 = two reference titles)')
open(p, 'w', encoding='utf-8', newline='').write(s)
