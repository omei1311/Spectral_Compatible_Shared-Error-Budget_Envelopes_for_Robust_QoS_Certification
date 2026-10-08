# -*- coding: utf-8 -*-
"""Fix build_v6_5_docx.py: names, input, and stricter litcompare table."""
p = 'paper/build_v6_5_docx.py'
s = open(p, encoding='utf-8').read()
s = s.replace('V6.4_EN.docx', 'V6.5_EN.docx')
s = s.replace("'paper_v6_4.md'", "'paper_v6_5.md'")
s = s.replace("'6.4'", "'6.5'")
old_hdr = ("['Ref', 'Original problem / assumptions', 'Original "
           "conclusion',\n               'Directly substitution yields our "
           "result?', 'Remaining steps'],")
new_hdr = ("['Ref', 'Known result', 'Retained info', 'Object optimized',\n"
           "               'Implies Eq.(6)?', 'Implies Thm 2?',\n"
           "               'beta0=alpha?', 'Structural difference'],")
assert old_hdr in s, "header not found"
s = s.replace(old_hdr, new_hdr)
open(p, 'w', encoding='utf-8', newline='').write(s)
print('builder patched')
