# -*- coding: utf-8 -*-
"""V6.5.4: sync citation numbers in builder tables (Table 1 [9]->[13];
Table 12 Kim row -> [13]; add DKW [12] row) and soften DKW wording in
the paper body."""
p = 'paper/build_v6_5_docx.py'
s = open(p, encoding='utf-8').read()
n = 0

# 1. Table 1: robust RIS/beamforming row [1],[4],[9] -> [1],[4],[13]
old = "'Robust RIS/beamforming [1],[4],[9]'"
assert s.count(old) == 1
s = s.replace(old, "'Robust RIS/beamforming [1],[4],[13]'")
n += 1

# 2. Table 12: Kim row label [9] -> [13]
old = "['[9] Kim et al.',"
assert s.count(old) == 1
s = s.replace(old, "['[13] Kim et al.',")
n += 1

# 3. Table 12: insert DKW row right after the Parrott row
anchor = ("'describes the completion SET; no closed-form three-branch '\n"
          "                'extremum, no QoS lifting, no boundary envelope "
          "follows '\n                'by substitution'],")
assert s.count(anchor) == 1, s.count(anchor)
dkw_row = (
    "'describes the completion SET; no closed-form three-branch '\n"
    "                'extremum, no QoS lifting, no boundary envelope "
    "follows '\n                'by substitution'],\n"
    "               ['[12] Davis–Kahan–Weinberger',\n"
    "                'operator completions / norm-preserving extensions "
    "under prescribed blocks or actions',\n"
    "                'prescribed partial operator action',\n"
    "                'admissible extensions / extension norms',\n"
    "                'No', 'No', 'No',\n"
    "                'characterizes admissible / norm-preserving "
    "extensions; does not evaluate the specific maximum output "
    "functional under (alpha, beta, beta0), nor the resulting "
    "six-scalar QoS class infimum'],")
s = s.replace(anchor, dkw_row)
n += 1

open(p, 'w', encoding='utf-8', newline='').write(s)
print(f'builder updated ({n} edits)')

# 4. Body: soften DKW wording in Related Work
q = 'paper/paper_v6_5.md'
t = open(q, encoding='utf-8').read()
old = ("with the norm-preserving dilation theory of Davis–Kahan–"
       "Weinberger [12], which characterizes the minimal norm of "
       "extensions with prescribed action —")
new = ("with the operator-completion / norm-preserving extension theory "
       "of Davis–Kahan–Weinberger [12], which characterizes admissible "
       "extensions under prescribed blocks or actions —")
assert t.count(old.replace('\n', ' ')) or old in t
t = t.replace(old, new)
if 'norm-preserving dilation theory' in t:
    import re
    pat = re.escape("with the norm-preserving dilation") + r'.*?' + \
        re.escape("prescribed action —")
    t2, k = re.subn(pat, new, t, count=1, flags=re.S)
    assert k == 1
    t = t2
open(q, 'w', encoding='utf-8', newline='').write(t)
print('body DKW wording softened')
