# -*- coding: utf-8 -*-
"""Patch build_v6_1_docx.py table captions/rows for V6.1 (run once)."""
p = 'paper/build_v6_1_docx.py'
s = open(p, encoding='utf-8').read()
orig = s

s = s.replace(
    "table(doc, 'Table 6. Realistic-link pool (seeds 301–305, distance-dependent path loss, '\n"
    "                   'Rician/Rayleigh components; 540 units, 256 cells). ρ is the relative radius.',",
    "table(doc, 'Table 6. Illustrative heterogeneous-link pool (seeds 301–305; distance-'\n"
    "                   'dependent path loss, Rician K=3 and Rayleigh components; 540 units, '\n"
    "                   '256 cells). r = ρ‖ĥ‖; 20log10(ρ) is a deterministic worst-case '\n"
    "                   'energy-ratio reading, not a statistical NMSE.',")

s = s.replace("Confirmation entry: 202 is the common final accepted set; per-method '\n"
              "                   'entries are reference-call totals.",
              "Separate-gain = comparison envelope; '\n"
              "                   'confirmation entry: 202 is the common final accepted set; per-method '\n"
              "                   'entries are reference-call totals.")

s = s.replace(
    "table(doc, 'Table 8. Controlled scaling: dimension varies, all else fixed '\n"
    "                   '(4 users, per-BS power 1, γ = 1, ρ = 0.15; medians over three seeds; '\n"
    "                   'ms per user problem; pipeline seconds per seed).',",
    "table(doc, 'Table 8. Controlled scaling: dimension varies, all else fixed '\n"
    "                   '(4 users, per-BS power 1, γ = 1, ρ = 0.15; medians over three seeds). '\n"
    "                   'Per-user-problem ms exclude summary construction; the reference reuses '\n"
    "                   'each user cached eigendecomposition; pipeline seconds include summary '\n"
    "                   'construction. Four dimension points do not establish a complexity order.',")

s = s.replace("['U ≥ q* (witness validity), max deficit', '1.3e-12 (fp noise)'],",
              "['U vs reference (witness validity)',\n"
              "                'agrees within tolerance; largest deficit 1.3e-12 (fp noise)'],")

s = s.replace("['Optimal recovery [8]', 'partial', 'finite data',\n"
              "                'worst-case radius', '—', '—', 'yes (framework)', '—'],",
              "['Optimal recovery [8]', 'partial', 'finite data',\n"
              "                'worst-case radius', '—', '—', 'radius concept (related)', '—'],")

s = s.replace("['This work', 'partial', 'six scalars',\n"
              "                'safe + sharp envelope', 'yes (L2)', 'yes', 'yes', 'yes']],",
              "['This work (six scalars)', 'partial', 'six scalars',\n"
              "                'safe + sharp envelope', 'yes', 'yes', 'yes (nondegenerate)', 'yes'],\n"
              "               ['This work (6+κ branch)', '+κ', 'seven scalars',\n"
              "                'quadratic-form bound', 'yes', 'no', 'not claimed', 'no']],")

lit = (
    "    elif name == 'litcompare':\n"
    "        table(doc, 'Table 12. Formula-level comparison: what is directly used, what is "
    "context, '\n"
    "                   'what is evaluated here.',\n"
    "              ['Ref', 'Original problem / assumptions', 'Original conclusion',\n"
    "               'Direct substitution yields our result?', 'Remaining steps'],\n"
    "              [['[5] Buzano; Bottazzi–Conde',\n"
    "                'inner-product bounds for |<Tx,y>| type quantities',\n"
    "                'two-sided C-S refinements',\n"
    "                'No: cross-term bound (4) proved by our 2x2 determinant argument',\n"
    "                'conceptual kinship only'],\n"
    "               ['[6] Grone et al.',\n"
    "                'partial Hermitian matrix, chordal specified graph',\n"
    "                'PSD completion iff specified blocks PSD',\n"
    "                'No: our 2x2 completion verified directly',\n"
    "                'general theorem in background'],\n"
    "               ['[7] Parrott',\n"
    "                'contraction completion / quotient norm, lifting',\n"
    "                'operator-ball structure of completions',\n"
    "                'No: lifting theorem not applied',\n"
    "                'completion philosophy cited as context'],\n"
    "               ['[8] Micchelli–Rivlin',\n"
    "                'optimal recovery from finite data',\n"
    "                'radius/center of information theory',\n"
    "                'No: B*(s) is an infimum of margins, not a recovery-error radius',\n"
    "                'viewpoint relation; no equivalence claimed'],\n"
    "               ['[9] Kim et al.',\n"
    "                'worst-case SINR, norm-bounded steering error',\n"
    "                'parallel/orthogonal error decomposition bounds',\n"
    "                'No: min{b0 r, a x + b y} not found as a stated envelope there',\n"
    "                'decomposition context; baseline assembled here']],\n"
    "              [.95, 1.6, 1.35, 1.7, 1.0])\n"
    "    else:\n"
    "        raise ValueError(name)")
s = s.replace("    else:\n        raise ValueError(name)", lit)

assert s != orig
open(p, 'w', encoding='utf-8', newline='').write(s)
print('builder patched:', sum(1 for a, b in
      [(1, 1)]) and 'ok')
