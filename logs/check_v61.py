# -*- coding: utf-8 -*-
"""V6.1 DOCX content verification per the revision checklist."""
import re

txt = open('logs/v61_check.txt', encoding='utf-8').read()
body = re.sub(r'^\[[^\]]*\]', '', txt, flags=re.M)

fails = []
def need(name, cond):
    print(('PASS ' if cond else 'FAIL ') + name)
    if not cond:
        fails.append(name)

# 1. inf/min correctness
need('EQ(2d) uses inf (B* = inf over class)', 'inf' in body and 'B^*' in body.replace('B*','B^*')
     or 'B*(s)=inf' in body.replace(' ', '') or 'B* = inf' in body or 'B*(s)= inf' in body)
need('Theorem 2 statement uses inf', 'inf_{(b,h' in body.replace(' ', '') or
     re.search(r'inf\s*\(b,\s*h', body) is not None or 'infimum of the true worst-case' in body)
need('attainment nuance present (not "never attained")',
     'Attainment of the infimum by a single instance is not claimed' in body)
# 2. six scalars vs external params
need('external params declared fixed and outside s',
     'fixed external parameters' in body and 'not part of the summary' in body)
# 3. kappa vs H branches
need('kappa branch defined', 'κ=∥A ĥ∥' in body or 'κ = ∥A ĥ∥' in body)
need('H-branch distinct and not merged',
     'not merged' in body and 'H=∥ĥ∥' in body.replace(' ', ''))
# 4. references complete (no placeholders)
need('no placeholder refs', 'to be finalized' not in body
     and 'details to be verified' not in body)
need('verified refs present (Parrott JFA 30, Filomat 37, LAA 58, TSP 56)',
     '30, no. 3, pp. 311–328' in body.replace('vol. ', '').replace('vol.', '')
     or '311–328' in body)
need('Kim et al. TSP 56(4) 1539–1547 present', '1539–1547' in body)
# 5. separate-gain naming, no baseless "previous method"
need('separate-gain comparison envelope named', 'separate-gain comparison' in body
     or 'separate-gain envelope' in body)
need('baseline attribution honest', 'assembled from the two standard norm constraints' in body
     or 'not attribute it to any specific paper' in body)
# 6. dB / NMSE separation
need('dB = deterministic energy ratio, not NMSE',
     'not a statistical NMSE' in body or 'not a statistical NMSE, which would require' in body)
need('no rho^2=NMSE equation', 'ρ² equals' not in body and 'rho^2 equals' not in body)
# 7. SDP scope limited to validation batch
need('SDP scope limited', 'no cross-validation claim is made beyond this batch' in body)
need('no "every certification statement"', 'every certification statement above' not in body)
need('U wording: agrees within tolerance', 'agrees with the reference within a preset' in body
     or 'agrees within tolerance' in body)
# 8. timing protocol
need('2.12x/2.01x tied to tolerance-stopped reference',
     'median paired speedup of 2.12 over the pure tolerance-stopped' in body)
need('200-iteration table labelled supplementary',
     'supplementary' in body and 'Table 10 reproduces the original 200 iteration' in body)
need('no mixed-round SDP speed claim (4–13x only, one round)',
     '4–13× slower than the tolerance-stopped reference' in body
     and '17×' not in body)
# 9. batch enumeration incl. easy zero
need('easy-pool zero gain enumerated', 'easy pool 0 (saturated)' in body
     or 'zero additional accepts' in body)
need('no "positive in every batch"', 'positive in every batch' not in body)
# 10. L_quad positioning
need('L_quad = 6+kappa branch, not full-matrix at evaluation',
     'not a full-matrix method at evaluation time' in body)
# 11. no stray terminology
need('no "six-scalar channel model"', 'six-scalar channel model' not in body)
# 12. tables 1..12 present
for n in range(1, 13):
    need(f'Table {n} caption present', f'Table {n}.' in body)
# 13. theory intactness markers
for m in ['Theorem 1', 'Proposition 1', 'Theorem 2', 'Theorem 3', 'Lemma 1',
          'Cases', 'cases', 'Appendix A', 'Reproducibility']:
    need(f'marker {m}', m in body)

print('\n' + ('ALL CHECKS PASSED' if not fails else f'{len(fails)} FAILURES: {fails}'))
