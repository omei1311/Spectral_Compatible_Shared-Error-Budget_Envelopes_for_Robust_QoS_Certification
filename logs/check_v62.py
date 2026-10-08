# -*- coding: utf-8 -*-
"""V6.2 consistency sweep on the extracted DOCX text (calibrated to the
actual OMML extraction format: brackets/carets stripped, superscripts
inlined)."""
import re

body = re.sub(r'^\[[^\]]*\]', '',
              open('logs/v62_check.txt', encoding='utf-8').read(), flags=re.M)
fails = []
flat = body.replace(' ', '')


def need(name, cond):
    print(('PASS ' if cond else 'FAIL ') + name)
    if not cond:
        fails.append(name)


# 1. Lemma 1 A-E branches (extraction strips [], ^, and fraction parens)
for tag, frag in [
        ('A zero operator', 'C=0,ĥ=ae1'),
        ('B parallel-only', 'CH=α,0,…,ĥ=ae1'),
        ('C alpha=0', 'CH=0,β,ĥ=a,I0/βT'),
        ('D single-formula hat', 'ĥ=a,I0−αa/δ,0T'),
        ('D outputs', 'CHĥ=I0,0T'),
        ('D three checks', '∥CHP∥2=maxδ,β=βand∥C∥2=maxα2+δ2,β=β0'),
        ('E boundary', 'CH=diagα,β,ĥ=a,J⊥/βT'),
        ('zero-op forces I0', 'α=β=β0=0,henceI0=0by(i)'),
        ('exhaustive claim', 'mutuallyexclusiveandexhaustive'),
        ('padding reason', 'enterneitherthedesirednortheinterferenceterm'),
]:
    need(f'Lemma1 {tag}', frag.replace(' ', '') in flat)

# 2. no unproven minimality
need('no "smallest summary set"', 'smallestsummaryset' not in flat)
need('neutral replacement present',
     'thethreespectralquantitiesofoneoperator' in flat
     and 'fixedsummarysize' in flat)
need('minimal-dimension disclaimer kept (legitimate)',
     'nominimaldimensionclaim' in flat)

# 3. Lemma 2 direction attainment (extraction format)
for tag, frag in [
        ('statement', 'maxC∈𝒟∥CHe∥=Rx,y'),
        ('class def', '∥CHv∥=α,∥CHP∥2=β,∥C∥2=β0'),
        ('char-poly root arg', 'vanishesatλ=β02becauseα2+δ2=β02'),
        ('construction', 'CH=d,ℓ,ℓ2v,u,u'),
        ('beta check', '∥ℓ∥=β02−k2≤β(ask≥ν)'),
        ('alpha=0 case', 'd=0,ℓ=β,0,0T,ℓ2=0,β,0T'),
        ('beta=0 case', 'CH=αvvH'),
        ('P action on z-perp-v', 'oninputsz⊥v'),
        ('separation from QoS class', 'directionproblemonly'),
]:
    need(f'Lemma2 {tag}', frag.replace(' ', '') in flat)

# 4. Theorem 2 attainment enumeration
for tag, frag in [
        ('three proven cases', 'Finiteattainmentisprovedinthreecases'),
        ('r=0', 'r=0(everyinstanceattainsit)'),
        ('alpha=0 beta>0', 'α=0withβ>0(theexplicitconstructionbelow'),
        ('boundary Thm3',
         'boundaryclassofTheorem3(attainedbyanexplicitinstance)'),
        ('not impossibility claim', 'notanimpossibilityclaim'),
        ('left open', 'leftopen'),
]:
    need(f'Thm2 {tag}', frag.replace(' ', '') in flat)
need('Thm2 no blanket unbounded claim',
     'onlyalongnominalchannelsofunboundednorm' not in flat)

# 5. SDP timing same-round
need('no blanket unbounded claim enforced',
     'onlyalongnominalchannelsofunboundednorm' not in flat)
need('17.3 ratio-of-medians', '17.3timesslower' in flat)
need('median pair ratios 15.6/7.2', '15.6' in body and '7.2' in body)
need('15/15 pairs', '15of15pairs' in flat)
need('single-round statement', 'single-roundfigures' in flat)
need('no old 4-13x', '4–13×' not in body and '4-13' not in body)
need('Table-11 scope for 0.032', 'amongthetwopipelinesofTable11' in flat)

# 6. standing invariants
need('inf in Thm2 (12)', 'infb,ĥ,C∈𝒞s' in flat)
need('no six-scalar channel model', 'six-scalarchannelmodel' not in flat)
need('no every-certification claim',
     'everycertificationstatementabove' not in flat)
need('no positive-in-every-batch', 'positiveineverybatch' not in flat)
need('Table count 12', all(f'Table {i}.' in body for i in range(1, 13)))
need('Theorems + Lemmas present',
     all(m in body for m in ['Theorem 1', 'Theorem 2', 'Theorem 3',
                             'Lemma 1', 'Lemma 2']))
need('no literal LaTeX leakage (cdot/frac/text)',
     '\\cdot' not in body and '\\frac' not in body
     and '\\text' not in body)

print('\n' + ('ALL CHECKS PASSED' if not fails
             else f'{len(fails)} FAILURES: {fails}'))
