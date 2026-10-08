# -*- coding: utf-8 -*-
import re

body = re.sub(r'^\[[^\]]*\]', '',
              open('logs/v651_check.txt', encoding='utf-8').read(), flags=re.M)
flat = body.replace(' ', '')
checks = [
    ('abstract: explore/may tighten',
     'explorehowboundingthenominal-channelnormmay' in flat),
    ('no certify-strictly-more',
     'can certify strictly more' not in body
     and 'certifystrictlymore' not in flat),
    ('7.7 bounded-H (exploratory)',
     '7.7Effectofaboundednominal-channelnorm(exploratory)' in flat),
    ('7.8 timing',
     '7.8Timingandimplementationsensitivity' in flat),
    ('7.6 cross-val unchanged',
     '7.6Independentexactcross-validation' in flat),
    ('5.2 pointer updated',
     'exploratorynumericalevidenceonthepossibleeffect' in flat),
    ('conclusion: remains open',
     'theunrestrictedconstructioncandependonlarge-normgeometries' in flat
     and 'restricted-classoptimumremainsopen' in flat),
    ('lower-bound honesty',
     'wouldrequirearestricted-classlowerboundaboveQG,whichwedonothave'
     in flat),
    ('discussion (ii) may-lead',
     'mayleadtoatighterclassoptimum' in flat),
    ('2.3 exact-fallback fixed',
     'richer,dimension-dependentinformation' in flat
     and 'spectralcache(λ,w,c)' in flat.replace(' ', '')),
    ('no full-quadratic-must claim',
     'requiresthefullquadraticdata' not in flat),
    ('6.1 Lemma 2 attribution',
     'thecompletioninLemma2attainsthescalarextremum' in flat),
    ('no stray quote artifacts',
     'quantifies empirically "' not in body
     and 'the search found no "' not in body),
]
fails = [nm for nm, c in checks if not c]
for nm, c in checks:
    print(('PASS ' if c else 'FAIL ') + nm)
print('ALL OK' if not fails else f'FAILS: {fails}')
