# -*- coding: utf-8 -*-
import re

body = re.sub(r'^\[[^\]]*\]', '', open('logs/v653_check.txt',
                                       encoding='utf-8').read(), flags=re.M)
flat = re.sub(r'\s+', '', body)   # newline-tolerant
checks = [
    ('7.7 excludes-construction-only sentence',
     'excludesthisexhibitedlarge-normnear-worst-caseconstruction' in flat),
    ('bracket-consistent wording',
     'consistentwithapotentiallymuchtighterrestricted-classoptimum' in flat),
    ('not-all-constructions caveat',
     'notallbounded-normconstructionsapproachingQG' in flat),
]
fails = [nm for nm, c in checks if not c]
for nm, c in checks:
    print(('PASS ' if c else 'FAIL ') + nm)
print('ALL OK' if not fails else f'FAILS: {fails}')
