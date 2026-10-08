# -*- coding: utf-8 -*-
import re

body = re.sub(r'^\[[^\]]*\]', '', open('logs/v655_check.txt',
                                       encoding='utf-8').read(), flags=re.M)
checks = [
    ('T1 row [1],[4],[13]', '[1],[4],[13]' in body),
    ('no stale RIS [9]', 'RIS/beamforming [1],[4],[9]' not in body),
    ('T12 Kim is [13]', re.search(r'\[13\]\s*Kim et al', body) is not None),
    ('T12 Kim [9] gone', re.search(r'\[9\]\s*Kim', body) is None),
    ('T12 DKW row', re.search(r'\[12\]\s*Davis', body) is not None),
    ('DKW known-result wording',
     'operator completions / norm-preserving extensions under prescribed '
     'blocks or actions' in body),
    ('DKW diff on objective',
     'does not evaluate the specific maximum output functional' in body),
    ('body DKW softened',
     'admissible extensions under prescribed blocks or actions' in body),
    ('no minimal-norm oversimplification',
     'minimal norm of extensions with prescribed action' not in body),
]
fails = [nm for nm, c in checks if not c]
for nm, c in checks:
    print(('PASS ' if c else 'FAIL ') + nm)
print('ALL OK' if not fails else f'FAILS: {fails}')
