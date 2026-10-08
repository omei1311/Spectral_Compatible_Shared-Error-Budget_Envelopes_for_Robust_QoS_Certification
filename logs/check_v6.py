# -*- coding: utf-8 -*-
import re

txt = open('logs/v6_check.txt', encoding='utf-8').read()
body = re.sub(r'^\[[^\]]*\]', '', txt, flags=re.M)
must = ['Sharp Robust QoS Certification from Six-Scalar Spectral Summaries',
        'information map', 'certification hierarchy', 'indistinguishable',
        'Physical interpretation', 'Related work', 'radius of information',
        'individually valid worst cases', 'additional channel information',
        '190', '90.05', '24.4', 'near parity', 'no stable runtime advantage',
        '7.5 Complementary', '7.6 Independent exact',
        'Parrott', 'Buzano', 'Grone', 'Micchelli', 'Zheng']
missing = [m for m in must if m not in body]
print('missing:', missing if missing else 'NONE')
print('negation present:', body.count('not a physical channel model'))
print('stray six-scalar-channel-model:', body.count('six-scalar channel model'))
for s in ['2.3 Six-Scalar Information', '2.4 Physical interpretation',
          '2.5 Related work', '3 Spectral compatible', '4 QoS envelopes',
          '5 Boundary behavior', '6.3 S procedure SDP',
          '7.4 Independent validation', '10 Conclusion', 'References']:
    print(('OK  ' if s in body else 'MISS') + ' ' + s)
print('defensive-phrase count:',
      sum(body.count(p) for p in
          ['not claimed', 'do not claim', 'does not establish']))
