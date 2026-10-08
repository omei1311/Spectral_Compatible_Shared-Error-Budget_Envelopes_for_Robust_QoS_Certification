# -*- coding: utf-8 -*-
"""V6.5.3: novelty-audit additions to paper_v6_5.md (DKW reference,
probabilistic-line citations, related-work sentences)."""
import re

p = 'paper/paper_v6_5.md'
s = open(p, encoding='utf-8').read()
n = [0]


def sp(old, new):
    global s
    pat = r'\s+'.join(re.escape(w) for w in old.split())
    new_s, k = re.subn(pat, new.replace('\\', '\\\\'), s, count=1,
                       flags=re.S)
    assert k == 1, (k, old[:60])
    s = new_s
    n[0] += 1


# 1. Related work: add DKW + probabilistic-line sentences
sp("A second line compresses information probabilistically, using "
   "distributional summaries with conservative approximations.",
   "A second line compresses information probabilistically: chance-"
   "constrained SINR designs convert probabilistic QoS requirements "
   "into conservative deterministic constraints via Bernstein-type "
   "inequalities for quadratic forms [10], with journal-length "
   "treatments and semi-closed forms in [11], and limited-feedback/"
   "covariance-feedback methods quantify what a fixed feedback budget "
   "buys [12]. These lines compress distributional knowledge, whereas "
   "our uncertainty is deterministic and the compressed object is a "
   "fixed set of spectral summaries of one operator.")

# 2. Add DKW to the operator-mathematics paragraph (closest framework)
sp("the chordal completion theorem [6] and "
   "the contraction-completion framework of Parrott [7] place it in a "
   "wider theory without our result being a direct substitution of "
   "either;",
   "the chordal completion theorem [6] and "
   "the contraction-completion framework of Parrott [7] — together "
   "with the norm-preserving dilation theory of Davis–Kahan–"
   "Weinberger [13], which characterizes the minimal norm of "
   "extensions with prescribed action — place it in a wider theory "
   "without our result being a direct substitution of any of them: "
   "these theorems describe the extension set or its minimal norm, "
   "whereas we evaluate the maximum interference output over that set "
   "in closed form;")

# 3. New references [10]-[13]
sp("[9] S.-J. Kim, A. Magnani, A. Mutapcic, S. P. Boyd, and Z.-Q. Luo,",
   "[10] K.-Y. Wang, T.-H. Chang, W.-K. Ma, A. Man-Cho So, and "
   "C.-Y. Chi, \"Probabilistic SINR Constrained Robust Transmit "
   "Beamforming: A Bernstein-Type Inequality Based Conservative "
   "Approach,\" in Proc. IEEE ICASSP, 2011; journal version: K.-Y. "
   "Wang, A. Man-Cho So, T.-H. Chang, W.-K. Ma, and C.-Y. Chi, "
   "\"Outage Constrained Robust Transmit Optimization for Broadcast "
   "Channels: Incorporating Semi-Closed Form Solutions,\" IEEE "
   "Journal of Selected Topics in Signal Processing, vol. 8, no. 1, "
   "pp. 28–42, 2014. (See [11].)\n\n"
   "[11] K.-Y. Wang, A. Man-Cho So, T.-H. Chang, W.-K. Ma, and "
   "C.-Y. Chi, \"Outage Constrained Robust Transmit Optimization for "
   "Broadcast Channels: Incorporating Semi-Closed Form Solutions,\" "
   "IEEE Journal of Selected Topics in Signal Processing, vol. 8, "
   "no. 1, pp. 28–42, 2014. DOI: 10.1109/JSTSP.2013.2290600.\n\n"
   "[12] D. J. Love, R. W. Heath, Jr., V. K. N. Lau, D. Gesbert, "
   "B. D. Rao, and M. Andrews, \"An Overview of Limited Feedback in "
   "Wireless Communication Systems,\" IEEE Journal on Selected Areas "
   "in Communications, vol. 26, no. 8, pp. 1341–1365, 2008. DOI: "
   "10.1109/JSAC.2008.080902.\n\n"
   "[13] C. Davis, W. M. Kahan, and H. F. Weinberger, \"Norm-"
   "preserving dilations and their applications to optimal error "
   "bounds,\" SIAM Journal on Numerical Analysis, vol. 19, no. 3, "
   "pp. 445–469, 1982.\n\n"
   "[14] S.-J. Kim, A. Magnani, A. Mutapcic, S. P. Boyd, and Z.-Q. Luo,")

# 4. Renumber old [9] -> [14] in the reference body and citations
s = s.replace("[9] S.-J. Kim, A. Magnani, A. Mutapcic, S. P. Boyd, "
              "and Z.-Q. Luo,",
              "[14] S.-J. Kim, A. Magnani, A. Mutapcic, S. P. Boyd, "
              "and Z.-Q. Luo,")
s = s.replace("Kim et al.", "Kim et al.")  # names unchanged
for old, new in [("[9],", "[14],"), ("[1],[4],[9]", "[1],[4],[14]"),
                 ("[9]", "[14]")]:
    s = s.replace(old, new)
n[0] += 1

# 5. Fix accidental duplicate ref block from step 3 (first [9] stub)
s = s.replace(
    "[10] K.-Y. Wang, T.-H. Chang, W.-K. Ma, A. Man-Cho So, and "
    "C.-Y. Chi, \"Probabilistic SINR Constrained Robust Transmit "
    "Beamforming: A Bernstein-Type Inequality Based Conservative "
    "Approach,\" in Proc. IEEE ICASSP, 2011; journal version: K.-Y. "
    "Wang, A. Man-Cho So, T.-H. Chang, W.-K. Ma, and C.-Y. Chi, "
    "\"Outage Constrained Robust Transmit Optimization for Broadcast "
    "Channels: Incorporating Semi-Closed Form Solutions,\" IEEE "
    "Journal of Selected Topics in Signal Processing, vol. 8, no. 1, "
    "pp. 28–42, 2014. (See [11].)\n\n"
    "[11] K.-Y. Wang",
    "[10] K.-Y. Wang, T.-H. Chang, W.-K. Ma, A. Man-Cho So, and "
    "C.-Y. Chi, \"Probabilistic SINR Constrained Robust Transmit "
    "Beamforming: A Bernstein-Type Inequality Based Conservative "
    "Approach,\" in Proc. IEEE ICASSP, 2011, pp. 3080–3083. DOI: "
    "10.1109/ICASSP.2011.5946976.\n\n"
    "[11] K.-Y. Wang")

sp("Version 6.5.2 • English revision • 8 October 2026",
   "Version 6.5.3 • English revision • 8 October 2026")

open(p, 'w', encoding='utf-8', newline='').write(s)
print(f'paper_v6_5.md updated ({n[0]} edits)')
