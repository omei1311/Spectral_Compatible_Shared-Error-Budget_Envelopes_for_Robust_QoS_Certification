# -*- coding: utf-8 -*-
"""V6.5.1 semantic patches on the cleaned paper_v6_5.md (whitespace-
tolerant matching)."""
import re

p = 'paper/paper_v6_5.md'
s = open(p, encoding='utf-8').read()
n = [0]


def sub(pattern, repl, flags=re.S, cnt=1):
    global s
    new, k = re.subn(pattern, repl, s, count=cnt, flags=flags)
    assert k == cnt, (k, pattern[:60])
    s = new
    n[0] += 1


W = r'\s+'


def sp(old_plain, new_plain):
    """Whitespace-tolerant literal replace."""
    pat = W.join(re.escape(w) for w in old_plain.split())
    sub(pat, new_plain.replace('\\', '\\\\'))


# 1. abstract ending
sp("Experiments verify the analytical formulas, exhibit strict but "
   "modest and batch-dependent empirical gains over the separate-gain "
   "envelope, and quantify how norm-restricted information classes "
   "can certify strictly more.",
   "Experiments verify the analytical formulas, exhibit strict but "
   "modest and batch-dependent empirical gains over the separate-gain "
   "envelope, and explore how bounding the nominal-channel norm may "
   "tighten the admissible information class.")

# 2. renumber sections
sp("### 7.7 Timing and implementation sensitivity",
   "### 7.8 Timing and implementation sensitivity")
sp("### 7.6 Effect of a bounded nominal-channel norm",
   "### 7.7 Effect of a bounded nominal-channel norm (exploratory)")

# 3. 5.2 pointer
sp("Section 7.6 quantifies empirically how much a norm restriction "
   "can tighten the class optimum.",
   "Section 7.7 provides exploratory numerical evidence on the "
   "possible effect of norm restrictions.")

# 4. 7.7 conclusion paragraph downgrade (match through 'future work.')
sub(re.escape("Together, the two studies") + r'.*?' + re.escape(
        "H$-restricted envelope is left as future work."),
    "Together, the two studies constitute exploratory numerical "
    "evidence on the possible effect of norm restrictions: norm "
    "restrictions may lead to a tighter class optimum, and the "
    "boundary-approach search indicates the potential gap can be "
    "large for boundary-approaching summaries, although the "
    "restricted optimum is not characterized; establishing any "
    "strict improvement would require a restricted-class lower bound "
    "above $Q_G$, which we do not have. A complete $H$-restricted "
    "envelope is left as future work.")

# 5. Conclusion sentence
sp("The class-optimal value is a statement of information "
   "sufficiency, not instance-wise exactness. The bounded-norm "
   "experiments illustrate that the unrestricted construction can "
   "depend on large-norm geometries, while the actual "
   "restricted-class optimum remains open." if False else
   "The class-optimal value is a statement of information "
   "sufficiency, not instance-wise exactness, and norm- restricted "
   "classes can certify strictly more, as the bounded-norm "
   "experiments quantify.",
   "The class-optimal value is a statement of information "
   "sufficiency, not instance-wise exactness. The bounded-norm "
   "experiments illustrate that the unrestricted construction can "
   "depend on large-norm geometries, while the actual "
   "restricted-class optimum remains open.")

# 6. Discussion limitation (ii)
sub(re.escape("(ii) Information-class sharpness holds for the") + r'.*?'
    + re.escape("large for some summaries."),
    "(ii) Information-class sharpness holds for the unrestricted "
    "class; additional physical restrictions — including a bounded "
    "nominal-channel norm — may lead to a tighter class optimum, and "
    "the boundary-approach experiments indicate the potential gap "
    "can be large for some summaries, although the restricted "
    "optimum is not characterized.")

# 7. Section 2.3 / exact-fallback contradiction fix (two spots)
sp("Finally, the full quadratic data support cheap bounds and exact "
   "trust-region or S-procedure solvers.",
   "Finally, exact certification requires richer, dimension-dependent "
   "information — full quadratic data for witness recovery, or a "
   "fixed-size spectral cache ($λ,w,c$) sufficient for optimal-value "
   "queries — supporting cheap bounds and exact trust-region or "
   "S-procedure solvers.")
sp("exact fallback, in contrast, requires the full quadratic data.",
   "exact fallback, in contrast, requires richer, dimension-dependent "
   "information (full quadratic data for witness recovery, or the "
   "spectral cache ($λ,w,c$) for optimal-value queries).")

# 8. 6.1 attribution -> Lemma 2
sp("the completion in Theorem 2 approaches the scalar extremum for "
   "each direction",
   "the completion in Lemma 2 attains the scalar extremum for each "
   "direction")

# 9. compact-client cross-ref in 7.5
sp("it is not a full-matrix method at evaluation time.",
   "it is not a full-matrix method at evaluation time. The compact "
   "exact client of Section 7.5 stores only ($λ,w,c$) per instance "
   "for optimal-value queries.")

# 10. version bump
sp("Version 6.5 • English revision • 8 October 2026",
   "Version 6.5.1 • English revision • 8 October 2026")

open(p, 'w', encoding='utf-8', newline='').write(s)
print(f'paper_v6_5.md updated ({n[0]} patches)')
