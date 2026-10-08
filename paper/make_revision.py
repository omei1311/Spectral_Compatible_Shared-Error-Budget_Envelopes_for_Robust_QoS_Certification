# -*- coding: utf-8 -*-
"""Assemble the revised paper DOCX (v4) from content parts. Keeps the
original v3 draft untouched."""
import docx_io
from paper_content_1 import PART1
from paper_content_2 import PART2

STYLE = {"T": "Title", "H1": "Heading1", "H2": "Heading2", "B": "BodyText"}
OUT = "谱兼容共享误差预算下的鲁棒QoS认证_修订版v4.docx"


def main():
    blocks = [(STYLE[k], t) for k, t in PART1 + PART2]
    docx_io.write_docx(OUT, blocks)
    n = len(blocks)
    print(f"wrote {OUT}: {n} paragraphs "
          f"({sum(1 for s, _ in blocks if s == 'Heading1')} sections)")


if __name__ == "__main__":
    main()
