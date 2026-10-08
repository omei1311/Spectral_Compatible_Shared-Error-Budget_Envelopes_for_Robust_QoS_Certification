# -*- coding: utf-8 -*-
"""Minimal stdlib DOCX reader/writer for this project (no python-docx).

Read: extract paragraph text (w:t and m:t nodes, in document order).
Write: build a .docx from a list of (style, text) blocks with WordprocessingML.
"""
import re
import zipfile

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
M = "http://schemas.openxmlformats.org/officeDocument/2006/math"


def read_paragraphs(path):
    """Return list of (style, text); math text is kept inline in order."""
    import xml.etree.ElementTree as ET
    with zipfile.ZipFile(path) as z:
        xml = z.read("word/document.xml")
    root = ET.fromstring(xml)
    body = root.find(f"{{{W}}}body")
    out = []
    for p in body.iter(f"{{{W}}}p"):
        style = ""
        pPr = p.find(f"{{{W}}}pPr")
        if pPr is not None:
            st = pPr.find(f"{{{W}}}pStyle")
            if st is not None:
                style = st.get(f"{{{W}}}val", "")
        parts = []
        for node in p.iter():
            tag = node.tag.split("}")[-1]
            if tag == "t" and node.text:
                parts.append(node.text)
        out.append((style, "".join(parts)))
    return out


def read_raw(path):
    with zipfile.ZipFile(path) as z:
        return {n: z.read(n) for n in z.namelist()}


# ---------------- writer ------------------------------------------------
_NS_TYPES = (
    "application/vnd.openxmlformats-officedocument.wordprocessingml"
    ".document.main+xml"
)
_CONTENT_TYPES = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
    '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
    '<Default Extension="rels" ContentType="application/vnd.openxmlformats-'
    'package.relationships+xml"/>'
    f'<Default Extension="xml" ContentType="application/xml"/>'
    f'<Override PartName="/word/document.xml" ContentType="{_NS_TYPES}"/>'
    "</Types>"
)
_RELS = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
    '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/'
    'relationships">'
    '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/'
    'officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>'
    "</Relationships>"
)
_W_RELS = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
    '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/'
    'relationships">'
    '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/'
    'officeDocument/2006/relationships/styles" Target="styles.xml"/>'
    "</Relationships>"
)
_STYLES = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
    f'<w:styles xmlns:w="{W}">'
    f'<w:style w:type="paragraph" w:styleId="Title"><w:name w:val="Title"/>'
    f'<w:pPr><w:jc w:val="center"/><w:spacing w:after="240"/></w:pPr>'
    f'<w:rPr><w:b/><w:sz w:val="36"/></w:rPr></w:style>'
    f'<w:style w:type="paragraph" w:styleId="Heading1"><w:name w:val="heading 1"/>'
    f'<w:pPr><w:outlineLvl w:val="0"/><w:spacing w:before="240" w:after="120"/>'
    f"</w:pPr>"
    f'<w:rPr><w:b/><w:sz w:val="30"/></w:rPr></w:style>'
    f'<w:style w:type="paragraph" w:styleId="Heading2"><w:name w:val="heading 2"/>'
    f'<w:pPr><w:outlineLvl w:val="1"/><w:spacing w:before="200" w:after="100"/>'
    f"</w:pPr>"
    f'<w:rPr><w:b/><w:sz w:val="26"/></w:rPr></w:style>'
    f'<w:style w:type="paragraph" w:styleId="BodyText"><w:name w:val="body"/>'
    f"<w:pPr><w:spacing w:line=\"360\" w:lineRule=\"auto\"/>"
    f'<w:ind w:firstLineChars="0"/></w:pPr>'
    f'<w:rPr><w:sz w:val="21"/></w:rPr></w:style>'
    "</w:styles>"
)


def _esc(s):
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            .replace('"', "&quot;"))


def write_docx(path, blocks):
    """blocks: list of (style, text). style in Title/Heading1/Heading2/Body."""
    ps = []
    for style, text in blocks:
        pPr = f'<w:pPr><w:pStyle w:val="{style}"/></w:pPr>'
        run = (f'<w:r><w:rPr><w:rFonts w:ascii="Cambria Math" '
               f'w:hAnsi="Cambria Math" w:eastAsia="宋体"/></w:rPr>'
               f'<w:t xml:space="preserve">{_esc(text)}</w:t></w:r>')
        ps.append(f"<w:p>{pPr}{run}</w:p>")
    doc = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
        f'<w:document xmlns:w="{W}"><w:body>'
        + "".join(ps)
        + ('<w:sectPr><w:pgSz w:w="11906" w:h="16838"/>'
           '<w:pgMar w:top="1440" w:right="1440" w:bottom="1440" '
           'w:left="1440"/></w:sectPr>')
        + "</w:body></w:document>"
    )
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", _CONTENT_TYPES)
        z.writestr("_rels/.rels", _RELS)
        z.writestr("word/_rels/document.xml.rels", _W_RELS)
        z.writestr("word/styles.xml", _STYLES)
        z.writestr("word/document.xml", doc)


if __name__ == "__main__":
    import sys
    for style, text in read_paragraphs(sys.argv[1]):
        if text.strip() or style:
            print(f"[{style:9s}] {text}")
