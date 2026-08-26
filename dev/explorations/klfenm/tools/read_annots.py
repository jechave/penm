#!/usr/bin/env python3
"""Print the annotations in a PDF: highlighted text and any attached comment."""
import sys
import fitz

KINDS = {0: "text-note", 8: "highlight", 9: "underline",
         10: "squiggly", 11: "strikeout", 4: "square", 2: "freetext"}

def quad_text(page, annot):
    """The text a markup annotation covers."""
    v = annot.vertices
    if not v:
        return ""
    quads = [fitz.Quad(v[i:i+4]) for i in range(0, len(v), 4)]
    return " ".join(page.get_textbox(q.rect).strip() for q in quads)

def main(path):
    doc = fitz.open(path)
    n = 0
    for pno, page in enumerate(doc, 1):
        for a in page.annots() or []:
            n += 1
            info = a.info
            kind = KINDS.get(a.type[0], a.type[1])
            covered = quad_text(page, a)
            comment = (info.get("content") or "").strip()
            print(f"\n--- p{pno}  [{kind}] ---")
            if covered:
                print(f"  TEXT: {covered}")
            if comment:
                print(f"  NOTE: {comment}")
            if not covered and not comment:
                print(f"  (empty, at {a.rect})")
    print(f"\n{n} annotation(s).")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "klfenm_report.pdf")
