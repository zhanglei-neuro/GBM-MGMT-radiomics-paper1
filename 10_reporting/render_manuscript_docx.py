# -*- coding: utf-8 -*-
r"""Paper 1 -- Markdown to journal-style Word (.docx).

Typography: A4, 2.54 cm margins, Times New Roman 12 pt, 1.5 line spacing,
three-line tables, centred figure blocks.

Usage
-----
    python render_manuscript_docx.py <input.md> <output.docx> [--figdir DIR]

`--figdir` is only needed when the Markdown embeds figures with the standard
`![caption](file.png)` syntax; relative paths are then resolved against that
folder.

The parser is deliberately small: it understands ATX headings, block quotes,
pipe tables, bullet lists, numbered reference lists and figure blocks. Anything
it does not recognise is emitted as a normal body paragraph, so a new section can
be added to the Markdown without touching this script.
"""
import argparse
import os
import re

from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt, RGBColor

BODY_PT = 12
SMALL_PT = 10.5
TABLE_PT = 9
BODY_FONT = 'Times New Roman'


def _style_run(run, size=BODY_PT, bold=False, italic=False, color=None):
    run.font.name = BODY_FONT
    run.font.size = Pt(size)
    run.bold = bold
    run.italic = italic
    if color:
        run.font.color.rgb = RGBColor(*color)
    # pin the East Asian slot to the same face so mixed text never falls back
    run._element.rPr.rFonts.set(qn('w:eastAsia'), BODY_FONT)


class Renderer(object):
    def __init__(self):
        self.doc = Document()
        sec = self.doc.sections[0]
        sec.page_width, sec.page_height = Cm(21), Cm(29.7)
        sec.top_margin = sec.bottom_margin = Cm(2.54)
        sec.left_margin = sec.right_margin = Cm(2.54)

        normal = self.doc.styles['Normal']
        normal.font.name = BODY_FONT
        normal.font.size = Pt(BODY_PT)

    # ---------------------------------------------------------------- blocks
    def para(self, text='', size=BODY_PT, bold=False, italic=False, align=None,
             space_after=6, line=1.5, first_indent=False, color=None):
        p = self.doc.add_paragraph()
        pf = p.paragraph_format
        pf.space_after = Pt(space_after)
        if line:
            pf.line_spacing_rule = WD_LINE_SPACING.MULTIPLE
            pf.line_spacing = line
        if align is not None:
            p.alignment = align
        if first_indent:
            pf.first_line_indent = Pt(size * 2)
        for chunk in re.split(r'(\*\*.+?\*\*|__.+?__)', text):
            if not chunk:
                continue
            if chunk.startswith('**') and chunk.endswith('**'):
                _style_run(p.add_run(chunk[2:-2]), size, bold=True, italic=italic, color=color)
            elif chunk.startswith('__') and chunk.endswith('__'):
                _style_run(p.add_run(chunk[2:-2]), size, bold=bold, italic=True, color=color)
            else:
                _style_run(p.add_run(chunk), size, bold=bold, italic=italic, color=color)
        return p

    def heading(self, text, level):
        size = {1: 16, 2: 14, 3: 12}.get(level, 12)
        self.para(text.lstrip('#').strip(), size=size, bold=True,
                  align=WD_ALIGN_PARAGRAPH.CENTER if level == 1 else None,
                  space_after=10 if level == 1 else 8)

    @staticmethod
    def _border(cell, top=None, bottom=None):
        tcPr = cell._tc.get_or_add_tcPr()
        borders = OxmlElement('w:tcBorders')
        for tag, size in (('top', top), ('bottom', bottom)):
            if not size:
                continue
            el = OxmlElement('w:' + tag)
            el.set(qn('w:val'), 'single')
            el.set(qn('w:sz'), str(size))
            el.set(qn('w:color'), '000000')
            el.set(qn('w:space'), '0')
            borders.append(el)
        tcPr.append(borders)

    def table(self, rows, caption=None):
        if caption:
            self.para(caption, size=SMALL_PT, bold=True,
                      align=WD_ALIGN_PARAGRAPH.CENTER, space_after=4, line=1.2)
        t = self.doc.add_table(rows=len(rows), cols=len(rows[0]))
        t.alignment = WD_TABLE_ALIGNMENT.CENTER
        for i, row in enumerate(rows):
            for j, value in enumerate(row):
                cell = t.cell(i, j)
                cell.text = ''
                p = cell.paragraphs[0]
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER
                p.paragraph_format.space_after = Pt(2)
                _style_run(p.add_run(str(value).replace('**', '')),
                           TABLE_PT, bold=(i == 0))
                if i == 0:
                    self._border(cell, top=12)
                elif i == 1:
                    self._border(cell, top=6)
                if i == len(rows) - 1:
                    self._border(cell, bottom=12)
        self.doc.add_paragraph().paragraph_format.space_after = Pt(4)

    def figure(self, path, caption='', width_cm=15):
        p = self.doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.add_run().add_picture(path, width=Cm(width_cm))
        if caption:
            self.para(caption, size=SMALL_PT, align=WD_ALIGN_PARAGRAPH.CENTER,
                      space_after=12, line=1.2)

    def page_break(self):
        self.doc.add_page_break()

    def save(self, path):
        self.doc.save(path)


def render(md_path, out_path, figdir=None):
    lines = open(md_path, encoding='utf-8').read().split('\n')
    r = Renderer()
    i, title_done, in_refs = 0, False, False

    while i < len(lines):
        ln = lines[i].rstrip()

        if ln.startswith('# ') and not title_done:
            r.heading(ln, 1)
            title_done = True
            i += 1
            continue
        if ln.startswith('### '):
            r.heading(ln, 3)
            i += 1
            continue
        if ln.startswith('## '):
            in_refs = ln[3:].strip().lower().startswith('reference')
            r.heading(ln, 2)
            i += 1
            continue
        if ln.startswith('> '):
            r.para(ln[2:].strip(), size=SMALL_PT, italic=True,
                   color=(0x59, 0x59, 0x59), space_after=3, line=1.3)
            i += 1
            continue
        if ln.strip() in ('', '---', '***'):
            i += 1
            continue

        # pipe table
        if ln.startswith('|'):
            rows = []
            while i < len(lines) and lines[i].strip().startswith('|'):
                cells = [c.strip() for c in lines[i].strip().strip('|').split('|')]
                if not all(re.fullmatch(r':?-{2,}:?', c) for c in cells):
                    rows.append(cells)
                i += 1
            if rows:
                r.table(rows)
            continue

        # figure block: ![caption](file)
        m = re.match(r'^!\[(.*?)\]\((.+?)\)\s*$', ln)
        if m:
            caption, target = m.group(1), m.group(2)
            if figdir and not os.path.isabs(target):
                target = os.path.join(figdir, target)
            if os.path.exists(target):
                r.figure(target, caption)
            else:
                r.para('[missing figure: %s]' % target, size=SMALL_PT, italic=True)
            i += 1
            continue

        # numbered entry inside the reference list -> hanging indent
        if in_refs and re.match(r'^\d+\.\s', ln):
            p = r.para(ln, size=SMALL_PT, space_after=3, line=1.3)
            p.paragraph_format.left_indent = Cm(0.75)
            p.paragraph_format.first_line_indent = Cm(-0.75)
            i += 1
            continue

        if ln.startswith('- '):
            r.para(u'\u2022 ' + ln[2:], size=BODY_PT - 1, space_after=3, line=1.4)
            i += 1
            continue

        r.para(ln, size=BODY_PT, first_indent=True)
        i += 1

    r.save(out_path)
    return out_path


def main():
    ap = argparse.ArgumentParser(description='Render Paper 1 Markdown to .docx')
    ap.add_argument('md', help='input Markdown file')
    ap.add_argument('out', help='output .docx file')
    ap.add_argument('--figdir', default=None, help='directory holding the figure PNGs')
    args = ap.parse_args()

    print('saved:', render(args.md, args.out, args.figdir))


if __name__ == '__main__':
    main()
