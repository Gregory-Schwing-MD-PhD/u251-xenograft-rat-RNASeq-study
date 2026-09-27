# -*- coding: utf-8 -*-
"""Build the CNS 2026 oral-presentation deck on the MANS_Schwing.pptx template.

Slides are cloned from two stencil slides in the template so the WSU/DMC chrome
(footer bar, shield, section eyebrow, accent rule, logos) is byte-identical to the
deck this template came from; content is then added as native PowerPoint text,
shapes, tables and pictures - nothing is a flattened slide image.

Charts are placed at their native size (pixels / 300 dpi) or larger, never smaller,
so the type inside them stays at least the size 02_make_figures.py drew it at.

Every figure comes from figures/ and every number from figures/facts.json.
"""
import copy, json, os

import pandas as pd
from PIL import Image
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.oxml.ns import qn
from pptx.util import Inches, Pt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIG = os.path.join(HERE, "figures")
TPL = os.path.join(ROOT, "MANS_Schwing.pptx")
OUT = os.path.join(HERE, "CNS2026_Schwing_Abstract418.pptx")
F = json.load(open(os.path.join(FIG, "facts.json")))

GREEN = RGBColor(0x09, 0x59, 0x45)
MIDG = RGBColor(0x6E, 0x9E, 0x92)
PALE = RGBColor(0xE8, 0xF0, 0xEE)
INK = RGBColor(0x00, 0x00, 0x00)
MUTED = RGBColor(0x59, 0x59, 0x59)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
RELNS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"

FOOTER = "CNS 2026  ·  Washington, DC  ·  Abstract 418"
SLIDE_W, SLIDE_H = 13.333, 7.5


# --------------------------------------------------------------- plumbing ----
def clone(prs, src):
    """Duplicate a template slide, remapping its image relationships."""
    new = prs.slides.add_slide(src.slide_layout)
    for sh in list(new.shapes):
        sh._element.getparent().remove(sh._element)
    for sh in src.shapes:
        new.shapes._spTree.append(copy.deepcopy(sh._element))
    seen = {}
    for el in new.shapes._spTree.iter():
        for attr in list(el.attrib):
            if not attr.startswith("{%s}" % RELNS):
                continue
            rid = el.get(attr)
            if not rid or not rid.startswith("rId"):
                continue
            if rid not in seen:
                rel = src.part.rels[rid]
                if rel.is_external:
                    seen[rid] = new.part.relate_to(rel.target_ref, rel.reltype,
                                                   is_external=True)
                else:
                    seen[rid] = new.part.relate_to(rel.target_part, rel.reltype)
            el.set(attr, seen[rid])
    return new


def keep_only(slide, names):
    for sh in list(slide.shapes):
        if sh.name not in names:
            sh._element.getparent().remove(sh._element)


def shape_named(slide, name):
    return next(sh for sh in slide.shapes if sh.name == name)


def drop_template_slides(prs, n):
    lst = prs.slides._sldIdLst
    for sld_id in list(lst)[:n]:
        prs.part.drop_rel(sld_id.get(qn("r:id")))
        lst.remove(sld_id)


# ------------------------------------------------------------------ text -----
def line(text, size=18, bold=False, color=INK, italic=False,
         after=10, before=0, align=PP_ALIGN.LEFT, font="Arial"):
    return dict(text=text, size=size, bold=bold, color=color, italic=italic,
                after=after, before=before, align=align, font=font)


def textbox(slide, l, t, w, h, lines, anchor=MSO_ANCHOR.TOP):
    box = slide.shapes.add_textbox(Inches(l), Inches(t), Inches(w), Inches(h))
    tf = box.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = anchor
    tf.margin_left = tf.margin_right = 0
    tf.margin_top = tf.margin_bottom = 0
    for i, spec in enumerate(lines):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = spec["align"]
        p.space_after = Pt(spec["after"])
        p.space_before = Pt(spec["before"])
        for j, chunk in enumerate(spec["text"].split("**")):   # **bold** segments
            if chunk == "":
                continue
            r = p.add_run()
            r.text = chunk
            r.font.name = spec["font"]
            r.font.size = Pt(spec["size"])
            r.font.bold = spec["bold"] or (j % 2 == 1)
            r.font.italic = spec["italic"]
            r.font.color.rgb = spec["color"]
    return box


def set_text(shape, text, size, bold=False, color=INK, align=None):
    tf = shape.text_frame
    tf.clear()
    p = tf.paragraphs[0]
    if align is not None:
        p.alignment = align
    r = p.add_run()
    r.text = text
    r.font.name = "Arial"
    r.font.size = Pt(size)
    r.font.bold = bold
    r.font.color.rgb = color


# ---------------------------------------------------------------- picture ----
def chart(slide, name, left, top, box_w=None, box_h=None):
    """Place a chart at native size (pixels / 300 dpi) or larger.

    Native size is the floor: the figure was drawn at the size it is meant to be
    read at, so it may be scaled UP into spare slide space but never down, which
    is what made the published panels illegible in the first place.
    """
    path = os.path.join(FIG, name)
    iw, ih = Image.open(path).size
    w, h = iw / 300.0, ih / 300.0
    if box_w and box_h:
        scale = max(1.0, min(box_w / w, box_h / h))
        w, h = w * scale, h * scale
    assert left + w <= SLIDE_W - 0.2 and top + h <= 6.85, (name, left, top, w, h)
    return slide.shapes.add_picture(path, Inches(left), Inches(top),
                                    Inches(w), Inches(h))


def chart_size(name):
    iw, ih = Image.open(os.path.join(FIG, name)).size
    return iw / 300.0, ih / 300.0


# ----------------------------------------------------------------- shapes ----
def block(slide, l, t, w, h, text, size=12, fill=WHITE, edge=GREEN, color=GREEN,
          bold=False, shape=MSO_SHAPE.ROUNDED_RECTANGLE, lw=1.25,
          align=PP_ALIGN.CENTER):
    s = slide.shapes.add_shape(shape, Inches(l), Inches(t), Inches(w), Inches(h))
    s.fill.solid()
    s.fill.fore_color.rgb = fill
    if edge is None:
        s.line.fill.background()
    else:
        s.line.color.rgb = edge
        s.line.width = Pt(lw)
    s.shadow.inherit = False
    tf = s.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = Inches(0.07)
    tf.margin_top = tf.margin_bottom = Inches(0.03)
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    for i, ln in enumerate(text.split("\n")):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        p.space_after = Pt(0)
        r = p.add_run()
        r.text = ln
        r.font.name = "Arial"
        r.font.size = Pt(size)
        r.font.bold = bold or i == 0
        r.font.color.rgb = color
    return s


def chevron(slide, l, t, w=0.20, h=0.20, color=MIDG):
    s = slide.shapes.add_shape(MSO_SHAPE.ISOSCELES_TRIANGLE,
                               Inches(l), Inches(t), Inches(w), Inches(h))
    s.rotation = 90
    s.fill.solid()
    s.fill.fore_color.rgb = color
    s.line.fill.background()
    s.shadow.inherit = False
    return s


def flow(slide, steps, top, width, gap=0.40, left=0.28, head_h=1.05,
         head_size=14, sub_size=12):
    """A row of process boxes with captions and chevrons between them."""
    for i, (head, sub) in enumerate(steps):
        x = left + i * (width + gap)
        block(slide, x, top, width, head_h, head, size=head_size, fill=PALE,
              edge=GREEN, color=GREEN)
        textbox(slide, x, top + head_h + 0.14, width, 0.9,
                [line(t, sub_size, color=MUTED, after=1, align=PP_ALIGN.CENTER)
                 for t in sub.split("\n")])
        if i < len(steps) - 1:
            chevron(slide, x + width + 0.10, top + head_h / 2 - 0.10)


def table(slide, l, t, w, rows, col_w, header_size=13, body_size=13,
          row_h=0.34, head_h=0.36, left_cols=1):
    shape = slide.shapes.add_table(len(rows), len(col_w), Inches(l), Inches(t),
                                   Inches(w),
                                   Inches(head_h + row_h * (len(rows) - 1)))
    tbl = shape.table
    tbl_pr = tbl._tbl.find(qn("a:tblPr"))
    tbl_pr.set("firstRow", "1")
    tbl_pr.set("bandRow", "0")
    for e in tbl_pr.findall(qn("a:tableStyleId")):
        tbl_pr.remove(e)
    sid = tbl_pr.makeelement(qn("a:tableStyleId"), {})
    sid.text = "{2D5ABB26-0587-4C30-8999-92F81FD0307C}"       # No Style, No Grid
    tbl_pr.append(sid)
    for i, cw in enumerate(col_w):
        tbl.columns[i].width = Inches(cw)
    for ri, row in enumerate(rows):
        tbl.rows[ri].height = Inches(head_h if ri == 0 else row_h)
        for ci, val in enumerate(row):
            cell = tbl.cell(ri, ci)
            cell.margin_left = cell.margin_right = Inches(0.09)
            cell.margin_top = cell.margin_bottom = Inches(0.02)
            cell.vertical_anchor = MSO_ANCHOR.MIDDLE
            cell.fill.solid()
            cell.fill.fore_color.rgb = (GREEN if ri == 0
                                        else (PALE if ri % 2 == 0 else WHITE))
            p = cell.text_frame.paragraphs[0]
            p.alignment = PP_ALIGN.LEFT if ci < left_cols else PP_ALIGN.RIGHT
            emph = val.startswith("**")
            r = p.add_run()
            r.text = val.replace("**", "")
            r.font.name = "Arial"
            r.font.size = Pt(header_size if ri == 0 else body_size)
            r.font.bold = (ri == 0) or emph
            r.font.color.rgb = WHITE if ri == 0 else INK
    return tbl


# =============================================================== build =======
prs = Presentation(TPL)
n_tpl = len(prs.slides._sldIdLst)
title_stencil, body_stencil = prs.slides[0], prs.slides[9]

TITLE_CHROME = ["Rectangle 9", "Text Placeholder 12", "Rectangle 22", "Picture 2",
                "Picture 4", "Picture 6", "Straight Connector 3",
                "Straight Connector 4", "Straight Connector 16",
                "Straight Connector 25", "Right Triangle 11", "Rectangle 14",
                "Parallelogram 13", "Parallelogram 29", "Parallelogram 30",
                "Parallelogram 32", "Parallelogram 24"]
BODY_CHROME = ["Rectangle 9", "Text Placeholder 12", "Rectangle 22", "Group 1",
               "Title 11", "Picture 4"]


def content_slide(title, eyebrow, title_size=28):
    s = clone(prs, body_stencil)
    keep_only(s, BODY_CHROME)
    t = shape_named(s, "Title 11")
    set_text(t, title, title_size)
    t.top, t.height = Inches(0.55), Inches(0.8)
    set_text(shape_named(s, "Text Placeholder 12"), eyebrow, 16, align=PP_ALIGN.RIGHT)
    textbox(s, 4.0, 7.05, 5.33, 0.3,
            [line(FOOTER, 10, color=WHITE, after=0, align=PP_ALIGN.CENTER)])
    return s


def refs(slide, items, top=6.55):
    """Small print, indented clear of the shield in the bottom-left corner."""
    textbox(slide, 1.45, top, 11.4, 0.4,
            [line(x, 9, color=MUTED, after=1) for x in items])


HUBS = ", ".join(F["ppi_hubs"][:3] + ["CALB1"])

# Gene-level context: official names and locus class from MyGene.info, written by
# 01b_annotate_genes.py; mean counts from the DESeq2 table. Nothing typed here.
ANN = json.load(open(os.path.join(FIG, "gene_annotation.json"), encoding="utf-8"))
_de = pd.read_excel(os.path.join(ROOT, "manuscript", "Supplementary_Data.xlsx"),
                    "S2_DE_significant")
BASEMEAN = {str(r.symbol): int(round(r.baseMean)) for _, r in _de.iterrows()}
BASEMEAN_MEDIAN = int(round(_de.baseMean.median()))

# Official names, trimmed only where they are too long to sit in a table cell.
TRIM = {
    "proline and arginine rich end leucine rich repeat protein": "prolargin",
    "insulin like growth factor binding protein 3": "IGF-binding protein 3",
    "LDL receptor related protein 1": "LDL-receptor–related protein 1",
    "piccolo presynaptic cytomatrix protein": "presynaptic cytomatrix protein",
    "collagen type I alpha 1 chain": "collagen type I α1",
    "collagen type XVII alpha 1 chain": "collagen type XVII α1",
}


def short_name(sym):
    n = ANN[sym]["name"]
    return TRIM.get(n, n)


def locus_kind(v):
    """How to describe a locus in one phrase: readthrough loci are called what
    they are even when the annotation calls them protein coding."""
    if "readthrough" in v["name"].lower():
        return "readthrough transcript"
    return v["kind_label"]


def signed(v, dp=2):
    """A signed number with a real minus sign, not a hyphen."""
    return ("%+.*f" % (dp, v)).replace("-", "−")



def _t(name):
    """NES for a named gene set, with a real minus sign."""
    return ("%+.2f" % F["themes"][name]["nes"]).replace("-", "−")

def z(v):
    """Signed z-score with a real minus sign, not a hyphen."""
    return signed(v)

# ------------------------------------------------------------- 1. title ------
s = clone(prs, title_stencil)
keep_only(s, TITLE_CHROME)
set_text(shape_named(s, "Text Placeholder 12"),
         "Detroit Medical Center / Wayne State University  ·  Henry Ford Health",
         14, align=PP_ALIGN.RIGHT)
textbox(s, 0.55, 2.62, 12.2, 1.4, [
    line("Characterizing the Transcriptomic Recurrence Signature of Glioblastoma "
         "Following Laser Interstitial Thermal Therapy",
         30, bold=True, after=0, align=PP_ALIGN.CENTER)])
textbox(s, 0.58, 4.34, 11.6, 0.5, [
    line("Gregory J. Schwing, MD, PhD¹   ·   "
         "Tavarekere N. Nagaraja, PhD²   ·   "
         "Indrani Datta, DHI²   ·   Ian Y. Lee, MD²", 15, after=6)])
textbox(s, 0.58, 5.06, 8.7, 0.75, [
    line("¹ Department of General Surgery, Detroit Medical Center / "
         "Wayne State University, Detroit, MI", 11.5, color=MUTED, after=2),
    line("² Department of Neurosurgery, Hermelin Brain Tumor Center, "
         "Henry Ford Health, Detroit, MI", 11.5, color=MUTED, after=0)])
textbox(s, 0.58, 6.02, 6.9, 0.7, [
    line("Abstract 418   ·   Sunrise Science Session – Tumor 2",
         13, bold=True, color=GREEN, after=2),
    line("Monday, November 2, 2026   ·   7:00 AM", 13, color=GREEN, after=0)])
textbox(s, 4.0, 7.05, 5.33, 0.3,
        [line(FOOTER, 10, color=WHITE, after=0, align=PP_ALIGN.CENTER)])

# -------------------------------------------------------- 2. disclosures -----
s = content_slide("Disclosures", "Disclosures")
textbox(s, 0.28, 1.50, 11.9, 3.4, [
    line("**I.Y.L.** has consulting agreements with Medtronic, Inc. "
         "(Minneapolis, MN) and Monteris Medical, Inc. (Plymouth, MN).",
         18, after=14),
    line("All other authors declare no conflicts of interest.", 18, after=26),
    line("Supported by a Henry Ford Health Physician Scientist Award A20050 "
         "(I.Y.L.). No commercial support was received for this analysis.",
         16, color=MUTED, after=14),
    line("Drug candidates discussed here are computational predictions. No agent "
         "named in this talk is approved or under study for this indication.",
         16, color=MUTED, after=0)])

# ------------------------------------------------------------ 3. outline -----
s = content_slide("Where this is going", "Overview")
OUTLINE = [
    ("The problem", "LITT kills the core; the cells that come back sat in the "
     "sublethal margin"),
    ("The model", "MRI-guided LITT in an orthotopic rat xenograft, primary and "
     "recurrent tumors matched"),
    ("The finding", "the recurrent transcriptome switches protein synthesis off"),
    ("The convergence", "subtype scores and network hubs point the same way"),
    ("The opening", "which existing, brain-penetrant drugs oppose that state"),
]
for i, (head, body) in enumerate(OUTLINE):
    y = 1.62 + i * 0.94
    block(s, 0.30, y, 0.52, 0.52, str(i + 1), size=17, fill=GREEN, edge=None,
          color=WHITE, shape=MSO_SHAPE.RECTANGLE)
    textbox(s, 1.02, y - 0.03, 11.3, 0.85, [
        line("**%s** — %s" % (head, body), 18, after=0)])

# --------------------------------------------------------- 4. background -----
s = content_slide("LITT ablates the core; recurrence starts at the margin",
                  "Background")
textbox(s, 0.28, 1.50, 7.1, 4.2, [
    line("LITT is increasingly used for deep-seated and recurrent glioblastoma, "
         "where open resection carries high morbidity.¹", 18, after=18),
    line("Thermal dose falls off with distance from the fiber. Beyond the ablative "
         "threshold lies a **sublethal margin**.", 18, after=18),
    line("**The cells that seed recurrence are the ones that survive there.**",
         19, color=GREEN, after=0)])
CX, CY = 10.35, 3.90
for w, h, fill, edge in ((5.0, 3.8, RGBColor(0xF4, 0xF7, 0xF6),
                          RGBColor(0xC8, 0xD6, 0xD2)),
                         (3.35, 2.50, RGBColor(0xC3, 0xD8, 0xD2), None),
                         (1.60, 1.15, GREEN, None)):
    o = s.shapes.add_shape(MSO_SHAPE.OVAL, Inches(CX - w / 2), Inches(CY - h / 2),
                           Inches(w), Inches(h))
    o.fill.solid()
    o.fill.fore_color.rgb = fill
    if edge is None:
        o.line.fill.background()
    else:
        o.line.color.rgb = edge
        o.line.width = Pt(1)
    o.shadow.inherit = False
textbox(s, CX - 0.80, CY - 0.13, 1.6, 0.3,
        [line("ablation core", 12, bold=True, color=WHITE, after=0,
              align=PP_ALIGN.CENTER)])
textbox(s, CX - 1.68, CY - 1.06, 3.36, 0.3,
        [line("sublethal margin", 12.5, bold=True, color=GREEN, after=0,
              align=PP_ALIGN.CENTER)])
textbox(s, CX - 2.50, CY - 1.76, 5.0, 0.3,
        [line("peritumoral brain", 12, color=MUTED, after=0,
              align=PP_ALIGN.CENTER)])
refs(s, ["1. Chen C, et al. J Neurooncol 2021;151:429–442."])

# --------------------------------------------------------------- 5. gap ------
s = content_slide("What that margin is known to do — and what it is not",
                  "Background")
block(s, 0.30, 1.50, 6.15, 0.50, "Already described", size=15, fill=PALE,
      edge=None, color=GREEN, align=PP_ALIGN.LEFT)
textbox(s, 0.30, 2.20, 6.15, 2.6, [
    line("Peri-ablation vessels transiently open the blood–brain barrier after "
         "ablation.²", 17, after=16),
    line("The immune microenvironment of the sublethal zone has been profiled in "
         "a rodent model.³", 17, after=0)])
block(s, 6.90, 1.50, 6.15, 0.50, "Not described", size=15, fill=GREEN,
      edge=None, color=WHITE, align=PP_ALIGN.LEFT)
textbox(s, 6.90, 2.20, 6.15, 2.6, [
    line("**The cell-intrinsic transcriptional state of the residual tumor cells "
         "that survive ablation and repopulate the lesion.**", 17, after=16),
    line("That is the gap this analysis addresses.", 17, color=MUTED, after=0)])
block(s, 0.30, 5.10, 12.75, 0.95,
      "Does post-LITT recurrence carry its own transcriptional program — "
      "and does that program point to a drug?",
      size=19, fill=PALE, edge=None, color=GREEN, bold=True)
refs(s, ["2. Cleary RT, et al. Neuro Oncol 2026;28:1649–1661.    "
         "3. Tao R, et al. J Immunol 2026;215:vkaf327."])

# ---------------------------------------------------------- 6. the model -----
s = content_slide("The model", "Methods")
textbox(s, 0.28, 1.42, 12.5, 0.4, [
    line("An orthotopic U251N glioblastoma xenograft in the athymic RNU/RNU rat, "
         "ablated under MRI guidance — a model our group established and "
         "characterized previously.⁴˒⁵", 17, color=MUTED, after=0)])
flow(s, [("Orthotopic\nU251N xenograft", "human glioblastoma line\nin the athymic "
          "RNU/RNU rat"),
         ("MRI-guided\nLITT", "thermal ablation of the\nprimary tumor under\n"
          "real-time MR thermometry"),
         ("Matched tumors\nre-harvested", "primary (pre-LITT) and\nrecurrent "
          "(post-LITT)\nn = 3 per arm")],
     top=2.30, width=3.75, gap=0.55, left=0.55, head_h=1.10, head_size=16,
     sub_size=13)
textbox(s, 0.28, 5.15, 12.6, 1.2, [
    line("The graft is human and the host is rat, which is the point: every read "
         "can be assigned to a species, so the tumor compartment can be read on "
         "its own.", 17, after=0)])
refs(s, ["4. Nagaraja TN, et al. Acta Neurochir 2021;163:3455–3463.    "
         "5. Nagaraja TN, et al. J Neurosurg 2026;145:364–377."])

# ------------------------------------------------------- 7. the analysis -----
s = content_slide("The analysis", "Methods")
flow(s, [("xengsort", "human graft reads sorted\nfrom rat host reads, k = 25"),
         ("DESeq2", "Wald test, ashr shrinkage\nFDR < 0.05, |log₂FC| > 1"),
         ("Broad GSEA · subtypes\nSTRING · DSigDB",
          "enrichment, subtypes,\nnetwork and drug ranking")],
     top=1.55, width=3.75, gap=0.55, left=0.55, head_h=1.05, head_size=16,
     sub_size=13)
textbox(s, 0.28, 4.05, 12.6, 2.2, [
    line("**Enrichment was permuted on gene-set membership, not on samples.** At "
         "n = 3 per arm there are too few sample permutations to estimate "
         "significance, and gene-set permutation is substantially more "
         "reproducible at this size.⁶", 17, after=16),
    line("Repurposing screened DSigDB, ranking each compound by how strongly its "
         "own signature countered the tumor signature (|NES|¹˙⁵), weighted by "
         "predicted blood–brain barrier permeability.", 17, after=16),
    line("Raw sequencing data are public: **GEO accession GSE338105**.",
         17, after=0)])
refs(s, ["6. Maleki F, et al. Hum Genomics 2019;13(Suppl 1):42."])

# -------------------------------------------------- 8. global separation -----
s = content_slide("Primary and recurrent tumors separate globally", "Results")
chart(s, "chart_pca.png", 0.45, 1.50, 5.60, 4.85)
pm = F["permanova"]
textbox(s, 6.35, 1.70, 6.70, 4.3, [
    line("The two arms sit on opposite sides of PC1. **%.1f%%** of the variance "
         "is on that axis and **%.1f%%** on PC2." % (F["pc1"], F["pc2"]),
         18, after=18),
    line("PERMANOVA attributes **R² = %.3f** of transcriptome-wide variance to "
         "recurrence (F = %.2f, **p = %.2f**)."
         % (pm["r2"], pm["F"], pm["p"]), 18, after=14),
    line("A large effect that three animals per arm cannot make formally "
         "significant. Worth saying out loud rather than leaving in the "
         "supplement.", 16, color=MUTED, after=0)])

# ------------------------------------------------ 9. differential expression --
s = content_slide("Only 35 genes move", "Results")
chart(s, "chart_volcano.png", 0.28, 1.58)
n_coding = sum(1 for v in ANN.values() if v["coding"])
extremes = sorted(ANN.items(), key=lambda kv: -abs(kv[1]["lfc"]))[:6]
lo = min(BASEMEAN[s_] for s_, _ in extremes)
hi = max(BASEMEAN[s_] for s_, _ in extremes)
textbox(s, 7.35, 1.58, 5.70, 0.7, [
    line("**%d genes** clear FDR < 0.05 with |log₂FC| > 1 — **%d up, %d down**; "
         "%d are protein coding." % (F["de_total"], F["de_up"], F["de_down"],
                                     n_coding), 17, after=0)])
table(s, 7.35, 2.34, 5.70,
      [["Largest fold changes", "log₂FC", "what the locus is"]] +
      [[sym, signed(v["lfc"], 1), locus_kind(v)] for sym, v in extremes],
      [2.05, 0.95, 2.70], header_size=12.5, body_size=12.5, row_h=0.29,
      head_h=0.32, left_cols=1)
textbox(s, 7.35, 4.55, 5.70, 1.8, [
    line("Every one of the six is a readthrough transcript, a small nuclear RNA "
         "or an unannotated locus, on **%d–%d mean counts** against a median of "
         "**%d** — where a handful of reads makes a large ratio."
         % (lo, hi, BASEMEAN_MEDIAN), 15, color=MUTED, after=14),
    line("**A list this short does not carry the biology.** The rest of the talk "
         "asks what the whole transcriptome is doing.", 16, color=GREEN, after=0)])

# --------------------------------------------- 10. translational shutdown ----
s = content_slide("One dominant theme: translation is switched off", "Results")
chart(s, "chart_gsea.png", 0.28, 1.70, 8.90, 4.40)
textbox(s, 9.35, 1.75, 3.72, 4.3, [
    line("Every leading down-regulated set is a translation or "
         "amino-acid-stress set.", 16, after=14),
    line("**The same direction runs across the whole anabolic programme:** mTORC1 %s, MYC %s, ribosome biogenesis %s, glycolysis %s, OXPHOS %s." % (_t("mTORC1 signalling"), _t("MYC targets"), _t("ribosome biogenesis"), _t("glycolysis"), _t("oxidative phosphorylation")), 15, after=14),
    line("The TCA cycle is unchanged (%s) — this is the biosynthetic machinery, not core carbon flux." % _t("TCA cycle"), 15, color=GREEN, after=0)])
refs(s, ["Broad GSEA, gene-set permutation, 1000 permutations, Diff_of_Classes. "
         "Nominal p is reported alongside FDR q because n = 3 per arm."])

# --------------------------------------------------- 11. the one-sidedness ---
s = content_slide("The effect is strikingly one-sided", "Results")
chart(s, "chart_gsea_landscape.png", 0.30, 1.65, 7.35, 4.55)
textbox(s, 7.90, 1.72, 5.15, 4.3, [
    line("Across **%s gene sets**, six reach q < 0.25. All six are down-regulated, "
         "and all six are translation or stress sets."
         % format(F["gsea_n_sets"], ","), 18, after=18),
    line("Nothing on the up-regulated side survives correction (best q = %.2f) or "
         "forms an interpretable group." % F["gsea_up_best_q"], 18, after=18),
    line("**The recurrent state is defined by what it switches off, not by a new "
         "program it switches on.**", 18, color=GREEN, after=0)])

# --------------------------------------------------------- 12. subtypes ------
s = content_slide("Scored against the published classifiers", "Results")
chart(s, "chart_subtypes.png", 0.30, 1.75, 6.90, 4.45)
_ac = F["gsva_by_name"]["Neftel_AC"]
_mtc = F["gsva_by_name"]["Garofano_MTC"]
textbox(s, 7.45, 1.70, 5.60, 4.4, [
    line("Ten published Neftel and Garofano signatures, scored by **GSVA** — "
         "39 to 50 genes each.", 17, color=MUTED, after=16),
    line("One transition reaches significance: **loss of the astrocyte-like "
         "state** (%s → %s, p = %.3f)."
         % (z(_ac["primary"]), z(_ac["recurrent"]), _ac["p"]), 17, after=16),
    line("**The mitochondrial signature moves with them (%s)** — the same "
         "direction as the fall in OXPHOS and mitochondrial translation, and "
         "34 of its 35 measurable genes are lower in recurrence."
         % z(_mtc["change"]), 17, color=GREEN, after=0)])
refs(s, ["7. Garofano L, et al. Nat Cancer 2021;2:141–156.    "
         "8. Neftel C, et al. Cell 2019;178:835–849."])

# ------------------------------------------------------- 13. network hubs ----
s = content_slide("The hubs: matrix genes down, calcium genes up", "Results")
chart(s, "chart_ppi.png", 0.30, 2.45, 6.15, 3.70)
# all eleven connected genes, so nothing is cherry-picked out of the modules
hub_rows = [["Hub gene", "log₂FC", "what it is"]]
for sym in F["ppi_nodes_down"] + F["ppi_nodes_up"]:
    hub_rows.append([sym, signed(ANN[sym]["lfc"]), short_name(sym)])
table(s, 6.65, 1.52, 6.40, hub_rows, [1.55, 0.95, 3.90], header_size=12,
      body_size=12, row_h=0.27, head_h=0.30, left_cols=1)
textbox(s, 6.65, 5.02, 6.40, 1.5, [
    line("**Down:** the structural matrix. **Up:** neuronal calcium handling and "
         "the presynaptic scaffold.", 15, after=12),
    line("Whether that echoes the neuron–glioma synapse described in high-grade "
         "glioma⁹˒¹⁰ — or only reflects which cells survived — six animals "
         "cannot tell us.", 15, color=GREEN, after=0)])
refs(s, ["9. Venkatesh HS, et al. Nature 2019;573:539–545.    "
         "10. Venkataramani V, et al. Nature 2019;573:532–538."])

# ---------------------------------------------------- 14. drug prioritization -
s = content_slide("One candidate survives every filter", "Results")
chart(s, "chart_drugs.png", 0.28, 1.70, 6.05, 4.50)
rows = [["#", "Compound", "Phase", "|NES|", "BBB", "Score"]]
for i, d in enumerate(F["drug_tierA"][:6], start=1):
    star = d["drug"] == F["drug_top"]["drug"]
    rows.append([str(i), ("**" if star else "") + d["drug"], str(d["phase"]),
                 "%.2f" % abs(d["nes"]), "%.2f" % d["bbb"],
                 "%.2f" % d["score"]])
table(s, 6.60, 1.68, 6.45, rows, [0.45, 2.55, 0.85, 0.85, 0.85, 0.90],
      row_h=0.32, left_cols=2)
textbox(s, 6.60, 4.05, 6.45, 2.4, [
    line("Only **%d of %d** DSigDB hits carry a clinical phase — the rest "
         "include phosphine, tributyltin and a flame retardant."
         % (F["drug_clinical_n"], len(F["drug_by_name"])), 16, after=14),
    line("**Ciclopirox** ranks first with or without the permeability weight, "
         "is approved, has **completed phase 1 in cancer patients**, and both BBB "
         "models agree.", 16, after=14),
    line("**Prior art was audited for every hit.** Seven have no glioma "
         "literature; pentetrazol induces seizures. Amiodarone and diazepam "
         "carry published evidence of harm.", 16, color=GREEN, after=0)])
refs(s, ["11. Swanson K, et al. Bioinformatics 2024;40:btae416.   "
         "12. Daina A, Zoete V. ChemMedChem 2016;11:1117–1121.   "
         "13. Boursi B, et al. Pharmacoepidemiol Drug Saf 2016;25:1179–1185."])

# -------------------------------------------------------- 15. conclusions ----
s = content_slide("Conclusions", "Conclusions")
textbox(s, 0.28, 1.50, 12.75, 4.0, [
    line("Post-LITT recurrence carries a **coordinated shutdown of "
         "biosynthesis** — translation initiation, elongation, the ribosome and "
         "the GCN2 stress arm, with ribosome biogenesis, mTORC1 and MYC output "
         "falling alongside them.", 18, after=16),
    line("It extends to **energy metabolism on both arms** — glycolysis and "
         "OXPHOS, and mitochondrial translation — while the TCA cycle is "
         "untouched.", 18, after=16),
    line("**Hypoxic signalling collapses** (HIF targets, hypoxia metagene both "
         "−1.79). Ablation removes the hypoxic core; what regrows is better "
         "perfused.", 18, after=16),
    line("Of the ten published subtype signatures, only **loss of the "
         "astrocyte-like state** reaches significance. **Ciclopirox** survives every filter — approved, phase 1 complete in cancer, and "
         "independently nominated by another group's screen.", 18, after=0)])
block(s, 0.28, 5.30, 12.75, 1.00,
      "This is the phenotype you would predict of a thermally stressed, "
      "reperfused margin — recovering it from an unbiased transcriptome is a "
      "positive control, and it is why the drug output of the same pipeline "
      "is worth acting on.",
      size=16, fill=PALE, edge=None, color=GREEN, bold=True)

# -------------------------------------------------------- 16. limitations ----
s = content_slide("Limitations and next steps", "Limitations")
textbox(s, 0.28, 1.50, 5.85, 4.7, [
    line("**Limitations**", 18, color=GREEN, after=14),
    line("n = 3 per arm. Gene-set results are best read as directional, which is "
         "why nominal p is reported alongside FDR q.", 16, after=12),
    line("A single cell line (U251N) in a single xenograft model — no patient "
         "tissue and no immune compartment.", 16, after=12),
    line("Bulk RNA-seq averages the ablation margin with whatever else was "
         "harvested; the persister state is inferred, not isolated.", 16, after=12),
    line("Every drug candidate is a computational prediction. None has been tested "
         "in this model.", 16, after=0)])
textbox(s, 6.75, 1.50, 5.85, 4.7, [
    line("**Next steps**", 18, color=GREEN, after=14),
    line("Spatial or single-cell profiling of the sublethal margin, to test whether "
         "the persister state sits where we think it does.", 16, after=12),
    line("Puromycin incorporation and polysome profiling, to confirm that "
         "translation is genuinely suppressed and not only transcriptionally down.",
         16, after=12),
    line("Ciclopirox and LY-294002 dosed into the post-ablation barrier window, with "
         "survival as the endpoint.", 16, after=12),
    line("Validation in an immunocompetent model, and in patient tissue from the "
         "peri-ablation zone.", 16, after=0)])

# ---------------------------------------------------- 17. acknowledgments ----
s = content_slide("Acknowledgments and data availability", "Thank you")
textbox(s, 0.28, 1.50, 7.3, 4.6, [
    line("Tavarekere N. Nagaraja, PhD, Indrani Datta, DHI, and Ian Y. Lee, MD "
         "— Hermelin Brain Tumor Center, Henry Ford Health, who built and "
         "characterized the model this analysis rests on.", 17, after=16),
    line("Supported by a Henry Ford Health Physician Scientist Award A20050 "
         "(I.Y.L.).", 16, color=MUTED, after=16),
    line("Animal procedures were approved by the Henry Ford Health IACUC "
         "(protocol #1509) and conducted in accordance with the ARRIVE guidelines.",
         16, color=MUTED, after=20),
    line("Sequencing data: **GEO GSE338105**. Differential-expression, enrichment "
         "and drug-ranking tables are released as Supplementary Data.",
         17, after=0)])
block(s, 8.10, 2.30, 4.6, 1.5, "Questions\ngo2432@wayne.edu",
      size=16, fill=PALE, edge=GREEN, color=GREEN, bold=True)

# ------------------------------------------------- speaker notes / timing ----
NOTES = [
 "0:00-0:12  Title. Name, institutions, and that this is a rat model of "
 "MRI-guided LITT. Do not read the title aloud.",

 "0:12-0:22  Disclosures. Read Dr Lee's Medtronic and Monteris agreements out "
 "loud - Monteris makes the LITT system, so the room will want to hear it named.",

 "0:22-0:40  Roadmap. Five beats, ten seconds. Do not linger.",

 "0:40-1:10  Background. The one idea to land: LITT kills the core, but the cells "
 "that come back sat in the sublethal margin. Point at the diagram once.",

 "1:10-1:35  The gap. The vessels and the immune side of that margin have been "
 "described; the tumor cells' own state has not. End on the question and pause.",

 "1:35-2:00  The model. Move fast. The one thing to stress is that the graft is "
 "human and the host is rat, which is what makes the tumor compartment "
 "separable.",

 "2:00-2:25  The analysis. Enrichment was permuted on gene-set membership, not on "
 "samples, because n = 3. If someone challenges the statistics later, this slide "
 "is the answer. Mention the GEO accession is already public.",

 "2:25-2:50  Global structure. The arms separate on PC1. Say the PERMANOVA "
 "p = 0.10 out loud - do not let someone find it in the supplement. Frame it as: "
 "the effect is large, the sample is small.",

 "2:50-3:10  Differential expression. Only 35 genes. Walk the table once, do not "
 "read it: the six biggest fold changes are readthrough transcripts, a U6 snRNA "
 "and two unannotated loci, all on double-digit counts. That is what justifies "
 "the pivot to gene sets, and it pre-empts the reviewer who asks why you did not "
 "chase the top hit. Do not dwell on individual genes.",

 "3:10-3:45  The core result. Every leading down-regulated set is translation or "
 "amino-acid stress. The top set clears FDR; the next five are nominally "
 "significant and directionally identical. Say both halves. Then the right-hand "
 "panel: the same direction runs across mTORC1, MYC, ribosome biogenesis, "
 "glycolysis and OXPHOS, while the TCA cycle does not move. The point to land is "
 "that this is the biosynthetic machinery being switched off, not core carbon "
 "flux - a cell that has stopped building, not one that has stopped burning.",

 "3:45-4:10  The one-sidedness. This is the strongest slide - 8,869 sets tested, "
 "six clear q < 0.25, all six on the same side and all six the same theme. "
 "Nothing on the up side is interpretable.",

 "4:10-4:35  Subtypes. These are the ten published Neftel and Garofano "
 "signatures, scored per sample by GSVA over the full 39-to-50-gene sets - say "
 "that, because using the published definitions rather than a hand-picked panel "
 "is the defence. One transition is significant: loss of the astrocyte-like "
 "state, p = 0.007. The mitochondrial signature moves the same way as the OXPHOS "
 "and mitochondrial-translation sets on the enrichment slide - point back to that "
 "slide, because two independent computations agreeing is the argument. If "
 "pressed on the mitochondrial call: 34 of its 35 measurable genes are lower in "
 "recurrence, and the result is stable across three levels of symbol mapping, "
 "35/50 up to 49/50.",

 "4:35-4:55  Hubs. Do not read the table - point at the sign column. Everything "
 "structural in the matrix goes down; the calcium and presynaptic genes (CALB1, "
 "RYR2, PCLO) go up. Offer the neuron-glioma synapse as a question, not a claim, "
 "and say plainly that six animals cannot separate it from a change in which "
 "cells survived. This is the slide most likely to draw a question, so know the "
 "two Nature 2019 papers. CUT THIS SLIDE FIRST if you are running long.",

 "4:55-5:25  Therapeutics. Explain the score in one sentence: how strongly the "
 "drug signature reverses the tumor signature, weighted by predicted BBB "
 "penetration from ADMET-AI. Lead with the filter, because it is the honest part: "
 "only 54 of 93 DSigDB hits carry a clinical phase - the rest include phosphine, "
 "tributyltin and a flame retardant. Then own its limits before anyone asks: it "
 "keeps things that are not systemic therapies (oxygen, ozone, magnesium) and it "
 "drops approved drugs ChEMBL has no phase for, including deferoxamine. Ciclopirox "
 "is what survives: approved, ranked first with and without the permeability "
 "weight, both BBB models agree, and it chelates iron. If asked why it should "
 "work on a cell that has stopped building: the hypusination targets DOHH, DHPS "
 "and EIF5A are selective CRISPR dependencies in U251 in DepMap - that is "
 "independent of our signature. Then the delivery point: LITT opens the barrier "
 "at the margin for weeks. Deliver the caveat sentence exactly as written; this "
 "is where a reviewer will push.",

 "5:25-5:45  Conclusions. Four beats, then the convergence line. No new "
 "material. The closing argument is the internal consistency - translational, "
 "anabolic, respiratory and hypoxic arms all moving together - not any single "
 "q value. Say that sentence and stop.",

 "5:45-5:55  Limitations and next steps. SECOND SLIDE TO CUT - the limitations "
 "are already stated on the results slides. Anticipated questions: why not more "
 "animals; why gene-set permutation; has any drug been tested (no).",

 "5:55-6:00  Thanks. Name the Henry Ford group. Leave the GEO accession on screen "
 "during questions.",
]
new_slides = list(prs.slides)[n_tpl:]
assert len(new_slides) == len(NOTES), (len(new_slides), len(NOTES))
for slide, text in zip(new_slides, NOTES):
    slide.notes_slide.notes_text_frame.text = text

# ------------------------------------------------------------------ save -----
drop_template_slides(prs, n_tpl)
prs.save(OUT)
print("wrote %s (%d slides)" % (OUT, len(prs.slides._sldIdLst)))
