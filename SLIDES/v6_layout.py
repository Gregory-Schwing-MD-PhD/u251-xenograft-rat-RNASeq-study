# -*- coding: utf-8 -*-
"""Slide layouts for deck v6: one panel per slide, the panel as big as the content area allows.

Greg, 2026-09-27: "a lot of the slides have figure too small and too much text and I can't interpret the plots since
they are 4 subplots - it needs to take up more of the slide". So v6 has one rule: every results slide carries ONE
single-panel figure, sized to the content area, with at most three short lines under it. The text height is what is
left over, not what the text feels like taking.

Content area of the CNS template (cns_template): x from L to R (10,150,000 EMU), y from TOP to BOTTOM
(4,830,000 EMU), with the reference line below at REFS_TOP.
"""
from __future__ import annotations

from cns_template import (BOTTOM, CAP, GAP, GREY, L, R, TOP, add_text, content_slide, fit, img_size, navy_bar,
                          notes, picture, refs)

GAP_X = GAP

PAD = 90_000          # between the figure and the first line of text
TXT = 18              # Greg's rule: slide text at 18-20 pt

# Reserving by paragraph count puts the last WRAPPED line through the reference footer, which is what v5 did.
# These constants are measured off the CNS template's own rendering at 18 pt across the content width.
EMU_PER_PT_CHAR = 6600      # 10,150,000 EMU / (18 pt) / 6,600 = about 85 characters on a line
LINE_H_PER_PT = 16_000      # 18 pt line, 1.2 spacing, in EMU
PT = 12_700


def _text(line):
    """The plain text of a line, which may be a list of (text, colour) runs."""
    if isinstance(line, (list, tuple)) and line and isinstance(line[0], (list, tuple)):
        return "".join(t for t, _ in line)
    return str(line)


def text_height(lines, width, size=TXT, gap=10):
    """Height the lines will actually occupy once wrapped, in EMU."""
    if not lines:
        return 0
    per_line = max(20, width / (size * EMU_PER_PT_CHAR))
    visual = sum(max(1, -(-len(_text(l)) // int(per_line))) for l in lines)
    return int(visual * size * LINE_H_PER_PT + (len(lines) - 1) * gap * PT + 0.35 * size * LINE_H_PER_PT)


TEXT_BOX = {0: 0, 1: 620_000, 2: 1_020_000, 3: 1_420_000}   # kept for callers that ask by count


def panel_slide(prs, title, image, lines=(), ref=None, note=None, size=TXT, gap=10, title_size=30):
    """One figure, as large as the content area allows, with at most three lines under it.

    `lines` may hold plain strings or the (text, colour) pairs cns_template.add_text accepts.
    """
    lines = list(lines)
    if len(lines) > 3:
        raise ValueError(f"{title!r}: {len(lines)} lines; v6 allows three")
    s = content_slide(prs, title, size=title_size)
    reserve = text_height(lines, R - L, size, gap) + (PAD if lines else 0)
    w, h = fit(image, R - L, BOTTOM - TOP - reserve)
    picture(s, image, L + (R - L - w) // 2, TOP, w, h)
    if lines:
        top = TOP + h + PAD
        add_text(s, L, top, R - L, BOTTOM - top, lines, size=size, gap=gap)
    if ref:
        refs(s, ref)
    if note:
        notes(s, note)
    return s


def side_panel_slide(prs, title, image, lines=(), ref=None, note=None, size=TXT, gap=14, img_share=0.62):
    """A squarish panel at full content height on the left, the reading beside it.

    A panel whose aspect is near 1:1 can only reach about 57 % of the content area when it sits above the text,
    because the area is 2.1:1. Put it at the left instead and it keeps the full height.
    """
    s = content_slide(prs, title)
    w, h = fit(image, int((R - L) * img_share), BOTTOM - TOP)
    picture(s, image, L, TOP + (BOTTOM - TOP - h) // 2, w, h)
    x2 = L + w + GAP_X
    if lines:
        add_text(s, x2, TOP + 120_000, R - x2, BOTTOM - TOP - 120_000, list(lines), size=size, gap=gap)
    if ref:
        refs(s, ref)
    if note:
        notes(s, note)
    return s


def auto_slide(prs, title, image, lines=(), ref=None, note=None, size=TXT, **kw):
    """Pick the layout that gives the panel the most area: wide panels above the text, squarish panels beside it."""
    w0, h0 = img_size(image)
    wide = (w0 / h0) >= 1.85
    fn = panel_slide if wide else side_panel_slide
    return fn(prs, title, image, lines=lines, ref=ref, note=note, size=size, **kw)


def statement_slide(prs, title, lines, bar=None, ref=None, note=None, size=20, gap=18):
    """No figure: a few lines and, optionally, one highlighted statement in the navy bar."""
    s = content_slide(prs, title)
    bar_h = 1_250_000 if bar else 0
    add_text(s, L, TOP + 100_000, R - L, BOTTOM - TOP - bar_h - 200_000, lines, size=size, gap=gap)
    if bar:
        navy_bar(s, L, BOTTOM - bar_h, R - L, bar_h, bar, size=19)
    if ref:
        refs(s, ref)
    if note:
        notes(s, note)
    return s


def panel_fit(image, n_lines=2, side=None):
    """(area share, fill) for the layout `auto_slide` would choose.

    `fill` is how much of the limiting dimension the panel uses. A portrait panel cannot cover more than about 48 %
    of a 2.1:1 content box however it is placed, so area alone would condemn a panel that is already as large as the
    slide allows; `fill` near 1.0 says the panel is limited by the slide, not by the layout.
    """
    area = (R - L) * (BOTTOM - TOP)
    w0, h0 = img_size(image)
    if side is None:
        side = (w0 / h0) < 1.85
    if side:
        box_w, box_h = int((R - L) * 0.62), BOTTOM - TOP
    else:
        reserve = (TEXT_BOX[n_lines] if not isinstance(n_lines, (list, tuple))
                   else text_height(list(n_lines), R - L)) + (PAD if n_lines else 0)
        box_w, box_h = R - L, BOTTOM - TOP - reserve
    w, h = fit(image, box_w, box_h)
    return (w * h) / area, max(w / box_w, h / box_h)


def figure_area_fraction(image, n_lines=2, side=None):
    """The share of the content area the panel occupies under the layout `auto_slide` would choose.

    The build prints this for every slide: a results slide whose panel covers less than about 45 % of the content
    area is the defect Greg reported, and the builder refuses to pass it silently.
    """
    area = (R - L) * (BOTTOM - TOP)
    w0, h0 = img_size(image)
    if side is None:
        side = (w0 / h0) < 1.85
    if side:
        w, h = fit(image, int((R - L) * 0.62), BOTTOM - TOP)
    else:
        reserve = TEXT_BOX[n_lines] + (PAD if n_lines else 0)
        w, h = fit(image, R - L, BOTTOM - TOP - reserve)
    return (w * h) / area
