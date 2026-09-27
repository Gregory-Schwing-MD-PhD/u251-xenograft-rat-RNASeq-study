#!/usr/bin/env python
"""Build the GEO metadata workbook for the eight EPIC arrays, templated on the RNA-seq submission already accepted.

The RNA submission is GSE338105 (public 17 Jul 2026). This script reuses its STUDY wording, its contributors, its
contact block and its protocols wherever they are the same experiment described twice, so the two GEO records read
as one study and can be joined sample by sample. It writes GEO/methylation/metadata_methylation.xlsx in the same
STUDY / SAMPLES / PROTOCOLS layout GEO's own template uses.

    python GEO/methylation/make_metadata.py

Nothing is invented here: the sample rows come from array_samplesheet.csv, the RNA cross-reference comes from
PIPELINE/00_fetch_reads/library_map.csv, and the reused wording is read out of the filled RNA workbook.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import openpyxl
from openpyxl.styles import Font

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
RNA_XLSX = REPO / "GEO" / "seq_template_u251_FILLED.xlsx"
SHEET = HERE / "array_samplesheet.csv"
MAP = REPO / "PIPELINE" / "00_fetch_reads" / "library_map.csv"
OUT = HERE / "metadata_methylation.xlsx"

PLATFORM = "Illumina Infinium MethylationEPIC BeadChip (EPIC v1.0)"
GPL = "GPL21145"          # Infinium MethylationEPIC v1.0 B5; confirm the exact B-version at submission
CHIP_SCANNED = "2026-04-07"


def rna_study_fields() -> dict:
    """The STUDY block of the accepted RNA submission: contributors, contact, and the protocols worth reusing."""
    wb = openpyxl.load_workbook(RNA_XLSX, data_only=True)
    ws = wb["Metadata"]
    out, key = {"contributor": []}, None
    for row in ws.iter_rows(min_row=1, max_row=200, values_only=True):
        cells = [c for c in row if c not in (None, "")]
        if not cells:
            continue
        k = str(cells[0]).strip().lstrip("*").lower()
        v = str(cells[1]).strip() if len(cells) > 1 else ""
        if k == "contributor" and v:
            out["contributor"].append(v)
        elif k in ("title", "summary (abstract)", "experimental design", "growth protocol",
                   "treatment protocol", "extract protocol"):
            out[k] = v
    return out


def main() -> int:
    if not RNA_XLSX.exists():
        raise SystemExit(f"missing the RNA submission workbook: {RNA_XLSX}")
    rna = rna_study_fields()
    rows = list(csv.DictReader(open(SHEET, encoding="utf-8-sig", newline="")))
    lib = {r["library"]: r for r in csv.DictReader(open(MAP, encoding="utf-8-sig", newline=""))}

    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Metadata"
    bold = Font(bold=True)

    def put(*cells, head=False):
        ws.append(list(cells))
        if head:
            for c in ws[ws.max_row]:
                c.font = bold

    put("# GEO metadata for the DNA methylation arm of the U251N LITT study.", head=True)
    put("# Built by GEO/methylation/make_metadata.py from array_samplesheet.csv and the accepted RNA submission.")
    put("# SUBMISSION STRATEGY: submit this as its own Series, then ask GEO to combine it with GSE338105 under a")
    put("#   SuperSeries. GEO's SuperSeries/SubSeries mechanism is how one study that used two assay types is held")
    put("#   together; a Series cannot mix an array platform with a sequencing platform.")
    put("")
    put("STUDY", head=True)
    put("*title", "DNA methylation profiling of U251N orthotopic glioblastoma xenografts in athymic RNU/RNU rats "
                 "before and after MRI-guided laser interstitial thermal therapy (LITT)")
    put("*summary (abstract)",
        "Companion DNA methylation dataset to GSE338105. Laser interstitial thermal therapy (LITT) is increasingly "
        "used for deep-seated and recurrent glioblastoma. Because the graft is human and the host is rat, every "
        "array probe can be assigned to tumour or host by sequence, so the tumour compartment of a regrowing lesion "
        "is measured rather than estimated. Eight Illumina Infinium MethylationEPIC arrays were run on one chip: the "
        "six orthotopic tumours that were also sequenced (three untreated, three regrown after LITT), the "
        "contralateral hemisphere of one tumour-bearing animal, and the parental U251N culture. The arrays support "
        "copy-number, MGMT-STP27 and tumour-fraction readouts, and are reported with the matched RNA libraries.")
    put("*experimental design",
        "Eight arrays on a single EPIC BeadChip (sentrix 205648300021, positions R01C01-R08C01), scanned "
        f"{CHIP_SCANNED}. Six tumours: primary/pre-LITT (IL67B, IL68B, IL69B) against recurrent/post-LITT (IL66B, "
        "IL70B, IL71B; the matched RNA libraries are named NL70B and NL71B). One contralateral hemisphere (N2; RNA "
        "library N269B) and one in vitro parental culture (C2B). Each array has a matched RNA-seq library in "
        "GSE338105; DNA and RNA were taken from different aliquots of the same piece. No probe reached genome-wide "
        "significance between arms; the arrays are reported for copy number, MGMT status and tumour fraction.")
    for c in rna.get("contributor", []):
        put("contributor", c)
    put("supplementary file", "U251N_EPIC_beta_noob.tsv.gz")
    put("supplementary file", "U251N_EPIC_detection_pvalues.tsv.gz")
    put("")
    put("SAMPLES", head=True)
    hdr = ["*Sample name", "*title", "*source name", "*organism", "characteristics: rat", "characteristics: arm",
           "characteristics: tissue", "characteristics: cell line", "characteristics: treatment",
           "characteristics: sentrix id", "characteristics: sentrix position",
           "characteristics: matched RNA library", "characteristics: matched RNA GEO sample",
           "*molecule", "*description", "*platform", "raw file: Grn", "raw file: Red", "processed data file"]
    put(*hdr, head=True)
    treat = {"Primary": "none (U251N implanted, not ablated)",
             "Recurrent": "MRI-guided laser interstitial thermal therapy (LITT), harvested at MRI-confirmed regrowth",
             "Contralateral": "none (U251N implanted in the opposite hemisphere)",
             "Culture": "none (in vitro reference)"}
    for r in rows:
        s, pos = r["sample"], r["sentrix_position"]
        rl = r["rna_library"]
        m = lib.get(rl, {})
        organism = "Homo sapiens; Rattus norvegicus" if r["arm"] != "Culture" else "Homo sapiens"
        src = ("in vitro cell culture" if r["arm"] == "Culture"
               else "brain, contralateral hemisphere" if r["arm"] == "Contralateral"
               else "brain (orthotopic xenograft)")
        put(f"{s}_EPIC",
            f"U251N {r['arm'].lower()} DNA methylation [{s}]",
            src, organism, r["rat"], r["arm"], r["tissue"], "U251N", treat.get(r["arm"], ""),
            r["sentrix_id"], pos, rl, m.get("geo_sample", "PENDING"),
            "genomic DNA",
            f"Matched RNA-seq library {rl} in GSE338105 ({m.get('geo_sample','PENDING')}); "
            f"{r['rna_human_pct']} % of that library's reads were assigned human.",
            f"{PLATFORM} ({GPL})",
            f"{r['sentrix_id']}_{pos}_Grn.idat", f"{r['sentrix_id']}_{pos}_Red.idat",
            "U251N_EPIC_beta_noob.tsv.gz")
    put("")
    put("PROTOCOLS", head=True)
    for k, label in (("growth protocol", "growth protocol"), ("treatment protocol", "treatment protocol")):
        if rna.get(k):
            put(label, rna[k])
    put("*extract protocol",
        "Genomic DNA was extracted from the same frozen pieces used for RNA (different aliquots). "
        "\\PENDING{extraction kit and operator from the laboratory record}")
    put("*label protocol", "Bisulfite conversion and Infinium chemistry as performed by the core facility. "
                           "\\PENDING{kit, input mass and conversion protocol from the core}")
    put("*hybridization protocol", "Illumina Infinium MethylationEPIC BeadChip, standard protocol. "
                                   "\\PENDING{core facility and instrument}")
    put("*scan protocol", f"Illumina iScan. Chip 205648300021 scanned {CHIP_SCANNED}.")
    put("*data processing",
        "IDATs were read with minfi; normalisation by noob; probes were filtered with the study's rat "
        "cross-hybridisation exclusion list before any human-compartment statistic. Copy number by conumee and "
        "DNAcopy against CopyNeutralIMA normals, with SeSAMe as an independent implementation; MGMT status by "
        "MGMT-STP27. The processed matrix supplied here is the noob beta matrix over all probes passing detection, "
        "with the matching detection P-value matrix.")
    put("*genome build/assembly", "hg19 (EPIC v1.0 manifest); copy-number segments lifted to hg38 where stated")
    put("*processed data files format and content",
        "U251N_EPIC_beta_noob.tsv.gz: probes (rows) by the eight arrays (columns), noob-normalised beta values. "
        "U251N_EPIC_detection_pvalues.tsv.gz: the same shape, detection P values.")
    put("")
    put("# The two file columns per sample are both required: GEO expects the Grn and Red IDAT of every array.")
    wb.save(OUT)
    print(f"wrote {OUT} ({len(rows)} samples)")
    print("PENDING fields to fill from the laboratory and core records: extraction kit, bisulfite kit and input, "
          "hybridisation facility, and the exact EPIC manifest B-version behind the platform accession.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
