#!/usr/bin/env python
"""Build the GEO metadata workbook for the eight EPIC arrays, templated on the RNA-seq submission already accepted.

The RNA submission is GSE338105 (public 17 Jul 2026). This script reuses its STUDY wording, its contributors, its
contact block and its protocols wherever they are the same experiment described twice, so the two GEO records read
as one study and can be joined sample by sample. It writes GEO/methylation/metadata_methylation.xlsx in the same
STUDY / SAMPLES / PROTOCOLS layout GEO's own template uses.

    python GEO/methylation/make_metadata.py

Nothing is invented here: the sample rows come from array_samplesheet.csv, the RNA cross-reference comes from
PIPELINE/00_fetch_reads/library_map.csv, and the reused wording is read out of the filled RNA workbook.

The array names are NOT changed to match the RNA names, and the RNA names are not changed either: the RNA names are
public and minted in SRA, and the array names are what the chip, the core's sheet and every analysis in this
repository use. Three pieces therefore carry two names each (NL70B/IL70B, NL71B/IL71B, N269B/N2). Every one of the
three is stated in the sample title, in four characteristics fields, and in the deposited
GEO/multiomics_sample_key.tsv, so no reader has to infer a pairing from a name. This matters because GEO mints a
BioSample for a sequencing sample but not for an array sample (checked against real multi-omics SuperSeries), so NCBI
will not link the two arms for anyone: the link exists only in the metadata written here.
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
GPL = "GPL21145"          # Infinium MethylationEPIC v1.0; the core annotated with ilm10b4 (manifest B4)
# 2022, not 2026. The RunInfo block of every IDAT records Scan on 4/7/2022 09:24-09:26 and Decode on 07/11/2021, and
# ANALYSIS/SAMPLE_KEY.md says the same. An earlier version of this file had 2026 and it reached the workbook.
CHIP_SCANNED = "2022-04-07"
CHIP_DECODED = "2021-07-11"


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
    put("#   together, one SubSeries per assay.")
    put("# SAMPLE NAMES: three pieces carry two names, because two facilities named the same animals differently.")
    put("#   The array side (this submission, the chip, the core's sheet) says IL70B, IL71B, N2. The RNA side")
    put("#   (GSE338105, SRA, the sequencing core's file names) says NL70B, NL71B, N269B. Neither is renamed. The")
    put("#   pairing is stated per sample below and in the deposited multiomics_sample_key.tsv.")
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
        "copy-number, MGMT-STP27 and tumour-fraction readouts, and are reported with the matched RNA libraries. "
        "Every array here has a matched RNA-seq library in GSE338105, taken from a different aliquot of the same "
        "piece of tissue; the pairing is given per sample and in the supplementary file multiomics_sample_key.tsv. "
        "Three pieces carry a different name in each arm because two facilities named the same animals differently: "
        "arrays IL70B, IL71B and N2 are RNA libraries NL70B, NL71B and N269B respectively. The remaining five "
        "(IL66B, IL67B, IL68B, IL69B, C2B) carry the same name in both arms. Two RNA libraries in GSE338105 have no "
        "array in this Series (IL64B and N168B).")
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
    put("supplementary file", "multiomics_sample_key.tsv")
    put("")
    put("SAMPLES", head=True)
    hdr = ["*Sample name", "*title", "*source name", "*organism", "characteristics: rat", "characteristics: arm",
           "characteristics: tissue", "characteristics: cell line", "characteristics: treatment",
           "characteristics: sentrix id", "characteristics: sentrix position",
           "characteristics: matched RNA library", "characteristics: matched RNA GEO sample",
           "characteristics: matched RNA SRA experiment", "characteristics: matched RNA BioSample",
           "characteristics: matched RNA series", "characteristics: array label on the core sheet",
           "characteristics: plate well", "characteristics: dna volume ul",
           "characteristics: qc ct HB-313", "characteristics: qc ct HB-365",
           "characteristics: detected probes at P<0.05", "characteristics: detected probes percent",
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
        # The array name and the RNA name differ for three of the eight pieces. Say so in the title, so that a
        # reader who only ever sees the sample list still sees the pairing.
        renamed = rl != s
        title = (f"U251N {r['arm'].lower()} DNA methylation [{s} = RNA library {rl}]" if renamed
                 else f"U251N {r['arm'].lower()} DNA methylation [{s}]")
        same = ("DNA and RNA come from the same culture"
                if r["arm"] == "Culture" else
                "DNA and RNA come from different aliquots of the same piece of tissue")
        alias = (f" This array is named {s} on the core's sample sheet and its RNA library is named {rl}: "
                 f"the same piece under two names, because the array core and the sequencing core named the "
                 f"animals differently. Neither name has been changed." if renamed else "")
        put(f"{s}_EPIC",
            title,
            src, organism, r["rat"], r["arm"], r["tissue"], "U251N", treat.get(r["arm"], ""),
            r["sentrix_id"], pos, rl, m.get("geo_sample", "PENDING"),
            m.get("experiment_accession", "PENDING"), m.get("sample_accession", "PENDING"), "GSE338105", s,
            r["well"], r["dna_volume_ul"], r["ct_hb313"], r["ct_hb365"], r["detected_probes"], r["detected_pct"],
            "genomic DNA",
            f"{same}. Matched RNA-seq library {rl} in GSE338105 ({m.get('geo_sample','PENDING')}, "
            f"{m.get('experiment_accession','PENDING')}, BioSample {m.get('sample_accession','PENDING')}).{alias} "
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
        "Genomic DNA was extracted by the laboratory from the same frozen pieces of tissue used for RNA, as separate "
        "aliquots. 45 uL of DNA per sample was submitted to the array core with no added water (plate 1724, one "
        "column of eight wells, A01 to H01). "
        "\\PENDING{the DNA extraction kit, and the DNA mass corresponding to that 45 uL, from the laboratory record}")
    put("*label protocol",
        "Bisulfite conversion and Infinium chemistry were performed by the array core on plate 1724 (\"Raj-8\"). The "
        "core checked each sample's DNA by qPCR before conversion; the per-sample Ct values for assays HB-313 and "
        "HB-365 are given in the sample characteristics. "
        "\\PENDING{the bisulfite conversion kit and the DNA input mass, from the core}")
    put("*hybridization protocol",
        "Illumina Infinium MethylationEPIC BeadChip v1.0, 8x5 format (SentrixDescriptor \"BeadChip 8x5\"); one chip, "
        "205648300021, carrying the eight samples in column 1. Standard Infinium protocol. "
        "\\PENDING{the name of the array core facility}")
    put("*scan protocol",
        f"Illumina iScan, scanner N141, iScan Control Software 3.4.8 (FPGA 4.0.20). Chip 205648300021 was decoded on "
        f"{CHIP_DECODED} and scanned on {CHIP_SCANNED} between 09:24 and 09:26; the instrument software extracted "
        "intensities with Extract Algorithm StandardWithBackground. Every value in this sentence is read from the "
        "RunInfo block of the deposited IDAT files themselves.")
    put("*data processing",
        "The deposited matrices were built from the deposited IDATs with minfi 1.56.0 and "
        "IlluminaHumanMethylationEPICmanifest: detection P values by detectionP, normalisation by preprocessNoob, "
        "beta values by getBeta, over all 866,238 probes on the array. Column names are the array sample names. "
        "Per array, 864,815 to 865,311 probes pass detection at P < 0.05 (99.879 to 99.937 %), as reported by the "
        "core and given per sample in the characteristics. In the study's own analyses, probes were then filtered "
        "with a rat cross-hybridisation exclusion list before any human-compartment statistic; copy number by "
        "conumee and DNAcopy against CopyNeutralIMA normals, with SeSAMe as an independent implementation; MGMT "
        "status by MGMT-STP27. sessionInfo_processed.txt records the exact versions used for the deposited matrices.")
    put("*genome build/assembly",
        "hg19. The array core's own analysis used the EPIC v1.0 B4 manifest "
        "(IlluminaHumanMethylationEPICanno.ilm10b4.hg19); copy-number segments are lifted to hg38 where stated.")
    put("*processed data files format and content",
        "U251N_EPIC_beta_noob.tsv.gz: probes (rows) by the eight arrays (columns), noob-normalised beta values. "
        "U251N_EPIC_detection_pvalues.tsv.gz: the same shape, detection P values.")
    put("")
    put("# The two file columns per sample are both required: GEO expects the Grn and Red IDAT of every array.")
    wb.save(OUT)
    print(f"wrote {OUT} ({len(rows)} samples)")
    print("Three fields still need a human:")
    print("  1. the DNA extraction kit, and the DNA mass behind the 45 uL   -- the laboratory")
    print("  2. the bisulfite conversion kit and the DNA input mass         -- the array core")
    print("  3. the name of the array core facility                         -- either")
    print("Everything else is filled. The scan protocol and the manifest version were read out of the data itself,")
    print("not asked for: the IDAT RunInfo block and the core's own dnamet.R.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
