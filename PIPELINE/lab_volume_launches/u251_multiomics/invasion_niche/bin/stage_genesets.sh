#!/bin/bash
# Copy every gene-set source file into results/genesets/ and write genesets/SOURCES.tsv (file, where it was copied
# from, source URL, md5, note). Nothing is typed in: the sets are built from these files by inv_stepC.R.
#   stage_genesets.sh <ivygap matrix> <spread signatures dir> <genesets_src dir> <hallmark gmt> <c2 gmt> <dest dir>
set -eu
IVY=$1; SPREAD=$2; SRC=$3; GMT_H=$4; GMT_C2=$5; DEST=$6
mkdir -p "$DEST"
S="$DEST/SOURCES.tsv.part"
printf 'file\tcopied_from\turl\tmd5\tnote\n' > "$S"
add() { # dest-name source url note
  cp -f "$2" "$DEST/$1"
  printf '%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$(md5sum "$DEST/$1" | cut -d' ' -f1)" "$4" >> "$S"
}
add ivygap2018_anatomic4.txt "$IVY" "https://www.science.org/doi/suppl/10.1126/science.aaf2666/suppl_file/aaf2666_table-s14.xlsx" \
  "Puchalski 2018 Science 360:660 (PMID 29748285), Table S14 sheet ciberGeneSignature_rerun, converted by ANALYSIS/cibersort/signatures_online (see its MANIFEST.md)"
for c in LE CT CTmvp CTpan; do
  [ -f "$SPREAD/ivygap_$c.txt" ] && add "spread_ivygap_$c.txt" "$SPREAD/ivygap_$c.txt" "derived" "frozen before any score (ANALYSIS/spread/signatures/MANIFEST.md, sha256 there); used only as a cross-check"
done
for s in MES AC OPC NPC; do
  [ -f "$SPREAD/neftel_$s.txt" ] && add "spread_neftel_$s.txt" "$SPREAD/neftel_$s.txt" "derived" "argmax of the repository's CIBERSORT signature built from the Neftel 2019 SS2 cells in U/ANALYSIS/refs/neftel_2019 (not the published meta-modules)"
done
[ -f "$SPREAD/MANIFEST.md" ] && add spread_signatures_MANIFEST.md "$SPREAD/MANIFEST.md" "derived" "rules and sha256 of the spread signatures"
# publisher supplements fetched on 2026-09-27 (scratchpad invasion_fetch, convert_supp.py); md5 checked against the fetch record
while IFS=$'\t' read -r f url md5 sheet rows dropped; do
  [ "$f" = "file" ] && continue
  if [ -f "$SRC/$f" ]; then
    got=$(md5sum "$SRC/$f" | cut -d' ' -f1)
    if [ "$got" != "$md5" ]; then echo "md5 mismatch for $f: $got vs fetch record $md5"; exit 3; fi
    add "$f" "$SRC/$f" "$url" "sheet '$sheet'; rows converted $rows; dropped cells: $dropped"
  fi
done < "$SRC/SOURCES.tsv"
add neftel2019_tableS2_metamodules.tsv "$SRC/neftel2019_tableS2_metamodules.tsv" \
  "https://ars.els-cdn.com/content/image/1-s2.0-S0092867419306877-mmc2.xlsx" "Neftel 2019 Cell 178:835, Table S2 meta-module gene lists, converted cell by cell by convert_supp.py"
add venkataramani2022_tableS2a_invasivity.tsv "$SRC/venkataramani2022_tableS2a_invasivity.tsv" \
  "https://ars.els-cdn.com/content/image/1-s2.0-S0092867422008479-mmc2.xlsx" "Venkataramani 2022 Cell 185:2899 (PMID 35914528), Table S2a invasivity score genes, converted by convert_supp.py"
add convert_supp.py "$SRC/convert_supp.py" "local" "the converter (openpyxl, no hand edits)"
add hs.hallmark.v2023.2.gmt "$GMT_H" "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/" "MSigDB 2023.2 hallmark (fallback and cross-check for msigdbr)"
add hs.c2_curated.v2023.2.gmt "$GMT_C2" "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/" "MSigDB 2023.2 C2 (fallback and cross-check for msigdbr)"
mv "$S" "$DEST/SOURCES.tsv"
echo "gene-set sources staged: $(($(wc -l < "$DEST/SOURCES.tsv") - 1)) files in $DEST"
