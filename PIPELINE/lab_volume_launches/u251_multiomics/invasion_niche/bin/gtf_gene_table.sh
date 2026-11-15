#!/bin/bash
# gene_id, chr, start, end, strand, gene_type, gene_name for every "gene" line of a GTF (plain or .gz), plus the
# GTF's header lines. Idempotent: the table is rebuilt only when the source's md5 differs from the one recorded.
#   gtf_gene_table.sh <gtf> <out.tsv>
set -eu
GTF=$1; OUT=$2
SRC_MD5=$(md5sum "$GTF" | cut -d' ' -f1)
if [ -s "$OUT" ] && [ -s "$OUT.src_md5" ] && [ "$(cat "$OUT.src_md5")" = "$SRC_MD5" ]; then
  echo "gene table up to date: $OUT (source md5 $SRC_MD5)"; exit 0
fi
case "$GTF" in *.gz) CAT="zcat" ;; *) CAT="cat" ;; esac
mkdir -p "$(dirname "$OUT")"
$CAT "$GTF" | awk -F'\t' '/^#/ { if (nh < 20) print > "/dev/stderr"; nh++; next } { exit }' 2> "$OUT.header.part" || true
$CAT "$GTF" | awk -F'\t' '
  BEGIN { OFS = "\t"; print "gene_id", "chr", "start", "end", "strand", "gene_type", "gene_name" }
  $3 == "gene" {
    id = ""; ty = ""; nm = ""
    n = split($9, a, ";")
    for (i = 1; i <= n; i++) {
      s = a[i]; sub(/^ +/, "", s)
      if (s ~ /^gene_id /)           { id = s; sub(/^gene_id "/, "", id); sub(/"$/, "", id) }
      else if (s ~ /^gene_type /)    { ty = s; sub(/^gene_type "/, "", ty); sub(/"$/, "", ty) }
      else if (s ~ /^gene_biotype /) { if (ty == "") { ty = s; sub(/^gene_biotype "/, "", ty); sub(/"$/, "", ty) } }
      else if (s ~ /^gene_name /)    { nm = s; sub(/^gene_name "/, "", nm); sub(/"$/, "", nm) }
    }
    print id, $1, $4, $5, $7, ty, nm
  }' > "$OUT.part"
n=$(($(wc -l < "$OUT.part") - 1))
dups=$(cut -f1 "$OUT.part" | sort | uniq -d | wc -l)
if [ "$n" -lt 1000 ] || [ "$dups" -ne 0 ]; then echo "gene table failed: $n genes, $dups duplicated ids"; exit 2; fi
mv "$OUT.header.part" "$OUT.header"
mv "$OUT.part" "$OUT"
echo "$SRC_MD5" > "$OUT.src_md5"
echo "gene table written: $OUT ($n genes; source md5 $SRC_MD5)"
