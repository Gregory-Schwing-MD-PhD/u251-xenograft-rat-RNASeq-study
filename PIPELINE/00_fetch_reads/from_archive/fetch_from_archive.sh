#!/usr/bin/env bash
# Brings the U251 work's raw reads back from the public archive into a folder of your choice, every file checked
# against the md5 the archive publishes, so the work can start again from GEO and nothing local.
#
#   bash fetch_from_archive.sh DEST [PATTERN]
#
# PATTERN is an extended regular expression matched against the manifest's study, run and library columns
# (default: every row), for example 'PRJNA1492510' for the ten U251 libraries or 'SRR1169975' for one run.
#
# The fetch is all or nothing and sized before any byte lands: each file's size is read from the archive (an HTTP
# HEAD) and must equal the manifest's, then the whole is admitted or refused. On Wayne State's grid a DEST on the lab
# volume (/rs/...) is admitted by the owner's size rule, the size rule of the lab grid repository (grid/size-rule.sh) (admit_fetch): run it as a
# Slurm job with SIZE_RULE set to that file and DATA to the data root the rule keeps its admissions under. Anywhere
# else the free space df reports for DEST must hold what is still to land. A file already in DEST with the right md5
# is kept, so a cut fetch is simply run again.
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MANIFEST=${MANIFEST:-$here/manifest.tsv}
dest=${1:?usage: fetch_from_archive.sh DEST [PATTERN]}
pattern=${2:-.}
mkdir -p "$dest"
todo=$(mktemp)
trap 'rm -f "$todo"' EXIT

md5_of() { md5sum < "$1" | cut -d' ' -f1; }

archive_bytes() {  # archive_bytes URL: the byte count the archive states for URL, read before any byte lands
  curl -sSfIL --retry 3 --max-time 120 "$1" | tr -d '\r' | awk 'tolower($1)=="content-length:" {n=$2} END {print n}'
}

size_everything() {  # writes the rows still to land into $todo; prints "<need> <total>"
  local need=0 total=0 study run library layout file url bytes md5 rest stated
  while IFS=$'\t' read -r study run library layout file url bytes md5 rest; do
    total=$((total + bytes))
    if [ -f "$dest/$file" ] && [ "$(md5_of "$dest/$file")" = "$md5" ]; then continue; fi
    stated=$(archive_bytes "$url")
    [ "$stated" = "$bytes" ] || { echo "refused: $file is $stated bytes at the archive, $bytes in the manifest" >&2; return 4; }
    need=$((need + bytes))
    printf '%s\t%s\t%s\t%s\n' "$file" "$url" "$bytes" "$md5" >> "$todo"
  done < <(awk -F'\t' -v p="$pattern" 'NR > 1 && ($1 ~ p || $2 ~ p || $3 ~ p)' "$MANIFEST")
  [ "$total" -gt 0 ] || { echo "refused: no manifest row matches $pattern" >&2; return 2; }
  echo "$need $total"
}

on_lab_volume() { case "$(readlink -f "$dest")" in /rs/*|/.rs/*) return 0 ;; *) return 1 ;; esac; }

admit() {  # admit NEED TOTAL: 0 when the whole fetch may land, 4 when it is refused
  local need=$1 total=$2 free
  if on_lab_volume; then
    [ -f "${SIZE_RULE:-}" ] && [ -n "${DATA:-}" ] || { echo "refused: DEST is on the lab volume; set SIZE_RULE and DATA" >&2; return 4; }
    [ -n "${SLURM_JOB_ID:-}" ] || { echo "refused: a fetch onto the lab volume runs as a Slurm job" >&2; return 4; }
    ( set +eu; . "$SIZE_RULE"; admit_fetch "${FETCH_ID:-u251_archive}" "$total" "$need" "u251 reads from the archive" )
    return
  fi
  free=$(df -B1 --output=avail "$dest" | tail -1 | tr -d ' ')
  [ "$free" -ge "$need" ] || { echo "refused: $need bytes to land, $free free at $dest" >&2; return 4; }
}

land_one() {  # land_one FILE URL BYTES MD5: the file fetched and kept only when its md5 is the archive's
  local file=$1 url=$2 bytes=$3 md5=$4 got
  curl -sSfL --retry 5 --retry-delay 30 -C - -o "$dest/$file.part" "$url"
  got=$(md5_of "$dest/$file.part")
  if [ "$got" != "$md5" ]; then
    mv "$dest/$file.part" "$dest/$file.bad"
    echo "md5 differs: $file $got, archive $md5" >&2
    return 1
  fi
  mv "$dest/$file.part" "$dest/$file"
  printf '%s\t%s\t%s\t%s\n' "$file" "$bytes" "$md5" "$(date -u +%FT%TZ)" >> "$dest/FETCHED.tsv"
}

read -r need total < <(size_everything) || exit $?
echo "to land: $need of $total bytes into $dest"
[ "$need" -eq 0 ] || admit "$need" "$total" || exit 4
failed=0
while IFS=$'\t' read -r file url bytes md5; do
  land_one "$file" "$url" "$bytes" "$md5" || failed=$((failed + 1))
done < "$todo"
[ "$failed" -eq 0 ] || { echo "$failed file(s) did not match the archive's md5" >&2; exit 1; }
echo "every file matches the archive's md5"
