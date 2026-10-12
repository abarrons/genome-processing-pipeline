#!/bin/bash
## Consolidate filtered assemblies from all 3 batches into one folder for downstream analysis
## and NCBI submission.
## Detects duplicate sampleIDs ACROSS batches (does NOT silently overwrite) and
## writes a duplicates report so you can decide which copy to use/submit.
##
## Usage:  bash 3.0_consolidate_assemblies.sh [OUTPUT_DIR]

set -euo pipefail

BASE=/dfs10/bio/abarrons/CA_curto/isolates
OUTDIR="${1:-/dfs10/bio/abarrons/CA_curto/NCBI_submission}"

DEST="$OUTDIR/assemblies"
DUPDIR="$OUTDIR/duplicates"
REPORT="$OUTDIR/consolidation_report.tsv"

mkdir -p "$DEST" "$DUPDIR"
: > "$REPORT"
printf "sampleID\tbatch\tsource_path\tstatus\n" >> "$REPORT"

for b in 1stbatch 2ndbatch 3rdbatch; do
  srcdir="$BASE/$b/spades_out/filtered_assemblies"
  [ -d "$srcdir" ] || { echo "WARN: no filtered_assemblies in $b" >&2; continue; }
  for f in "$srcdir"/*_filtered_scaffolds.fasta; do
    [ -e "$f" ] || continue
    base=$(basename "$f")
    sample=${base%_filtered_scaffolds.fasta}
    target="$DEST/$base"
    if [ -e "$target" ]; then
      ## Collision: this sampleID already came from an earlier batch.
      ## Keep BOTH by stashing the new one under a batch-tagged name in DUPDIR.
      dupname="${sample}__${b}_filtered_scaffolds.fasta"
      cp "$f" "$DUPDIR/$dupname"
      printf "%s\t%s\t%s\tDUPLICATE (kept in duplicates/ as %s)\n" \
             "$sample" "$b" "$f" "$dupname" >> "$REPORT"
      echo "DUPLICATE sampleID: $sample also in $b -> stashed as $dupname"
    else
      cp "$f" "$target"
      printf "%s\t%s\t%s\tconsolidated\n" "$sample" "$b" "$f" >> "$REPORT"
    fi
  done
done

echo
echo "Consolidated: $(ls "$DEST"/*_filtered_scaffolds.fasta 2>/dev/null | wc -l) genomes in $DEST"
ndup=$(ls "$DUPDIR"/*.fasta 2>/dev/null | wc -l)
echo "Duplicates flagged: $ndup (in $DUPDIR)"
echo "Report: $REPORT"
if [ "$ndup" -gt 0 ]; then
  echo
  echo "== Duplicate sampleIDs to resolve =="
  awk -F'\t' '$4 ~ /DUPLICATE/ {print "  "$1" (from "$2")"}' "$REPORT"
  echo "The FIRST-seen copy (1stbatch>2ndbatch>3rdbatch order) is in assemblies/;"
  echo "the colliding copy is in duplicates/. Decide which to keep before submission."
fi
