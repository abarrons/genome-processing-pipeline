#!/bin/bash

#SBATCH --job-name=btk_filter       ## Name of the job.
#SBATCH -A JMARTINY_LAB             ## Account to charge; personal or lab
#SBATCH -p standard                 ## Partition/queue name
#SBATCH --nodes=1                   ## number of nodes
#SBATCH --ntasks=1                  ## number of processes
#SBATCH --cpus-per-task=8           ## CPUs per task
#SBATCH -o btkfilt_%A_%a.out        ## STDOUT; %A=jobid, %a=array index
#SBATCH -e btkfilt_%A_%a.err        ## STDERR; %A=jobid, %a=array index
## NOTE: --array is set at submission time by submit_btkfilter.sh (auto-counted).

## BlobToolKit (BlobTools2) per-isolate QC + filter, with a
## Curtobacterium-conditional GC filter.
##
## Flow (v2):
##   create BlobDir -> add cov -> add hits(+taxdump) -> view (plots)
##   -> genus check -> filter:
##        length + cov ALWAYS applied
##        gc applied ONLY if isolate is predominantly Curtobacterium (genus)
##
## Requires: taxdump downloaded via 00_setup_taxdump.sh.
## (Completeness/contamination via CheckM is handled as a separate step.)
## BATCH passed in via --export by the launcher.

## ------------------------------------------------------------
: "${BATCH:?BATCH not set -- submit via ./submit_btkfilter.sh, not sbatch directly}"

## Initialize conda (cluster-specific: module + init script)
module load miniconda3/25.11.1
. ~/.mycondainit-25.11.1

BASE=/dfs10/bio/abarrons/CA_curto/isolates/${BATCH}
ASSEMBLIES=$BASE/spades_out/assemblies
BAMDIR=$BASE/bowtie2_output
BLASTDIR=$BASE/blast_out
BLOBDIRS=$BASE/blobtools_out/blobdirs
PLOTS=$BASE/blobtools_out/plots
TABLES=$BASE/blobtools_out/tables
FILTERED=$BASE/spades_out/filtered_assemblies

TAXDUMP=/dfs10/bio/abarrons/databases/taxdump

## Filter thresholds (v2 field names; GC is a proportion 0-1)
## NOTE: MINCOV is a light low-coverage floor to drop near-zero-cov artifacts.
## These isolates sequence at ~single-digit depth (~7x on real Curto contigs),
## so cov=3 removes junk without discarding the real genome. (The old value
## of 30, carried from a different coverage metric, wrongly cut real contigs.)
MINLEN=5000
MINGC=0.5
MINCOV=3

## Genus considered the "clean" target (gets the GC filter when dominant)
TARGET_GENUS=Curtobacterium

THREAD=8

mkdir -p "$BLOBDIRS" "$PLOTS" "$TABLES" "$FILTERED"

## taxdump sanity
if [ ! -d "$TAXDUMP" ] || [ -z "$(ls -A "$TAXDUMP" 2>/dev/null)" ]; then
  echo "ERROR: taxdump not found/empty at $TAXDUMP -- run 00_setup_taxdump.sh first" >&2
  exit 1
fi

## ---- Per-task isolate selection ----
cd "$ASSEMBLIES"
FASTA=$(ls *_scaffolds.fasta | head -n "$SLURM_ARRAY_TASK_ID" | tail -n 1)
ISOLATE_ID=${FASTA%_scaffolds.fasta}
ASM="$ASSEMBLIES/${ISOLATE_ID}_scaffolds.fasta"
DS="$BLOBDIRS/${ISOLATE_ID}"

echo "[$(date)] Batch=$BATCH Task=$SLURM_ARRAY_TASK_ID Isolate=$ISOLATE_ID"

## ---- conda env (conda already initialized at top via module + init script) ----
conda activate btk

## ================= STAGE 1: build BlobDir =================
## IMPORTANT: BlobTools2 will NOT overwrite existing fields; a stale BlobDir
## from a previous run makes create/add silently skip (cov, hits, etc. never
## get built). Remove any prior BlobDir for this isolate so it is built fresh.
## The BlobDir is derived data only -- assemblies/BAMs/blast are untouched.
rm -rf "$DS"

blobtools create --fasta "$ASM" "$DS"

## blobtools' cov parser expects (and then deletes) a CSI-format BAM index.
## Our BAMs were indexed as .bai by samtools, so blobtools crashes trying to
## remove a nonexistent .sorted.bam.csi. Pre-build the .csi so it exists.
## Needs samtools; the btk env usually has it, else load the module.
BAM="$BAMDIR/${ISOLATE_ID}.sorted.bam"
if command -v samtools >/dev/null 2>&1; then
  samtools index -c "$BAM"        # -c writes <bam>.csi
else
  module load samtools/1.15.1 && samtools index -c "$BAM"
fi

## Add coverage with an explicit library name so the resulting field is
## predictably 'reads_cov' (blobtools appends _cov to the given name),
## rather than depending on the BAM basename.
blobtools add --replace --cov "${BAM}=reads" --threads $THREAD "$DS"

blobtools add --replace \
  --hits "$BLASTDIR/${ISOLATE_ID}.blast.out" \
  --taxrule bestsumorder \
  --taxdump "$TAXDUMP" \
  "$DS"

## ================= STAGE 2: plots =================
blobtools view --plot --out "$PLOTS" "$DS" \
  || echo "WARN: plot step failed for $ISOLATE_ID (filtering continues)"

## ================= STAGE 3: genus check =================
## Dump a table with the genus call AND length per contig, then measure the
## fraction of classified assembly LENGTH (bp) assigned to TARGET_GENUS.
## Length-weighting is more robust than contig counts: many tiny contaminant
## contigs won't outvote one large Curto contig.
## ================= STAGE 3: dominant-genus check =================
## Identify the DOMINANT genus by summed contig LENGTH (bp) among classified
## contigs (ignoring no-hit). Length-weighting is robust to many tiny contigs.
## The dominant genus becomes the target to KEEP; everything else is dropped.
GENUS_TABLE="$TABLES/${ISOLATE_ID}.genus.tsv"
blobtools filter \
  --table "$GENUS_TABLE" \
  --table-fields bestsumorder_genus,length \
  "$DS"

DOMINANT_GENUS=""
OTHER_GENERA=""
if [ -f "$GENUS_TABLE" ]; then
  ## Find genus + length columns by header, sum length per genus, report the max.
  DOMINANT_GENUS=$(awk '
    NR==1 {
      for (i=1;i<=NF;i++) {
        if ($i ~ /genus/)  gcol=i
        if ($i == "length") lcol=i
      }
      next
    }
    {
      g=$gcol; len=$lcol+0
      if (g!="" && g!="no-hit" && g!="None" && g!="undef") sum[g]+=len
    }
    END {
      best=""; bestlen=0
      for (k in sum) if (sum[k] > bestlen) { bestlen=sum[k]; best=k }
      print best
    }' "$GENUS_TABLE")

  ## Comma-separated list of all classified genera EXCEPT the target
  ## (used to strip contaminants from clean-Curto isolates while keeping no-hit).
  OTHER_GENERA=$(awk -v tg="$TARGET_GENUS" '
    NR==1 { for (i=1;i<=NF;i++) if ($i ~ /genus/) gcol=i; next }
    { g=$gcol
      if (g!="" && g!="no-hit" && g!="None" && g!="undef" && g!=tg) seen[g]=1 }
    END { out=""; for (k in seen) out=(out=="")?k:out","k; print out }' "$GENUS_TABLE")

  echo "  dominant genus (length-weighted): '${DOMINANT_GENUS:-<none classified>}'"
  [ -n "$OTHER_GENERA" ] && echo "  other genera present: ${OTHER_GENERA}"
else
  echo "WARN: genus table not produced -- cannot determine dominant genus" >&2
fi

## ================= STAGE 4: filter -> FASTA =================
## Policy depends on whether the dominant genus is the target (Curtobacterium):
##
##  A) Curto-dominant ("clean") isolate:
##       keep Curtobacterium + no-hit, STRIP other (contaminant) genera,
##       then length + cov + gc.
##       -> preserves unclassified contigs that are likely real Curto.
##
##  B) Non-Curto-dominant ("mixed"/mislabeled) isolate:
##       keep ONLY the dominant genus (drop no-hit + everything else),
##       then length + cov (no gc).
##       -> isolates the actual organism out of a contaminated assembly.
##
## blobtools' --invert is global, so the taxonomic step and the threshold
## step run as TWO chained filters (threshold Mins must not be inverted).
COV_FIELD=reads_cov

if [ -z "$DOMINANT_GENUS" ]; then
  echo "WARN: no dominant genus for $ISOLATE_ID; length+cov only, no taxonomic filter" >&2
  blobtools filter \
    --param "length--Min=${MINLEN}" \
    --param "${COV_FIELD}--Min=${MINCOV}" \
    --fasta "$ASM" \
    --table "$TABLES/${ISOLATE_ID}.filtered.table.tsv" \
    --summary "$TABLES/${ISOLATE_ID}.summary.json" \
    "$DS"

elif [ "$DOMINANT_GENUS" = "$TARGET_GENUS" ]; then
  ## ---- A) CLEAN Curto isolate ----
  FILT_DS="$BLOBDIRS/${ISOLATE_ID}_taxfilt"
  rm -rf "$FILT_DS"
  if [ -n "$OTHER_GENERA" ]; then
    ## Exclude contaminant genera (default --Keys behaviour); keeps Curto + no-hit.
    echo "  step 1: CLEAN Curto -- keep Curtobacterium + no-hit, strip [${OTHER_GENERA}]"
    blobtools filter \
      --param "bestsumorder_genus--Keys=${OTHER_GENERA}" \
      --output "$FILT_DS" \
      "$DS"
  else
    ## Nothing to strip; carry the dataset forward unchanged.
    echo "  step 1: CLEAN Curto -- no other genera present, nothing to strip"
    cp -r "$DS" "$FILT_DS"
  fi
  echo "  step 2: length + ${COV_FIELD} + gc"
  blobtools filter \
    --param "length--Min=${MINLEN}" \
    --param "${COV_FIELD}--Min=${MINCOV}" \
    --param "gc--Min=${MINGC}" \
    --fasta "$ASM" \
    --table "$TABLES/${ISOLATE_ID}.filtered.table.tsv" \
    --summary "$TABLES/${ISOLATE_ID}.summary.json" \
    "$FILT_DS"

else
  ## ---- B) NON-Curto (mixed / mislabeled) isolate ----
  FILT_DS="$BLOBDIRS/${ISOLATE_ID}_keep_${DOMINANT_GENUS}"
  rm -rf "$FILT_DS"
  echo "  NOTE: $ISOLATE_ID dominant genus is '${DOMINANT_GENUS}', not ${TARGET_GENUS}"
  echo "  step 1: keep ONLY genus '${DOMINANT_GENUS}' (drop no-hit + others)"
  blobtools filter \
    --param "bestsumorder_genus--Keys=${DOMINANT_GENUS}" \
    --invert \
    --output "$FILT_DS" \
    "$DS"
  echo "  step 2: length + ${COV_FIELD}  (GC filter not applied for non-Curto)"
  blobtools filter \
    --param "length--Min=${MINLEN}" \
    --param "${COV_FIELD}--Min=${MINCOV}" \
    --fasta "$ASM" \
    --table "$TABLES/${ISOLATE_ID}.filtered.table.tsv" \
    --summary "$TABLES/${ISOLATE_ID}.summary.json" \
    "$FILT_DS"
fi

## Move the produced filtered FASTA to the tidy output dir
SRC_FILT="$ASSEMBLIES/${ISOLATE_ID}_scaffolds.filtered.fasta"
if [ -f "$SRC_FILT" ]; then
  mv "$SRC_FILT" "$FILTERED/${ISOLATE_ID}_filtered_scaffolds.fasta"
  echo "  filtered FASTA -> $FILTERED/${ISOLATE_ID}_filtered_scaffolds.fasta"
else
  echo "WARN: expected filtered FASTA not found: $SRC_FILT" >&2
fi

conda deactivate
echo "[$(date)] Done: $ISOLATE_ID"
