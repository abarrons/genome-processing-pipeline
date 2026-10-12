#!/bin/bash

#SBATCH --job-name=checkm2          ## Name of the job.
#SBATCH -A JMARTINY_LAB             ## Account to charge
#SBATCH -p standard                 ## Partition/queue
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16          ## CheckM2/DIAMOND scales with threads
#SBATCH --mem=64G                   ## DIAMOND against the CheckM2 DB wants memory
#SBATCH -o checkm2_%j.out
#SBATCH -e checkm2_%j.err

## Runs CheckM2 on a folder of assemblies in one pass, producing a single
## quality_report.tsv with completeness + contamination for every genome.
##
## Usage:
##   sbatch 4.0_checkM2.sh                 # run on the consolidated set (default)
##   sbatch 4.0_checkM2.sh <input_dir>     # run on any folder (e.g. the duplicates)
##
## The output dir is derived from the input dir name, so a second run on the
## duplicates folder writes to its own output and won't clobber the main report.
##
## Downstream transition note:
##   Genomes meeting quality criteria in quality_report.tsv (e.g., completeness >= 90%,
##   contamination < 5%) are curated and placed into $BASE/all_good_genomes/ as .fna files,
##   along with any reference/outgroup genomes used for phylogenetic rooting or comparative
##   genomics (stages 5 through 7).

## ------------------------------------------------------------
## Initialize conda + activate the checkm2 env
module load miniconda3/25.11.1
. ~/.mycondainit-25.11.1
conda activate checkm2

## CheckM2 DIAMOND database
export CHECKM2DB=/dfs10/bio/abarrons/checkM2DB/CheckM2_database/uniref100.KO.1.dmnd

THREAD=16

## ------------------------------------------------------------
## Input folder: first arg, or the consolidated set by default.
INDIR="${1:-/dfs10/bio/abarrons/CA_curto/NCBI_submission/assemblies}"

## Derive a distinct output dir from the input folder's name, so main-set and
## duplicates runs land in separate places under NCBI_submission/.
NBASE=/dfs10/bio/abarrons/CA_curto/NCBI_submission
INNAME=$(basename "$INDIR")
OUTDIR="$NBASE/checkm2_out_${INNAME}"

[ -d "$INDIR" ] || { echo "ERROR: input dir not found: $INDIR" >&2; exit 1; }

## Sanity: the DB file must exist (a folder or missing file is the usual failure)
if [ ! -f "$CHECKM2DB" ]; then
  echo "ERROR: CheckM2 DB file not found: $CHECKM2DB" >&2
  echo "       (CHECKM2DB must be the .dmnd FILE, not a directory)" >&2
  exit 1
fi

## ---- Clean any prior output for this input so the run starts fresh ----
if [ -d "$OUTDIR" ]; then
  echo "[CLEAN] removing prior output: $OUTDIR"
  rm -rf "$OUTDIR"
fi
mkdir -p "$OUTDIR"


## Count genomes present (handles both *_filtered_scaffolds.fasta and plain .fasta)
NGEN=$(ls "$INDIR"/*.fasta 2>/dev/null | wc -l)
echo "[$(date)] CheckM2 on $NGEN genomes in $INDIR"
echo "  DB:     $CHECKM2DB"
echo "  Output: $OUTDIR"

checkm2 predict \
  --threads $THREAD \
  --input "$INDIR" \
  --extension .fasta \
  --output-directory "$OUTDIR" \
  --force

echo "[$(date)] Done. Report: $OUTDIR/quality_report.tsv"
