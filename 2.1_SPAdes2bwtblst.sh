#!/bin/bash

#SBATCH --job-name=SPAdes           ## Name of the job.
#SBATCH -A JMARTINY_LAB             ## Account to charge; personal or lab
#SBATCH -p standard                 ## Partition/queue name
#SBATCH --nodes=1                   ## (-N) number of nodes the job will use
#SBATCH --cpus-per-task=8           ## number of CPUs to be used per task
#SBATCH -o SPAdes_%A_%a.out         ## STDOUT; %A=jobid, %a=array index
#SBATCH -e SPAdes_%A_%a.err         ## STDERR; %A=jobid, %a=array index
## NOTE: --array is set at submission time by the launcher (submit_batch.sh),
##       so it is intentionally NOT hardcoded here.

## This script runs SPAdes over the bbduk quality-/polyG-filtered reads,
## maps raw reads back to the ASSEMBLY, and blasts contigs for taxa ID
## to prepare data for blobtools.

## ============================================================
## BATCH and CLEAN are passed in by submit_batch.sh via --export.
## Do NOT hardcode them here -- submit through the launcher.
##   BATCH  : 1stbatch | 2ndbatch | 3rdbatch
##   CLEAN  : 1 = wipe prior outputs before running (raw reads untouched)
## ------------------------------------------------------------
: "${BATCH:?BATCH not set -- submit via ./submit_batch.sh, not sbatch directly}"
CLEAN="${CLEAN:-0}"


BASE=/dfs10/bio/abarrons/CA_curto/isolates/${BATCH}
BBMAPDIR=/dfs10/bio/abarrons/programs/bbmap/resources
DB=/dfs10/bio/abarrons/databases
BIN=/dfs10/bio/abarrons/programs/blobology-master
raw_seqs=$BASE/raw_seqs/seqs

THREAD=8

## Load programs
module load SPAdes/3.15.4
module load bbmap/38.87
module load bowtie2/2.4.5
module load samtools/1.15.1
module load ncbi-blast/2.13.0

## ------------------------------------------------------------
## Optional FULL-BATCH cleanup. Off unless CLEAN=1.
## Runs ONLY on array task 1 so the tasks don't race to delete
## (and so it can't wipe another task's fresh output mid-run).
## Deletes ONLY derived output dirs -- raw_seqs is never touched.
## ------------------------------------------------------------
if [ "$CLEAN" = "1" ] && [ "${SLURM_ARRAY_TASK_ID:-1}" = "1" ]; then
  echo "[CLEAN] CLEAN=1 on task 1 -- wiping prior outputs under $BASE"
  rm -rf "$BASE/spades_out" \
         "$BASE/bowtie_dbs" \
         "$BASE/blast_out" \
         "$BASE/bowtie2_output"
  echo "[CLEAN] Done. raw_seqs left intact: $raw_seqs"
fi


## Output dirs (mkdir -p is safe to re-run / safe across array tasks)
mkdir -p $BASE/spades_out
mkdir -p $BASE/bowtie_dbs
mkdir -p $BASE/blast_out
mkdir -p $BASE/bowtie2_output
mkdir -p $BASE/spades_out/assemblies

cd $raw_seqs

## Per-task sample selection
FFILE=`ls *READ1-Sequences.txt.gz | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
RFILE=${FFILE//READ1-Sequences.txt.gz/READ2-Sequences.txt.gz}
ISOLATE_ID=$(echo ${FFILE} | cut -d "-" -f 4)

echo "[$(date)] Batch=$BATCH  Task=$SLURM_ARRAY_TASK_ID  Isolate=$ISOLATE_ID"
echo "  R1=$FFILE"
echo "  R2=$RFILE"

## ---- Read cleaning: adapters + quality + poly-G (two-color chemistry) ----
bbduk.sh in1=$FFILE in2=$RFILE \
  out1=${ISOLATE_ID}.R1.clean.fq.gz out2=${ISOLATE_ID}.R2.clean.fq.gz \
  ref=$BBMAPDIR/adapters.fa ktrim=r k=25 mink=11 hdist=1 tbo tpe \
  qtrim=rl trimq=10 trimpolyg=10 minlen=50 ordered=t

## ---- Assemble with SPAdes ----
## Wipe this isolate's own output dir first. SPAdes refuses to start clean in a
## non-empty dir -- it tries to RESUME and crashes in generate_configs when the
## prior attempt's configs are incomplete. Clearing it per task guarantees a
## fresh start regardless of any leftover/half-built dir. (Derived output only;
## raw reads are untouched.)
rm -rf $BASE/spades_out/${ISOLATE_ID}
spades.py \
  -1 ${ISOLATE_ID}.R1.clean.fq.gz -2 ${ISOLATE_ID}.R2.clean.fq.gz \
  -k 31,41,51,61,71,81,91,95,101,105,111 --careful -t $THREAD --cov-cutoff auto \
  -o $BASE/spades_out/${ISOLATE_ID}

rm -f ${ISOLATE_ID}.R1.clean.fq.gz
rm -f ${ISOLATE_ID}.R2.clean.fq.gz

## ---- Map raw reads back to the assembly ----
cd $BASE/spades_out

bowtie2-build \
  ${ISOLATE_ID}/scaffolds.fasta \
  $BASE/bowtie_dbs/${ISOLATE_ID} \
  --threads $THREAD --verbose

bowtie2 \
  -x $BASE/bowtie_dbs/${ISOLATE_ID} \
  -U <($BIN/shuffleSequences_fastx.pl 4 <(zcat $raw_seqs/*${ISOLATE_ID}*READ1*.txt.gz) <(zcat $raw_seqs/*${ISOLATE_ID}*READ2*.txt.gz)) \
  --very-fast-local -k 1 -t --reorder --mm \
  --threads $THREAD |\
  samtools view -bh -@ $THREAD | samtools sort -@ $THREAD > $BASE/bowtie2_output/${ISOLATE_ID}.sorted.bam

samtools index -@ $THREAD $BASE/bowtie2_output/${ISOLATE_ID}.sorted.bam

## ---- BLAST contigs for taxa ID ----
blastn \
  -task megablast \
  -query $BASE/spades_out/${ISOLATE_ID}/scaffolds.fasta \
  -db $DB/nt \
  -outfmt '6 qseqid staxids bitscore std' \
  -max_target_seqs 1 \
  -max_hsps 1 \
  -evalue 1e-5 \
  -num_threads $THREAD \
  -out $BASE/blast_out/${ISOLATE_ID}.blast.out

## ---- Tidy: move assembly to a collected folder ----
mv $BASE/spades_out/${ISOLATE_ID}/scaffolds.fasta \
   $BASE/spades_out/assemblies/${ISOLATE_ID}_scaffolds.fasta

echo "[$(date)] Done: $ISOLATE_ID"
