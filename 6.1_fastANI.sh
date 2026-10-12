#!/bin/bash

#SBATCH --job-name=fastANI  ## Name of the job.
#SBATCH -A JMARTINY_LAB                ## Account to charge; personal or lab
#SBATCH -p standard                    ## Partition/queue name
#SBATCH --time=24:00:00
#SBATCH --nodes=1                     ## (-N) number of nodes the job will use
#SBATCH --mem-per-cpu=4GB             ##
#SBATCH --cpus-per-task=32            ## number of CPUs to be used per task
#SBATCH -o fastANI_%j.out             ## File to which STDOUT will be written, %j inserts jobid
#SBATCH -e fastANI_%j.err             ## File to which STDERR will be written, %j inserts jobid


#we will use fastANI to calculate ANI for all pairs of genomes

# Set folder paths
BASE=/dfs10/bio/abarrons/CA_curto/isolates
WD=$BASE/all_good_genomes
PROGRAM_PATH=/dfs10/bio/abarrons/programs

# move to working directory
cd $WD

# generate list of genomes to analyze, excluding Frigobacterium and any isolates in isolates_to_remove.txt
EXCLUDE_TMP=$(mktemp)
if [ -s "$WD/isolates_to_remove.txt" ]; then
  cut -f2 "$WD/isolates_to_remove.txt" > "$EXCLUDE_TMP"
fi

find . -maxdepth 1 -name '*.fna' ! -name 'Frigobacterium*' | sed 's|^\./||' | sort > all_files.txt

if [ -s "$EXCLUDE_TMP" ]; then
  grep -v -F -f "$EXCLUDE_TMP" all_files.txt > files.txt
else
  cp all_files.txt files.txt
fi
rm -f "$EXCLUDE_TMP" all_files.txt

# make output directory
mkdir -p "$WD/fastANI_out_reduced"

# now run fastANI
$PROGRAM_PATH/fastANI --ql files.txt --rl files.txt -t 32 --matrix -o "$WD/fastANI_out_reduced/fastANI.out"
