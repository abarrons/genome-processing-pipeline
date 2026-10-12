#!/bin/bash

#SBATCH --job-name=fastAAI  ## Name of the job.
#SBATCH -A JMARTINY_LAB                ## Account to charge; personal or lab
#SBATCH -p standard                    ## Partition/queue name
#SBATCH --time=24:00:00
#SBATCH --nodes=1                     ## (-N) number of nodes the job will use
#SBATCH --mem-per-cpu=4GB             ##
#SBATCH --cpus-per-task=4            ## number of CPUs to be used per task
#SBATCH -o fastAAI_%j.out             ## File to which STDOUT will be written, %j inserts jobid
#SBATCH -e fastAAI_%j.err             ## File to which STDERR will be written, %j inserts jobid


#we will use fastAAI to calculate AAI for all pairs of genomes

# Set folder paths
BASE=/dfs10/bio/abarrons/CA_curto/isolates
WD=$BASE/all_good_genomes
PROGRAM_PATH=/dfs10/bio/abarrons/programs
THREADS=4

#load required modules

module load python/3.10.2

# move to working directory
cd $WD

# prepare curto_genomes using symlinks so original files are never moved or deleted
rm -rf "$WD/curto_genomes"
mkdir -p "$WD/curto_genomes"

EXCLUDE_TMP=$(mktemp)
if [ -s "$WD/all_genomes_BUSCO_score.txt" ]; then
  grep "NO" "$WD/all_genomes_BUSCO_score.txt" | grep -v "28G4B03" | cut -d " " -f1 | awk '{print $0 ".fna"}' > "$EXCLUDE_TMP"
fi
cat << 'EOF' >> "$EXCLUDE_TMP"
25G3B06.fna
23G4B11.fna
S.43.fna
Curtobacterium_GCF_000710345.2_CURSPS6_genomic.fna
EOF

find . -maxdepth 1 -name '*.fna' ! -name 'Frigobacterium*' | sed 's|^\./||' | sort > all_files.txt
grep -v -F -f "$EXCLUDE_TMP" all_files.txt > files.txt
rm -f "$EXCLUDE_TMP" all_files.txt

while IFS= read -r file; do
  [ -f "$file" ] && ln -sf "$PWD/$file" "$WD/curto_genomes/$file"
done < files.txt
rm -f files.txt

# remove previous runs files
rm -rf fastaai_build
rm -rf FastAAI
rm -rf FastAAI_matrix

# activate FastAAI environment
source $PROGRAM_PATH/fastAAI/bin/activate

# now run fastAAI

fastaai build_db \
  --genomes curto_genomes/ \
  --threads $THREADS \
  --verbose \
  --output fastaai_build \
  --database fastaai_db.db \
  --compress

# clean up symlinks
rm -rf "$WD/curto_genomes"

fastaai db_query \
  --query fastaai_build/database/fastaai_db.db \
  --target fastaai_build/database/fastaai_db.db \
  --threads $THREADS \
  --verbose \
  --do_stdev

fastaai db_query \
  --query fastaai_build/database/fastaai_db.db \
  --target fastaai_build/database/fastaai_db.db \
  --threads $THREADS \
  --verbose \
  --output FastAAI_matrix \
  --output_style matrix
