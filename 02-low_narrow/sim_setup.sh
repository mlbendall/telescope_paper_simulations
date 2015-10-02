#! /bin/bash

### Prerequisites: R and python
module load python2.7/2.7.6
module load R/3.1.1

### HERV annotations
HERV_FAS=$(dirname $PWD)/annotations/knownHERV.hg19.fna
HERV_GTF=$(dirname $PWD)/annotations/knownHERV.hg19.gtf

### Identify candidate loci for simulation. 
# Candidates are between 7000 and 10000 bps; HML2 and HML5 are included in HERVK and 
# HERVFRD is combined with HERVW

python get_candidates.py < $HERV_GTF > candidates.txt

### Simulate expression
# Randomly draw 1 locus from each family for each simulation. Each locus is equally expressed.

Rscript sim_expression.R 123456

### Simulate sequencing

# We are using a development version of the Polyester package by Alyssa C. Frazee
# Preprint is published in bioRxiv:
# Frazee, A.C. et al., 2014. Polyester: simulating RNA-seq datasets with differential 
# transcript expression, Cold Spring Harbor Labs Journals.

# Set POLYESTER_DIR environment variable
export POLYESTER_DIR=$(dirname $PWD)/polyester.git

Rscript sim_reads.R $HERV_FAS simulation_proportions.txt 40 123456    # 10 per sample
Rscript sim_reads.R $HERV_FAS simulation_proportions.txt 80 123456    # 20
Rscript sim_reads.R $HERV_FAS simulation_proportions.txt 160 123456   # 40
Rscript sim_reads.R $HERV_FAS simulation_proportions.txt 320 123456   # 80
Rscript sim_reads.R $HERV_FAS simulation_proportions.txt 640 123456   # 160
Rscript sim_reads.R $HERV_FAS simulation_proportions.txt 1280 123456  # 320
Rscript sim_reads.R $HERV_FAS simulation_proportions.txt 2560 123456  # 640

for sd in simdata*; do
    echo $sd
    pigz $sd/*.fasta
    tail -n+2 $sd/samples.txt | cut -f2 | while read samp; do
        mkdir $sd/$samp
        mv $sd/${samp}_1.fasta.gz $sd/$samp
        mv $sd/${samp}_2.fasta.gz $sd/$samp
    done
    echo -e '*\n!.gitignore' > $sd/.gitignore
done
