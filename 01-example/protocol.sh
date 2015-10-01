#! /bin/bash

# Prerequisites: R and python
module load python2.7/2.7.6
module load R/3.1.1

# We are using a development version of the Polyester package by Alyssa C. Frazee
# Preprint is published in bioRxiv:
# Frazee, A.C. et al., 2014. Polyester: simulating RNA-seq datasets with differential 
# transcript expression, Cold Spring Harbor Labs Journals.

# Set POLYESTER_DIR environment variable
export POLYESTER_DIR=$(dirname $PWD)/polyester.git

# HERV annotations
HERV_FAS=$(dirname $PWD)/annotations/knownHERV.hg19.fna
HERV_GTF=$(dirname $PWD)/annotations/knownHERV.hg19.gtf

# HERV_FAS=/lustre/groups/cbi/Repeatsome/data/annotations/knownHERV.hg19/knownHERV.hg19.fna
# HERV_GTF=/lustre/groups/cbi/Repeatsome/data/annotations/knownHERV.hg19/knownHERV.hg19.gtf

mkdir -p simdata
echo -e '*\n!.gitignore' > simdata/.gitignore

Rscript simulate.R $HERV_FAS simulation_proportions.txt 400 123456
cd simdata

# One folder per sample
samp=sample_01
mkdir $samp
mv ${samp}_1.fasta ${samp}
mv ${samp}_2.fasta ${samp}

### Run RepEnrich ########################################################################
module load repeatsome

# Index for bowtie
bt1idx=/lustre/groups/cbi/shared/References/Homo_sapiens/UCSC/hg19/Sequence/BowtieIndex/genome

# RepEnrich references
rmskall=/lustre/groups/cbi/Repeatsome/data/RepEnrich/rmsk_ALL.hg19

fasta1=${samp}/${samp}_1.fasta
fasta2=${samp}/${samp}_2.fasta

rp_repenrich_PE $bt1idx $(ls $rmskall/*.repenrich.txt) $rmskall/index $fasta1 $fasta2 > ${samp}/RepEnrich.log

### Run Telescope ########################################################################
module load repeatsome

# Align with Bowtie 2
bt2idx=/lustre/groups/cbi/shared/References/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome

fasta1=${samp}/${samp}_1.fasta
fasta2=${samp}/${samp}_2.fasta

rp_align_PE $bt2idx $fasta1 $fasta2 ${samp}/align_bt2 &> ${samp}/align_PE.log

# Reassign using Telescope
tsopts="--verbose --out_matrix --updated_sam --thetaPrior 200000 --maxIter 1000"
rp_telescope ${samp}/align_bt2.sam $HERV_GTF $samp/tsout knownHERV "$tsopts" &> ${samp}/tsout.log

# Tag alignments
module load telescope
module load samtools

telescope tag --gtffile $HERV_GTF < ${samp}/align_bt2.sam > ${samp}/align_bt2.tagged.sam
samtools view -uS ${samp}/align_bt2.tagged.sam | samtools sort - ${samp}/align_bt2.tagged
samtools index ${samp}/align_bt2.tagged.bam

# Generate pi_0, pi, x_init, and x_hat
telescope load --outparam pi_0 ${samp}/tsout/knownHERV-matrix_final.p > ${samp}/tsout/knownHERV-pi_0.txt
telescope load --outparam pi ${samp}/tsout/knownHERV-matrix_final.p > ${samp}/tsout/knownHERV-pi.txt
telescope load --outparam x_init ${samp}/tsout/knownHERV-matrix_final.p > ${samp}/tsout/knownHERV-x_init.txt
telescope load --outparam x_hat ${samp}/tsout/knownHERV-matrix_final.p > ${samp}/tsout/knownHERV-x_hat.txt

### Create the summary table #############################################################
python ../summarize_results.py > ../summary_tables.txt

cd ..
