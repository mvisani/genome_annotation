#!/usr/bin/env bash

#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=EDTA_auto
#SBATCH --mail-user=marco.visani@students.unibe.ch
#SBATCH --mail-type=fail,end
#SBATCH --output=/data/users/mvisani/annotation_genome/scripts/TE_annotation_auto_%j.o
#SBATCH --error=/data/users/mvisani/annotation_genome/scripts/TE_annotation_auto_%j.e
#SBATCH --partition=pall
#SBATCH --array=0,1

ASSEMBLIES=('/data/users/mvisani/genome_assembly/polishing/canu/canu_pilon.fasta' \
  '/data/users/mvisani/genome_assembly/polishing/flye/flye_pilon.fasta')

CONTAINER=/data/courses/assembly-annotation-course/containers2
WORK=('canu' 'flye')
WORKDIR=/data/users/mvisani


OUTPUT=/data/users/mvisani/annotation_genome/transpos_element_annotation/auto
OUTDIR=$OUTPUT/${WORK[$SLURM_ARRAY_TASK_ID]}
if ! [ -d $OUTPUT ]; then mkdir $OUTPUT; fi
if ! [ -d $OUTDIR ]; then mkdir $OUTDIR; fi
rm -rf $OUTDIR/*

cd $OUTDIR

wget https://www.arabidopsis.org/download_files/Sequences/TAIR10_blastsets/TAIR10_cds_20110103_representative_gene_model
CDS=$OUTDIR/TAIR10_cds_20110103_representative_gene_model

singularity exec \
  --bind $CONTAINER \
  --bind $WORKDIR \
$CONTAINER/EDTA_v1.9.6.sif \
EDTA.pl \
  --genome ${ASSEMBLIES[${SLURM_ARRAY_TASK_ID}]} \
  --species others \
  --step all \
  --cds $CDS \
  --anno 1  \
  --threads $SLURM_CPUS_PER_TASK

