#!/usr/bin/env bash

#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=4G
#SBATCH --time=05:00:00
#SBATCH --job-name=salmon
#SBATCH --mail-user=marco.visani@students.unibe.ch
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=pall
#SBATCH --array=0,1

module load UHTS/Analysis/busco/4.1.4;

WORK=('canu' 'flye')
PROJECTDIR=/data/users/mvisani/annotation_genome
OUTDIR=$PROJECTDIR/salmon
WORKDIR=$OUTDIR/${WORK[$SLURM_ARRAY_TASK_ID]}
BUSCO=$PROJECTDIR/busco/${WORK[$SLURM_ARRAY_TASK_ID]}
SALMONDIR=/data/users/mvisani/annotation_genome/salmon-1.9.0_linux_x86_64/bin

MAKEOUT=$PROJECTDIR/coding_region_annotation/${WORK[$SLURM_ARRAY_TASK_ID]}
SALMONINDEX=$WORKDIR/salmon_index
READS=/data/users/mvisani/genome_assembly/participant_2/Illumina

mkdir $PROJECTDIR/busco
if ! [ -d $OUTDIR ]; then mkdir $OUTDIR; fi 
if ! [ -d $WORKDIR ]; then mkdir $WORKDIR; fi 
if ! [ -d $BUSCO ]; then mkdir $BUSCO; fi 

cd $BUSCO
rm -rf $BUSCO/*
rm -rf $OUTDIR/*
busco \
    -i $MAKEOUT/${WORK[$SLURM_ARRAY_TASK_ID]}_pilon.all.maker.proteins.fasta \
    -l brassicales_odb10 \
    -m proteins \
    -c 4 \
    --out busco_out

$SALMONDIR/salmon index \
    -t $MAKEOUT/${WORK[$SLURM_ARRAY_TASK_ID]}_pilon.all.maker.transcripts.fasta \
    -i $SALMONINDEX \
    -k 31

$SALMONDIR/salmon quant \
    -i $SALMONINDEX \
    -l A \
    -1 $READS/ERR3624579_1.fastq.gz \
    -2 $READS/ERR3624579_2.fastq.gz \
    --validateMappings \
    -o $OUTDIR/${WORK[$SLURM_ARRAY_TASK_ID]}_transcripts_quant