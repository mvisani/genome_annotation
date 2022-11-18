#!/usr/bin/env bash

#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=star
#SBATCH --mail-user=marco.visani@students.unibe.ch
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=pall
#SBATCH --array=0,1

module load UHTS/Assembler/cufflinks/2.2.1;

WORK=('canu' 'flye')
PROJECTDIR=/data/users/mvisani/annotation_genome
GFF=$PROJECTDIR/coding_region_annotation/${WORK[$SLURM_ARRAY_TASK_ID]}/maker_functional.gff
OUTDIR=$PROJECTDIR/star
WORKDIR=$OUTDIR/${WORK[$SLURM_ARRAY_TASK_ID]}

if ! [ -d $OUTDIR ]; then mkdir $OUTDIR; fi 
if ! [ -d $WORKDIR ]; then mkdir $WORKDIR; fi 

cd $WORKDIR
rm -rf $WORKDIR/*

gffread -E $GFF -T -o  ${WORK[$SLURM_ARRAY_TASK_ID]}_functional.gtf


READS=/data/users/mvisani/genome_assembly/participant_2/RNAseq
ASSEMBLIES=('/data/users/mvisani/genome_assembly/polishing/canu/canu_pilon.fasta' \
    '/data/users/mvisani/genome_assembly/polishing/flye/flye_pilon.fasta')

module add UHTS/Aligner/STAR/2.7.9a


zcat $READS/ERR754061_1.fastq.gz > read1.fasta
zcat $READS/ERR754061_2.fastq.gz > read2.fasta


STAR \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --runMode genomeGenerate \
    --genomeDir $WORKDIR \
    --genomeFastaFiles ${ASSEMBLIES[$SLURM_ARRAY_TASK_ID]} \
    --sjdbGTFfile ${WORK[$SLURM_ARRAY_TASK_ID]}_functional.gtf \
    --genomeSAindexNbases 12


STAR \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --genomeDir $WORKDIR \
    --readFilesIn read1.fasta read2.fasta \
    --outFileNamePrefix star \
    --outSAMtype BAM SortedByCoordinate \
    --outFilterMultimapNmax 10 \
    --outFilterMismatchNoverLmax 0.01 \
    --alignIntronMax 60000
