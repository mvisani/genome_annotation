#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00
#SBATCH --job-name=TE_phylogeny
#SBATCH --mail-user=marco.visani@students.unibe.ch
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=pall
#SBATCH --array=0,1

module load UHTS/Analysis/BEDTools/2.29.2;

#set variables 
WORK=('canu' 'flye')
genome=('/data/users/mvisani/annotation_genome/transpos_element_annotation/auto/canu/canu_pilon.fasta.mod.EDTA.intact.gff3' \
    '/data/users/mvisani/annotation_genome/transpos_element_annotation/auto/flye/flye_pilon.fasta.mod.EDTA.intact.gff3')

ASSEMBLIES=('/data/users/mvisani/genome_assembly/polishing/canu/canu_pilon.fasta' \
    '/data/users/mvisani/genome_assembly/polishing/flye/flye_pilon.fasta')

PROJECTDIR=/data/users/mvisani/annotation_genome

OUTDIR=$PROJECTDIR/transpos_element_annotation/dating/
WORKDIR=$OUTDIR/${WORK[$SLURM_ARRAY_TASK_ID]}

if ! [ -d $OUTDIR ]; then mkdir $OUTDIR; fi
if ! [ -d $WORKDIR ]; then mkdir $WORKDIR; fi

cd $WORKDIR
rm -rf $WORKDIR/*

awk '$3~/retrotransposon/' ${genome[$SLURM_ARRAY_TASK_ID]} \
    > $WORKDIR/${WORK[$SLURM_ARRAY_TASK_ID]}.mod.EDTA.intact.gff3_LTR

sed -i 's/ID.\+Name=//' $WORKDIR/${WORK[$SLURM_ARRAY_TASK_ID]}.mod.EDTA.intact.gff3_LTR
sed -i 's/;.\+//' $WORKDIR/${WORK[$SLURM_ARRAY_TASK_ID]}.mod.EDTA.intact.gff3_LTR

#rename pil to pilon because bedtools is shit
sed -i 's/_pi\s/_pilon\t/g' $WORKDIR/${WORK[$SLURM_ARRAY_TASK_ID]}.mod.EDTA.intact.gff3_LTR
sed -i 's/_pil\s/_pilon\t/g' $WORKDIR/${WORK[$SLURM_ARRAY_TASK_ID]}.mod.EDTA.intact.gff3_LTR
sed -i 's/_pilo\s/_pilon\t/g' $WORKDIR/${WORK[$SLURM_ARRAY_TASK_ID]}.mod.EDTA.intact.gff3_LTR

awk '{print($1,$2,$9,$4,$5,$6,$7,$8,$3)}' $WORKDIR/${WORK[$SLURM_ARRAY_TASK_ID]}.mod.EDTA.intact.gff3_LTR |sed 's/ /\t/g' \
    > $WORKDIR/${WORK[$SLURM_ARRAY_TASK_ID]}.mod.EDTA.intact.gff3_LTR_ref

bedtools getfasta \
    -fi ${ASSEMBLIES[$SLURM_ARRAY_TASK_ID]} \
    -bed $WORKDIR/${WORK[$SLURM_ARRAY_TASK_ID]}.mod.EDTA.intact.gff3_LTR_ref \
    -name \
    > $WORKDIR/${WORK[$SLURM_ARRAY_TASK_ID]}.mod.EDTA.intact.LTR.fa

sed -i 's/:/_/' $WORKDIR/${WORK[$SLURM_ARRAY_TASK_ID]}.mod.EDTA.intact.LTR.fa
sed -i 's/>/>Ler_/' $WORKDIR/${WORK[$SLURM_ARRAY_TASK_ID]}.mod.EDTA.intact.LTR.fa
sed -i 's/:/_/g' $WORKDIR/${WORK[$SLURM_ARRAY_TASK_ID]}.mod.EDTA.intact.LTR.fa

#perl the stuff
PERL=$WORKDIR/perl
if ! [ -d $PERL ]; then mkdir $PERL; fi
cp /data/courses/assembly-annotation-course/CDS_annotation/scripts/LTR $PERL
cp /data/courses/assembly-annotation-course/CDS_annotation/scripts/split_flat $PERL
cp /data/courses/assembly-annotation-course/CDS_annotation/scripts/date_pair $PERL

cd $PERL

module add Emboss/EMBOSS/6.6.0;
perl split_flat $WORKDIR/${WORK[$SLURM_ARRAY_TASK_ID]}.mod.EDTA.intact.LTR.fa
perl LTR Ler_ N
perl date_pair > $WORKDIR/date_pair.txt