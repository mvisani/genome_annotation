#!/usr/bin/env bash

#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=dup_gen
#SBATCH --mail-user=marco.visani@students.unibe.ch
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=pall
#SBATCH --array=0,1

WORK=('canu' 'flye')
PROJECTDIR=/data/users/mvisani/annotation_genome
GFF=$PROJECTDIR/coding_region_annotation/${WORK[$SLURM_ARRAY_TASK_ID]}/${WORK[$SLURM_ARRAY_TASK_ID]}_pilon.all.gff
NNU=/data/courses/assembly-annotation-course/CDS_annotation/NNU_mRNA_single_model.gff
nnu_pep=/data/courses/assembly-annotation-course/CDS_annotation/NNU.pep.fa.ref.single_model
OUTDIR=$PROJECTDIR/out_08_dup_gen
WORKDIR=$OUTDIR/${WORK[$SLURM_ARRAY_TASK_ID]}
dup_gen=/data/users/mvisani/annotation_genome/DupGen_finder


if ! [ -d $OUTDIR ]; then mkdir $OUTDIR; fi 
if ! [ -d $WORKDIR ]; then mkdir $WORKDIR; fi 

cd $WORKDIR
rm -rf $WORKDIR/*

awk '{sub(/ID=/, "", $9); sub(/;.*/, "", $9); sub(/:.*/, "", $9); printf ("%s\t%s\t%s\t%s\n", $1, $9, $4, $5)}' $GFF | \
    grep "\-RA" \
    > $WORKDIR/Ath.gff

cat $WORKDIR/Ath.gff $NNU > $WORKDIR/Ath_Nnu.gff 

module load Blast/ncbi-blast/2.9.0+
makeblastdb \
    -in $PROJECTDIR/coding_region_annotation/${WORK[$SLURM_ARRAY_TASK_ID]}/${WORK[$SLURM_ARRAY_TASK_ID]}_pilon.all.maker.proteins.fasta \
    -dbtype prot \
    -out $WORKDIR/blast_db/Ath

blastp \
    -query $PROJECTDIR/coding_region_annotation/${WORK[$SLURM_ARRAY_TASK_ID]}/${WORK[$SLURM_ARRAY_TASK_ID]}_pilon.all.maker.proteins.fasta \
    -db $WORKDIR/blast_db/Ath \
    -num_threads $SLURM_CPUS_PER_TASK \
    -outfmt 6 \
    -evalue 1e-10 \
    -max_target_seqs 5 \
    > $WORKDIR/Ath.blast

blastp \
    -query $nnu_pep \
    -db $WORKDIR/blast_db/Ath \
    -num_threads $SLURM_CPUS_PER_TASK \
    -outfmt 6 \
    -evalue 1e-10 \
    -max_target_seqs 5 \
    > $WORKDIR/Ath_Nnu.blast

mkdir dup_gen
cp $WORKDIR/Ath_Nnu.blast ./dup_gen/
cp $WORKDIR/Ath.blast ./dup_gen/
cp $WORKDIR/Ath_Nnu.gff ./dup_gen/
cp $WORKDIR/Ath.gff ./dup_gen/

$dup_gen/DupGen_finder.pl -i $WORKDIR/dup_gen -t Ath -c Nnu -o results