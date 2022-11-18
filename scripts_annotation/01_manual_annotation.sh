#!/usr/bin/env bash

#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00
#SBATCH --job-name=manual_annotation
#SBATCH --mail-user=marco.visani@students.unibe.ch
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=pall
#SBATCH --array=0,1

module load UHTS/Analysis/SeqKit/0.13.2
module load Blast/ncbi-blast/2.9.0+
module load SequenceAnalysis/MultipleSequenceAlignment/clustalw2/2.1
module load Emboss/EMBOSS/6.6.0

ASSEMBLIES=('/data/users/mvisani/genome_assembly/polishing/canu/canu_pilon.fasta' '/data/users/mvisani/genome_assembly/polishing/flye/flye_pilon.fasta')
WORK=('canu' 'flye')
PROJECTDIR=/data/users/mvisani/annotation_genome

GENOME=${ASSEMBLIES[$SLURM_ARRAY_TASK_ID]}
OUTDIR=$PROJECTDIR/transpos_element_annotation/manual/
WORKDIR=$OUTDIR/${WORK[$SLURM_ARRAY_TASK_ID]}
if ! [ -d $OUTDIR ]; then mkdir $OUTDIR; fi
if ! [ -d $WORKDIR ]; then mkdir $WORKDIR; fi
rm -rf $WORKDIR/*

cd $WORKDIR

seqkit fx2tab -l -n $GENOME > $WORKDIR/id.txt
sort -n -r -k 2 $WORKDIR/id.txt > $WORKDIR/id_sorted.txt
awk 'NR==1{print $1}' $WORKDIR/id_sorted.txt > $WORKDIR/id_max.txt

seqkit grep \
    -n -f $WORKDIR/id_max.txt \
    $GENOME \
    -o $WORKDIR/contig.fasta

if ! [ -f $WORKDIR/contig.fasta ]; then exit 0; fi

seqkit sliding \
    -s 500 \
    -W 500 \
    -g $WORKDIR/contig.fasta \
    -o $WORKDIR/contig_windows.fasta

if ! [ -f $WORKDIR/contig_windows.fasta ]; then exit 0; fi

mkdir $WORKDIR/blast_db

makeblastdb \
    -in $WORKDIR/contig_windows.fasta  \
    -dbtype nucl \
    -out $WORKDIR/blast_db/contig_windows

if ! [ -d $WORKDIR/blast_db ]; then exit 0; fi

blastn \
    -query $WORKDIR/contig_windows.fasta \
    -db $WORKDIR/blast_db/contig_windows \
    -num_threads 10 \
    -outfmt 6 \
    -perc_identity 80 \
    -max_hsps 1 \
    > $WORKDIR/contig_windows.blastn


awk '{print $1}' $WORKDIR/contig_windows.blastn > $WORKDIR/blast_id.txt
sort $WORKDIR/blast_id.txt > $WORKDIR/blast_id_sorted.txt
uniq -c $WORKDIR/blast_id_sorted.txt > $WORKDIR/blast_id_cnt.txt
sort -n -r -k 1 $WORKDIR/blast_id_cnt.txt > $WORKDIR/blast_id_cnt_sorted.txt
awk '1<=NR && NR<=50{print $2}' $WORKDIR/blast_id_cnt_sorted.txt > $WORKDIR/blast_50max.txt


seqkit \
    grep -n -f $WORKDIR/blast_50max.txt \
    $WORKDIR/contig_windows.fasta\
    -o $WORKDIR/contig_windows_TOP50.fa


#Blast a window containing a candidate TE sequence against the genome
makeblastdb \
    -in $GENOME \
    -dbtype nucl \
    -out $WORKDIR/blast_db/${WORK[$SLURM_ARRAY_TASK_ID]}


#I've chose the fifth contif because there were many sequences on flye and a couple on canu (arbitrary)
WINDOW=('tig00000001_pilon_sliding:9501-10000' 'scaffold_601_pilon_sliding:3008001-3008500')
awk -v seq=${WINDOW[$SLURM_ARRAY_TASK_ID]} \
    -v RS='>' '$1 == seq {print RS $0}' $WORKDIR/contig_windows_TOP50.fa > $WORKDIR/window.fa



blastn \
    -query $WORKDIR/window.fa \
    -db $WORKDIR/blast_db/${WORK[$SLURM_ARRAY_TASK_ID]} \
    -outfmt 6 \
    -num_threads $SLURM_CPUS_PER_TASK \
    -perc_identity 80 \
    -qcov_hsp_perc 80 \
    -max_hsps 1 \
    > $WORKDIR/window.blastn

awk '$10 < $9 {print($2"\t"$10-1"\t"$9)}' $WORKDIR/window.blastn > $WORKDIR/window.bed
awk '$10 > $9 {print($2"\t"$9-1"\t"$10)}' $WORKDIR/window.blastn >> $WORKDIR/window.bed

seqkit subseq \
    --bed $WORKDIR/window.bed $GENOME \
    -u 2000 \
    -d 2000 \
    > $WORKDIR/window.bed.fa

clustalw2 -INFILE=$WORKDIR/window.bed.fa -OUTFILE=$WORKDIR/window.bed.aln

cons \
    -sequence $WORKDIR/window.bed.aln \
    -outseq $WORKDIR/TE5_family.cons \
    -name TE5_cons



makeblastdb \
    -in /data/courses/assembly-annotation-course/Brassicaceae_repbase_all_march2019.fasta \
    -dbtype nucl \
    -out $WORKDIR/blast_db/Brassicaceae

blastn \
    -query $WORKDIR/TE5_family.cons \
    -db $WORKDIR/blast_db/Brassicaceae \
    -outfmt 6 \
    -num_threads $SLURM_CPUS_PER_TASK \
    -max_hsps 1 \
    > $WORKDIR/TE5_consensus_blasted.blastn