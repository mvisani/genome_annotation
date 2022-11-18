#!/usr/bin/env bash

#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=PAML
#SBATCH --mail-user=marco.visani@students.unibe.ch
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=pall
#SBATCH --array=0,1

module load UHTS/Analysis/SeqKit/0.13.2

WORK=('canu' 'flye')
PROJECTDIR=/data/users/mvisani/annotation_genome
OUTDIR=$PROJECTDIR/out_09
WORKDIR=$OUTDIR/${WORK[$SLURM_ARRAY_TASK_ID]}
ALGO=/data/courses/assembly-annotation-course/CDS_annotation/scripts


if ! [ -d $OUTDIR ]; then mkdir $OUTDIR; fi 
if ! [ -d $WORKDIR ]; then mkdir $WORKDIR; fi 

cd $WORKDIR
rm -rf $WORKDIR/*

cut -f 1 $PROJECTDIR/out_08_dup_gen/${WORK[$SLURM_ARRAY_TASK_ID]}/results/Ath.wgd.genes > $WORKDIR/ID.txt

seqkit grep \
    -f $WORKDIR/ID.txt \
    $PROJECTDIR/coding_region_annotation/${WORK[$SLURM_ARRAY_TASK_ID]}/${WORK[$SLURM_ARRAY_TASK_ID]}_pilon.all.maker.transcripts.fasta \
    -o $WORKDIR/Ath.wgd.genes.fa

seqkit translate Ath.wgd.genes.fa -o Ath.wgd.genes.fa.proteins

sed -i 's/_frame=1/_p/' Ath.wgd.genes.fa.proteins

cut -f 1,3 $PROJECTDIR/out_08_dup_gen/${WORK[$SLURM_ARRAY_TASK_ID]}/results/Ath.wgd.pairs > $WORKDIR/wgd_pairs.csv
sed -i 's/RA/RA_p/g' $WORKDIR/wgd_pairs.csv

mkdir result
cp $WORKDIR/wgd_pairs.csv ./result/
cp $WORKDIR/Ath.wgd.genes.fa ./result/
cp $WORKDIR/Ath.wgd.genes.fa.proteins ./result/

cd result
$ALGO/split_flat Ath.wgd.genes.fa
$ALGO/split_flat Ath.wgd.genes.fa.proteins

# 5. Align proteins with bestflash_from_list. Output: protein alignments (.pair)
module add Emboss/EMBOSS/6.6.0
$ALGO/bestflash_from_list wgd_pairs.csv 

# 6. Generate codon-by-codon CDS alignments based on the protein alignments with pair_to_CDS_paml_pair.
#Output: CDS alignment in PAML format (.CDS_paml_pair)

$ALGO/pair_to_CDS_paml_pair LER pair

#7. Use yn00 from PAML (PAML_yn00_from_CDS_pair script) on all WGD pair alignments to calculate Ka/Ks (omega) values
module load Phylogeny/paml/4.9j
$ALGO/PAML_yn00_from_CDS_pair LER > $WORKDIR/PAML_yn00_results

cd $WORKDIR
awk '{print($1,$1,$6,$7,$5)}' PAML_yn00_results \
    |sed 's/ /\t/g' \
    |sed 's/__x__/\t/' \
    |sed 's/_p//g' \
    |cut -f 1,2,4,5,6 \
    |sed 's/dN=//' \
    |sed 's/dS=//' \
    |sed 's/omega=//' \
    |awk '$4<5' \
    > $WORKDIR/Ler.wgd.kaks

cp $PROJECTDIR/out_08_dup_gen/${WORK[$SLURM_ARRAY_TASK_ID]}/results/Ath.collinearity $WORKDIR/Ler.collinearity
perl $ALGO/add_ka_ks_to_collinearity_file.pl Ler

perl $ALGO/compute_ks_for_synteny_blocks.pl Ler.collinearity.kaks

python $ALGO/plot_syntenic_blocks_ks_distri.py Ler.synteny.blocks.ks.info 2 Ler