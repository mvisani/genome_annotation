#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00
#SBATCH --job-name=TE_phylogeny
#SBATCH --mail-user=marco.visani@students.unibe.ch
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=pall
#SBATCH --array=0,1

module load UHTS/Analysis/SeqKit/0.13.2
module load SequenceAnalysis/MultipleSequenceAlignment/clustal-omega/1.2.4
module load Phylogeny/FastTree/2.1.10

#set variables 
WORK=('canu' 'flye')
genome=('canu_pilon.fasta' 'flye_pilon.fasta')
ASSEMBLIES=('/data/users/mvisani/genome_assembly/polishing/canu/canu_pilon.fasta' \
    '/data/users/mvisani/genome_assembly/polishing/flye/flye_pilon.fasta')

PROJECTDIR=/data/users/mvisani/annotation_genome


OUTDIR=$PROJECTDIR/transpos_element_annotation/phylogeny/auto/
WORKDIR=$OUTDIR/${WORK[$SLURM_ARRAY_TASK_ID]}

if ! [ -d $OUTDIR ]; then mkdir $OUTDIR; fi
if ! [ -d $WORKDIR ]; then mkdir $WORKDIR; fi

cd $WORKDIR
rm $WORKDIR/*


more $PROJECTDIR/transpos_element_annotation/auto/${WORK[$SLURM_ARRAY_TASK_ID]}/parsed.fa.rexdb-plant.cls.pep | \
    grep Ty1-RT | sed 's/>//' | sed 's/ .\+//' > copia_arab.txt

more $PROJECTDIR/transpos_element_annotation/auto/${WORK[$SLURM_ARRAY_TASK_ID]}/parsed.fa.rexdb-plant.cls.pep | \
    grep Ty3-RT | sed 's/>//' | sed 's/ .\+//' > gypsy_arab.txt 

more $PROJECTDIR/transpos_element_annotation/auto/${WORK[$SLURM_ARRAY_TASK_ID]}/Brassicaceae.cls.pep | \
    grep Ty1-RT | sed 's/>//' | sed 's/ .\+//' > copia_brassic.txt

more $PROJECTDIR/transpos_element_annotation/auto/${WORK[$SLURM_ARRAY_TASK_ID]}/Brassicaceae.cls.pep | \
    grep Ty3-RT | sed 's/>//' | sed 's/ .\+//' > gypsy_brassic.txt

seqkit grep \
    -f copia_arab.txt \
    $PROJECTDIR/transpos_element_annotation/auto/${WORK[$SLURM_ARRAY_TASK_ID]}/parsed.fa.rexdb-plant.cls.pep \
    -o $WORKDIR/copia_arab.fa 

seqkit grep \
    -f gypsy_arab.txt \
    $PROJECTDIR/transpos_element_annotation/auto/${WORK[$SLURM_ARRAY_TASK_ID]}/parsed.fa.rexdb-plant.cls.pep \
    -o $WORKDIR/gypsy_arab.fa 

seqkit grep \
    -f copia_brassic.txt \
    $PROJECTDIR/transpos_element_annotation/auto/${WORK[$SLURM_ARRAY_TASK_ID]}/Brassicaceae.cls.pep \
    -o $WORKDIR/copia_brassic.fa

seqkit grep \
    -f gypsy_brassic.txt \
    $PROJECTDIR/transpos_element_annotation/auto/${WORK[$SLURM_ARRAY_TASK_ID]}/Brassicaceae.cls.pep \
    -o $WORKDIR/gypsy_brassic.fa

#concatenate files 
cat $WORKDIR/copia_arab.fa $WORKDIR/copia_brassic.fa > $WORKDIR/copia.fa
cat $WORKDIR/gypsy_arab.fa  $WORKDIR/gypsy_brassic.fa > $WORKDIR/gypsy.fa

#cancel useless files 
rm gypsy_*
rm copia_*

#for Copia TE
sed -i 's/#.\+//' $WORKDIR/copia.fa
sed -i 's/:/_/g' $WORKDIR/copia.fa

#for Gypsy TE
sed -i 's/#.\+//' $WORKDIR/gypsy.fa
sed -i 's/:/_/g' $WORKDIR/gypsy.fa

clustalo -i $WORKDIR/copia.fa -o $WORKDIR/copia_aligned.fa
clustalo -i $WORKDIR/gypsy.fa -o $WORKDIR/gypsy_aligned.fa

FastTree -out copia.tre $WORKDIR/copia_aligned.fa
FastTree -out gypsy.tre $WORKDIR/gypsy_aligned.fa

Brassicaceae_tsv=$PROJECTDIR/transpos_element_annotation/auto/${WORK[$SLURM_ARRAY_TASK_ID]}/Brassicaceae.cls.tsv
arab_tsv=$PROJECTDIR/transpos_element_annotation/auto/${WORK[$SLURM_ARRAY_TASK_ID]}/parsed.fa.rexdb-plant.cls.tsv

#for i in ${Brassicaceae_tsv} ${arab_tsv} ; do
#    for j in 'Athila RT' 'CRM RT' 'Reina RT' 'Retand RT' 'Tekay RT' 'Galadriel RT' ; do 
#        grep -e $j $i | cut -f 1 | sed 's/:/_/' | sed 's/$/ #9606/' > $WORKDIR/Retand_ID.txt
#    done
#done
