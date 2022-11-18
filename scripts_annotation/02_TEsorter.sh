#!/usr/bin/env bash

#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00
#SBATCH --job-name=TE_sorter
#SBATCH --mail-user=marco.visani@students.unibe.ch
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=pall
#SBATCH --array=0,1

#load modules
module load UHTS/Analysis/BEDTools/2.29.2

#set variables 
WORK=('canu' 'flye')
genome=('canu_pilon.fasta' 'flye_pilon.fasta')
ASSEMBLIES=('/data/users/mvisani/genome_assembly/polishing/canu/canu_pilon.fasta' \
    '/data/users/mvisani/genome_assembly/polishing/flye/flye_pilon.fasta')

PROJECTDIR=/data/users/mvisani/annotation_genome
CONTAINER=/data/courses/assembly-annotation-course/containers2

OUTDIR=$PROJECTDIR/transpos_element_annotation/auto/
WORKDIR=$OUTDIR/${WORK[$SLURM_ARRAY_TASK_ID]}

if ! [ -d $OUTDIR ]; then mkdir $OUTDIR; fi
if ! [ -d $WORKDIR ]; then mkdir $WORKDIR; fi

cd $WORKDIR

rm $WORKDIR/*_edited
rm $WORKDIR/parsed*
awk '$3~/retrotransposon/' $WORKDIR/${genome[$SLURM_ARRAY_TASK_ID]}.mod.EDTA.TEanno.gff3 \
    > $WORKDIR/${genome[$SLURM_ARRAY_TASK_ID]}.mod.EDTA.TEanno.gff3_edited

grep -v LTR $WORKDIR/${genome[$SLURM_ARRAY_TASK_ID]}.mod.EDTA.TEanno.gff3 \
    >> $WORKDIR/${genome[$SLURM_ARRAY_TASK_ID]}.mod.EDTA.TEanno.gff3_edited

sed 's/\(ID=\)\(\w*\)\(;.*\)/\2/g' $WORKDIR/${genome[$SLURM_ARRAY_TASK_ID]}.mod.EDTA.TEanno.gff3_edited | \
    awk '{printf ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $9, $4, $5, $6, $7, $8, $3)}' \
    > $WORKDIR/parsed.tab

sed -i 's/_pi\s/_pilon\t/g' $WORKDIR/parsed.tab
sed -i 's/_pil\s/_pilon\t/g' $WORKDIR/parsed.tab
sed -i 's/_pilo\s/_pilon\t/g' $WORKDIR/parsed.tab

bedtools getfasta \
    -s \
    -name \
    -fi ${ASSEMBLIES[$SLURM_ARRAY_TASK_ID]} \
    -bed $WORKDIR/parsed.tab \
    -fo $WORKDIR/parsed.fa


singularity exec \
    --bind $CONTAINER \
    --bind $PROJECTDIR \
$CONTAINER/TEsorter_1.3.0.sif \
TEsorter $WORKDIR/parsed.fa \
    -db rexdb-plant \
    -p $SLURM_CPUS_PER_TASK

#on whole genome 
singularity exec \
    --bind $CONTAINER \
    --bind $PROJECTDIR \
    --bind /data/courses/assembly-annotation-course \
$CONTAINER/TEsorter_1.3.0.sif \
TEsorter /data/courses/assembly-annotation-course/Brassicaceae_repbase_all_march2019.fasta \
    -pre Brassicaceae \
    -p $SLURM_CPUS_PER_TASK \
    -db rexdb-plant
    
    