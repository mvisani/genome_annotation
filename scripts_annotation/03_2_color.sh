#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=512MB
#SBATCH --time=01:00:00
#SBATCH --job-name=te_dynamics.phylogenetic.color
#SBATCH --mail-user=marco.visani@students.unibe.ch
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=pall
#SBATCH --array=0,1

ASSDIR=$1
RTs=('Ty1-RT' 'Ty3-RT')
RT=${RTs[$SLURM_ARRAY_TASK_ID]}

INPDIR=/data/users/mvisani/annotation_genome/transpos_element_annotation/auto/$ASSDIR
OUTDIR=/data/users/mvisani/annotation_genome/transpos_element_annotation/phylogeny/auto/$ASSDIR/$RT

mkdir $OUTDIR
# COLORBLIND PALETTE
#PALETTE=("#3366FF" "#99CCFF" "#CC99FF" "#CCCCCC" "#6699CC" "#0099CC" "#33CCCC" \
#  "#006699" "#003366" "#000066" "#333399" "#333333" "#666699" "#99CCCC" "#669999" "#006666" \
#  "#003333" "#000033" "#333366" "#000000")

PALETTE=('#3cb44b' '#ffe119' '#4363d8' '#f58231' '#911eb4' '#42d4f4' '#f032e6' \
      '#bfef45' '#fabed4' '#469990' '#dcbeff' '#9A6324' '#fffac8' '#800000' '#aaffc3' \
      '#808000' '#ffd8b1' '#000075' '#a9a9a9' '#000000' '#ff0000')

REF_CLADES=$( awk '! /#/ {print $4}' ${INPDIR}/Brassicaceae.cls.tsv | sort | uniq )
CUR_CLADES=$( awk '! /#/ {print $4}' ${INPDIR}/parsed.fa.rexdb-plant.cls.tsv | sort | uniq )
CLADES=($(printf "%s\n" "${REF_CLADES[@]}" "${CUR_CLADES[@]}" | sort | uniq))

for i in "${!CLADES[@]}"; do
  if [ "${CLADES[$i]}" == "unknown" ]; then
    unset 'CLADES[$i]'
  fi
done

rm -rf $OUTDIR/color_clade.map
for i in "${!CLADES[@]}"; do
  echo "${PALETTE[$i]} ${CLADES[$i]}" >> $OUTDIR/color_clade.map
done

sed -i -e "s/$RT//g" $OUTDIR/*.tre
cp /data/users/agiottonini/assembly_annotation_pipeline/src/te_dynamics/phylogenetic/color_template.txt $OUTDIR/colors.txt
for i in "${!CLADES[@]}"; do
    grep -e ${CLADES[$i]} ${INPDIR}/parsed.fa.rexdb-plant.cls.tsv \
        | cut -f 1 \
        | sed 's/[:()]/_/g' \
        | sed 's/__/_/g' \
        | sed -e "s/$/ ${PALETTE[$i]}/g" \
        | grep '#' >> $OUTDIR/colors.txt
    cat $OUTDIR/colors.txt
    grep -e ${CLADES[$i]} ${INPDIR}/Brassicaceae.cls.tsv \
        | cut -f 1 \
        | sed 's/[:()]/_/g' \
        | sed 's/__/_/g' \
        | sed -e "s/$/ ${PALETTE[$i]}/g" \
        | grep '#' >> $OUTDIR/colors.txt
done
cp /data/users/agiottonini/assembly_annotation_pipeline/src/te_dynamics/phylogenetic/labels_template.txt $OUTDIR/labels.txt
for i in "${!CLADES[@]}"; do
    grep -e ${CLADES[$i]} $INPDIR/parsed.fa.rexdb-plant.cls.tsv \
        | cut -f 1 \
        | sed 's/[:()]/_/g' \
        | sed 's/__/_/g' \
        | sed -e "s/$/ ${CLADES[$i]} clade/g" \
        >> $OUTDIR/labels.txt
    grep -e ${CLADES[$i]} ${INPDIR}/Brassicaceae.cls.tsv \
        | cut -f 1 \
        | sed 's/[:()]/_/g' \
        | sed 's/__/_/g' \
        | sed -e "s/$/ ${CLADES[$i]} clade/g" \
        >> $OUTDIR/labels.txt
done
# LEGENDS
CLADES_STR=$(printf "%s " "${CLADES[@]}")
PALETTE_STR=$(printf "%s " "${PALETTE[@]}")
SHAPE_STR=$(printf "%0.s1 " $(seq 1 20))
sed -i -e "s/#LEGEND_TITLE/LEGEND_TITLE TEST/g" $OUTDIR/colors.txt
sed -i -e "s/#LEGEND_POSITION_X.*//g" $OUTDIR/colors.txt
sed -i -e "s/#LEGEND_POSITION_Y.*//g" $OUTDIR/colors.txt
sed -i -e "s/#LEGEND_SHAPES.*/LEGEND_SHAPES $SHAPE_STR/g" $OUTDIR/colors.txt
sed -i -e "s/#LEGEND_LABELS.*/LEGEND_LABELS $CLADES_STR/g" $OUTDIR/colors.txt
sed -i -e "s/#LEGEND_COLORS.*/LEGEND_COLORS $PALETTE_STR/g" $OUTDIR/colors.txt

sed -i -E "s/_contig\S+/_/g" $OUTDIR/colors.txt
sed -i -E "s/_scaffold\S+/_/g" $OUTDIR/colors.txt
sed -i -E "s/_tig\S+/_/g" $OUTDIR/colors.txt

sed -i -E "s/_contig\S+/_/g" $OUTDIR/labels.txt
sed -i -E "s/_scaffold\S+/_/g" $OUTDIR/labels.txt
sed -i -E "s/_tig\S+/_/g" $OUTDIR/labels.txt