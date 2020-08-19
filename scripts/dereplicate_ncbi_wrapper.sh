#!/bin/bash

LIN=$1
TL_MAIN=$2
TL_SUB=$3
ANI_main=$4
ANI_sub=$5
DEREP=$6
TAX=$7
DIR=$8
CPU=$9

# if sub lineage specified, split taxonomy file and run dereplication separately
if [[ $LIN != "" ]]
then
  grep -v "$LIN" $TAX > $DIR/assembly_taxonomy_main.txt
  grep "$LIN" $TAX > $DIR/assembly_taxonomy_sub.txt
else
  cp $TAX $DIR/assembly_taxonomy_main.txt
fi

# process main taxonomy genomes
if [[ $TL_MAIN -eq 7 ]]
then
  mv $DIR/assembly_taxonomy_main.txt $DIR/assembly_taxonomy_main_tl.txt
else
  cut -f2 $DIR/assembly_taxonomy_main.txt | cut -d';' -f -$TL_MAIN | paste <(cut -f1 $DIR/assembly_taxonomy_main.txt) - > $DIR/assembly_taxonomy_main_tl.txt
fi

# run dereplication for main lineages
$DEREP --threads $CPU --threshold $ANI_main $DIR/genomes/ $DIR/derep_genomes/ $DIR/assembly_taxonomy_main_tl.txt

# process sub taxonomy genomes
if [[ $LIN != "" ]]
then
  if [[ $TL_SUB -eq 7 ]]
  then
    mv $DIR/assembly_taxonomy_sub.txt $DIR/assembly_taxonomy_sub_tl.txt
  else
    cut -f2 $DIR/assembly_taxonomy_sub.txt | cut -d';' -f -$TL_SUB | paste <(cut -f1 $DIR/assembly_taxonomy_sub.txt) - > $DIR/assembly_taxonomy_sub_tl.txt
  fi
  # run dereplication for sub lineages
  $DEREP --threads $CPU --threshold $ANI_sub $DIR/genomes/ $DIR/derep_genomes/ $DIR/assembly_taxonomy_sub_tl.txt
fi
