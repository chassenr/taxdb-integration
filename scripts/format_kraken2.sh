#!/bin/bash

ACCNOS=$(echo $1 | sed 's/\s/\t/' | cut -f1)
TAXID=$(echo $1 | sed 's/\s/\t/' | cut -f2)
GENDIR=$2

FILE=$(ls $GENDIR/$ACCNOS*)

zcat $FILE | sed "/^>/s/$/\|kraken:taxid\|$TAXID/" | dustmasker -in - -outfmt fasta

