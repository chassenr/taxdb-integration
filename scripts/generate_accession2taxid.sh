#!/bin/bash

ACCNOS=$(echo $1 | sed 's/\s/\t/' | cut -f1)
TAXID=$(echo $1 | sed 's/\s/\t/' | cut -f2)
GENDIR=$2

if ls $GENDIR/$ACCNOS* 1> /dev/null 2>&1
then
  FILE=$(ls $GENDIR/$ACCNOS*)
  zgrep '^>' $FILE | sed -e "s/^>//" -e "s/$/\t$TAXID/"
fi


