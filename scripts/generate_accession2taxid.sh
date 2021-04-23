#!/bin/bash

ACCNOS=$1
TAXID=$2
GENDIR=$3

if ls $GENDIR/$ACCNOS* 1> /dev/null 2>&1
then
  FILE=$(ls $GENDIR/$ACCNOS*)
  zgrep '^>' $FILE | sed -e "s/^>//" -e "s/$/\t$TAXID/"
fi


