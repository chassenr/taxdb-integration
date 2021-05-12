#!/bin/bash

ACCNOS=$(echo $1 | sed 's/\s/\t/' | cut -f1)
TAXID=$(echo $1 | sed 's/\s/\t/' | cut -f2)
LIST=$2

FILE=$(grep "$ACCNOS" $LIST)
if [[ $FILE != "" ]]
then
  zgrep '^>' $FILE | sed -e "s/^>//" -e "s/\s.*$//" -e "s/$/\t$TAXID/"
fi


