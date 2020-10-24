#!/bin/bash

FAA=$1
MAP=$2


# get taxid
TAXID=$(echo $FAA | sed -e 's/_/\t/g' -e 's/\t/_/' | cut -f1 | grep -F -f - $MAP | cut -f2)

# replace fasta header with index_taxid
sed "s/^>/>$TAXID /" $FAA

