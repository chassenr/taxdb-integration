#!/bin/bash

TAXON=$1
TAXFILE=$2
PDIR=$3
TMPDIR=$4
COV=$5
CMODE=$6
MINID=$7

TAXONHASH=$(echo "${TAXON}" | sha1sum | head -c 40)

ACCNOS=$(echo "${TAXON}" | sed 's/[^;]$/&;/' | grep -F -f - <(sed 's/[^;]$/&;/' ${TAXFILE}) | cut -f1 | tr '\n' ' ')

for i in $ACCNOS
do
  INFILE=$(echo "${i}" | grep -F -f - <(find ${PDIR}))
  zcat $INFILE | sed -e "/^>/s/\s.*$//" -e "/^>/s/^>/>${i}_/"
done >> ${TMPDIR}/in_${TAXONHASH}.faa

mmseqs easy-cluster ${TMPDIR}/in_${TAXONHASH}.faa ${TMPDIR}/${TAXONHASH} ${TMPDIR}/tmp_${TAXONHASH} -c ${COV} --cov-mode ${CMODE} --min-seq-id ${MINID} --threads 1 >${TMPDIR}/${TAXONHASH}.log 2>&1

cat ${TMPDIR}/${TAXONHASH}_rep_seq.fasta

