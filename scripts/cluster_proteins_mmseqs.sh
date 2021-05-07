#!/bin/bash

TAXON=$1
TAXFILE=$2
PDIR=$3
TMPDIR=$4
COV=$5
CMODE=$6
MINID=$7

TAXONHASH=$(echo "${TAXON}" | sha1sum | head -c 40)

INFILES=$(echo "${TAXON}" | grep -F -f - ${TAXFILE} | cut -f1 | grep -F -f - <(find ${PDIR}) | tr '\n' ' ')

mmseqs easy-cluster ${INFILES} ${TMPDIR}/${TAXONHASH} ${TMPDIR}/tmp_${TAXONHASH} -c ${COV} --cov-mode ${CMODE} --min-seq-id ${MINID} --threads 1 >${TMPDIR}/${TAXONHASH}.log 2>&1

cat ${TMPDIR}/${TAXONHASH}_rep_seq.fasta

