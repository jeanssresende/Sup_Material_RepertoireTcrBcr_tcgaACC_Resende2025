#!/bin/bash
TRUST4="path/to/TRUST4/directory"
SAMPLES="02-trimming/paired"

while read SAMP \n
    do
    echo "preprocessing sample ${SAMP}"
    $TRUST4/./run-trust4 -f $TRUST4/human_IMGT+C.fa --ref $TRUST4/human_IMGT+C.fa \
    -1 $SAMPLES/${SAMP}_1_trimmomatic_paired.fastq -2 $SAMPLES/${SAMP}_2_trimmomatic_paired.fastq \
    -o ${SAMP} --od trust4
    done < sample-names.txt

cp trust4/*_report.tsv cdr3