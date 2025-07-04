#!/bin/bash
INPUT="/path/to/FASTQ/files"
OUTPUT="01-qualityControl/fastqc"
time while read SAMP \n
        do
            echo "Processing sample ${SAMP}"
            fastqc $INPUT/${SAMP}/* -o $OUTPUT
        done < sample-names.txt

multiqc $OUTPUT/* -o 01-qualityControl/multiqc