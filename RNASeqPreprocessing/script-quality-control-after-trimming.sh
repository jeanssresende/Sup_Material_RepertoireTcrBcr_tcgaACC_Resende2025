#!/bin/bash
INPUT="02-trimming/paired"
OUTPUT="03-qualityControl/fastqc"
time while read SAMP \n
        do
            echo "Processing sample ${SAMP}"
            fastqc $INPUT/${SAMP}/* -o $OUTPUT
        done < sample-names.txt

multiqc $OUTPUT/* -o 03-qualityControl/multiqc