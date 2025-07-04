#!/bin/bash
INPUT="/path/to/FASTQ/files"
OUTPUTP="02-trimming/paired"
OUTPUTUNP="02-trimming/unpaired"

time while read SAMP \n
            do
            java -jar trimmomatic.39/trimmomatic-0.39.jar PE $INPUT/${SAMP}_R1.fastq $INPUT/${SAMP}_R2.fastq \
            $OUTPUTP/${SAMP}_1_trimmomatic_paired.fastq \
            $OUTPUTUNP/${SAMP}_1_trimmomatic_unpaired.fastq \
            $OUTPUTP/${SAMP}_2_trimmomatic_paired.fastq \
            $OUTPUTUNP/${SAMP}_2_trimmomatic_unpaired.fastq \
            ILLUMINACLIP:trimmomatic.39/adapters/TruSeq3-PE-2.fa:2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        done < sample-names.txt