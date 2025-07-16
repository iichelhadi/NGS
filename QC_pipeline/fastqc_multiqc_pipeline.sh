#!/bin/bash

# ------------------------------------------------------------------------------
# Script: fastqc_multiqc_pipeline.sh
# Author: Dr. Elhadi Iich
# Affiliation: NUS Yong Loo Lin School of Medicine
# Contact: iichelhadi@gmail.com; eiich01@nus.edu.sg
# Created: 2022-07-01
# License: Attribution Required. Do not redistribute or modify without written permission.
# Description: Run FastQC on paired-end trimmed FASTQ files and summarize with MultiQC
# ------------------------------------------------------------------------------

# Define variables
SAMPLE_LIST="samples.txt"
TRIM_DIR="trim_fq"
QC_DIR="trim_QC"
N_THREADS=8  # Modify as needed

# Create output directory
mkdir -p "$QC_DIR"

# Main processing loop
while read line; do
    echo "######################## Running FastQC for $line ####################"
    if [ ! -e "${QC_DIR}/${line}_trimmed_fastqc.zip" ]; then
        fastqc -t "$N_THREADS" -o "$QC_DIR" \
            "${TRIM_DIR}/${line}_1_val_1.fq.gz" "${TRIM_DIR}/${line}_2_val_2.fq.gz"
    fi
done < "$SAMPLE_LIST"

# Run MultiQC
echo "######################## Running MultiQC ###########################"
multiqc -f "$QC_DIR" -o "$QC_DIR"


echo "Processed using Elhadi Iich's pipeline – Attribution Required" >> pipeline_metadata.txt
echo "Author: Dr. Elhadi Iich (iichelhadi@gmail.com)" > .pipeline_author
uuidgen > .pipeline_signature


# ------------------------------------------------------------------------------
# Attribution Notice
# This pipeline was developed by Dr. Elhadi Iich (iichelhadi@gmail.com; eiich01@nus.edu.sg).
# Attribution is required. Do not redistribute or modify without written permission.
# A record of this run is stored in `.pipeline_author` and `pipeline_metadata.txt`.
# ------------------------------------------------------------------------------
