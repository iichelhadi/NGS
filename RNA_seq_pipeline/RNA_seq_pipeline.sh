#!/bin/bash

# ------------------------------------------------------------------------------
# Script: RNA_seq_pipeline
# Author: Dr. Elhadi Iich
# Affiliation: NUS Yong Loo Lin School of Medicine
# Contact: iichelhadi@gmail.com; eiich01@nus.edu.sg
# Created: 2022-07-01
# License: Attribution Required. Do not redistribute or modify without written permission.
# Description: Paired-end RNA-seq pipeline for trimming, alignment, and quantification
# ------------------------------------------------------------------------------


# Define variables (adjust as needed)
SAMPLE_LIST="samples.txt"
RAW_DATA_DIR="/path/to/raw_fastq"
TRIM_GALORE="/path/to/trim_galore"
STAR_INDEX="/path/to/STAR_index"
GTF_FILE="/path/to/genes.gtf"
N_THREADS=24

# Create output directories
mkdir -p trim_fq/ SAM_files/ feature_counts/ BAM_files/

# Main processing loop
while read line; do

    echo "################# Running TrimGalore for $line ###################"
    if [ ! -e trim_fq/${line}_1_val_1.fq.gz ] || [ ! -e trim_fq/${line}_2_val_2.fq.gz ]; then
        "$TRIM_GALORE" --paired \
            "$RAW_DATA_DIR/${line}/${line}_1.fq.gz" \
            "$RAW_DATA_DIR/${line}/${line}_2.fq.gz" \
            -o trim_fq/
    fi

    echo "################# Generating BAM for $line ###################"
    if [ ! -e BAM_files/${line}.bam ]; then
        STAR --runThreadN "$N_THREADS" \
            --genomeDir "$STAR_INDEX" \
            --readFilesIn trim_fq/${line}_1_val_1.fq.gz trim_fq/${line}_2_val_2.fq.gz \
            --readFilesCommand zcat \
            --sjdbGTFfile "$GTF_FILE" \
            --outSAMtype SAM \
            --outSAMorder PairedKeepInputOrder \
            --outFileNamePrefix SAM_files/${line}_

        samtools view -Sb SAM_files/${line}_Aligned.out.sam > BAM_files/${line}.bam
    fi

    echo "################# Running featureCounts for $line ###################"
    if [ ! -e feature_counts/${line}_counts.txt ]; then
        featureCounts -T "$N_THREADS" -p -s 2 -t exon -g gene_id \
            -a "$GTF_FILE" \
            -o feature_counts/${line}_counts.txt BAM_files/${line}.bam
    fi

    rm -f SAM_files/${line}_Aligned.out.sam

done < "$SAMPLE_LIST"


echo "Processed using Elhadi Iich's pipeline – Attribution Required" >> pipeline_metadata.txt
echo "Author: Dr. Elhadi Iich (iichelhadi@gmail.com)" > .pipeline_author
uuidgen > .pipeline_signature


# ------------------------------------------------------------------------------
# Attribution Notice
# This pipeline was developed by Dr. Elhadi Iich (iichelhadi@gmail.com; eiich01@nus.edu.sg).
# Attribution is required. Do not redistribute or modify without written permission.
# A record of this run is stored in `.pipeline_author` and `pipeline_metadata.txt`.
# ------------------------------------------------------------------------------
