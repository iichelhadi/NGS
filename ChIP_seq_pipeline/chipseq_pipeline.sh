#!/bin/bash

# ------------------------------------------------------------------------------
# Script: chipseq_pipeline.sh
# Author: Dr. Elhadi Iich
# Affiliation: NUS Yong Loo Lin School of Medicine
# Contact: iichelhadi@gmail.com; eiich01@nus.edu.sg
# Created: 2022-07-01
# License: Attribution Required. Do not redistribute or modify without written permission.
# Description: ChIP-seq processing pipeline for trimming, alignment, BAM filtering, and bigWig generation.
# ------------------------------------------------------------------------------

# Configurable variables (modify as needed)
SAMPLE_LIST="samples.txt"
RAW_DATA_DIR="/path/to/raw_fastq"
TRIM_GALORE="/path/to/trim_galore"
BOWTIE2_INDEX="/path/to/bowtie2/genome"
GTF_FILE="/path/to/genes.gtf"
N_THREADS=24

# Create output directories
mkdir -p trim_fq/ SAM_files/ BAM_files/ bigWig/ bigWig/norm/ vis/

# Main loop
while read line; do

    echo "########## Running TrimGalore for $line ##########"
    if [ ! -e trim_fq/${line}_1_val_1.fq.gz ] || [ ! -e trim_fq/${line}_2_val_2.fq.gz ]; then
        "$TRIM_GALORE" --paired \
            "$RAW_DATA_DIR/${line}/${line}_1.fq.gz" "$RAW_DATA_DIR/${line}/${line}_2.fq.gz" \
            -o trim_fq/
    fi

    echo "########## Running Bowtie2 ##########"
    if [ ! -e SAM_files/${line}.sam ]; then
        bowtie2 -p "$N_THREADS" --local -x "$BOWTIE2_INDEX" \
            -1 trim_fq/${line}_1_val_1.fq.gz \
            -2 trim_fq/${line}_2_val_2.fq.gz \
            -S SAM_files/${line}.sam
    fi

    echo "########## Converting SAM to BAM ##########"
    if [ ! -e BAM_files/${line}.bam ]; then
        samtools view -Shb SAM_files/${line}.sam > BAM_files/${line}.bam
    fi

    echo "########## Sorting BAM ##########"
    if [ ! -e BAM_files/${line}_sorted.bam ]; then
        sambamba sort -t "$N_THREADS" -o BAM_files/${line}_sorted.bam BAM_files/${line}.bam
    fi

    echo "########## Filtering Aligned BAM ##########"
    if [ ! -e BAM_files/${line}_aln.bam ]; then
        sambamba view -h -t "$N_THREADS" \
            -f bam \
            -F "[XS] == null and not unmapped and not duplicate" \
            BAM_files/${line}_sorted.bam > BAM_files/${line}_aln.bam
    fi

    echo "########## Indexing BAM ##########"
    samtools index BAM_files/${line}_aln.bam

    echo "########## Generating bigWig for $line ##########"
    if [ ! -e bigWig/${line}.bw ]; then
        bamCoverage -b BAM_files/${line}_aln.bam \
            -o bigWig/${line}.bw \
            --binSize 20 \
            --normalizeUsing BPM \
            --smoothLength 60 \
            --extendReads 150 \
            --centerReads \
            -p 12 \
            2> bigWig/${line}_bamCoverage.log
    fi

done < "$SAMPLE_LIST"

# ------------------------------------------------------------------------------
# Optional downstream steps (commented for manual control)
# Uncomment and customize if needed
# ------------------------------------------------------------------------------

# Run MultiQC on FastQC outputs (if available)
# multiqc -f trim_QC/ -o trim_QC/

# MACS2 peak calling (example template)
# macs2 callpeak -t BAM_files/sample1_aln.bam -c BAM_files/control1_aln.bam -f BAM -n sample1_vs_control1 --outdir macs2/

# ComputeMatrix / plotProfile / plotHeatmap — example commented templates included below

# computeMatrix reference-point --referencePoint TSS -b 1000 -a 1000 \
#     -R /path/to/genes.bed \
#     -S bigWig/sample1.bw \
#     --skipZeros -o bigWig/sample1_TSS.gz \
#     -p 12 --outFileSortedRegions bigWig/sample1_TSS.bed

# plotProfile -m bigWig/sample1_TSS.gz -out vis/sample1_profile.png --refPointLabel "TSS"

# plotHeatmap -m bigWig/sample1_TSS.gz -out vis/sample1_heatmap.png --colorMap RdBu


echo "Processed using Elhadi Iich's pipeline – Attribution Required" >> pipeline_metadata.txt
echo "Author: Dr. Elhadi Iich (iichelhadi@gmail.com)" > .pipeline_author
uuidgen > .pipeline_signature


# ------------------------------------------------------------------------------
# Attribution Notice
# This pipeline was developed by Dr. Elhadi Iich (iichelhadi@gmail.com; eiich01@nus.edu.sg).
# Attribution is required. Do not redistribute or modify without written permission.
# A record of this run is stored in `.pipeline_author` and `pipeline_metadata.txt`.
# ------------------------------------------------------------------------------
