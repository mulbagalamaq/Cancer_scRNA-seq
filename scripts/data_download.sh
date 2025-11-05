#!/bin/bash

# Script to download scRNA-seq data for cancer cell line analysis
# This script downloads data from SRA/ENA for multiple cancer cell lines

# Create data directory if it doesn't exist
mkdir -p ../data/raw
mkdir -p ../data/processed

# Set base directory
DATA_DIR="../data/raw"

echo "Starting data download for single-cell RNA-seq analysis..."
echo "=========================================================="

# Note: Replace these with actual SRA accession numbers for your analysis
# Example SRA numbers (these are placeholders - replace with real accessions)
# SRA_IDS=("SRR12345678" "SRR12345679" "SRR12345680")

# Download using sra-tools (requires prefetch and fasterq-dump)
# Uncomment and modify the following section once you have SRA IDs:

# for SRA_ID in "${SRA_IDS[@]}"; do
#     echo "Downloading $SRA_ID..."
#     prefetch $SRA_ID
#     fasterq-dump $SRA_ID -O $DATA_DIR
# done

# Alternative: Download from GEO using wget or curl
# Example for GEO dataset GSEXXXXXX:
# wget -O $DATA_DIR/GSEXXXXXX.tar.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSEXXXXXX&format=file"

echo ""
echo "Data download script ready."
echo "Please update this script with actual SRA/GEO accession numbers."
echo ""
echo "For public datasets, consider:"
echo "1. Cancer Cell Line Encyclopedia (CCLE) scRNA-seq data"
echo "2. GEO datasets with 'single cell' and 'cancer cell line' keywords"
echo "3. 10x Genomics public datasets"
echo ""
echo "Example GEO search:"
echo "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSEXXXXXX"

