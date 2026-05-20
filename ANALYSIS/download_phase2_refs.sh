#!/bin/bash

# Define the target directory
TARGET_DIR="refs"

# 1. Create the directory if it doesn't exist
mkdir -p "$TARGET_DIR"

echo "Starting Phase 2 downloads..."

# 2. Download Human Gene Annotation (GTF)
# Required for: Mapping Ensembl IDs to Gene Symbols in the pipeline
echo "Downloading Human GRCh38 Annotation (GTF)..."
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.annotation.gtf.gz \
  -O "$TARGET_DIR/GRCh38.primary_assembly.annotation.gtf.gz"

# 3. Download MSigDB Hallmark Gene Sets (GMT)
# Required for: GSEA Pathway Analysis (e.g., Hypoxia, Apoptosis, EMT)
echo "Downloading MSigDB Hallmark Gene Sets (GMT)..."
wget https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt \
  -O "$TARGET_DIR/h.all.v2023.2.Hs.symbols.gmt"

echo "Phase 2 downloads complete. Files saved to $TARGET_DIR"
