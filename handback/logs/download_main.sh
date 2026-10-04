#!/bin/bash
set -u
cd ~/Documents/genept-leakage
dl() { echo "=== $(date -u +%Y-%m-%dT%H:%M:%SZ) START $1 -> $2"; curl -fSL --retry 3 --retry-delay 5 -o "$2" "$1"; echo "=== $(date -u +%Y-%m-%dT%H:%M:%SZ) END rc=$? $2"; }
dl 'https://zenodo.org/records/10833191/files/GenePT_emebdding_v2.zip?download=1' data/downloads/GenePT_emebdding_v2.zip
dl 'https://current.geneontology.org/annotations/goa_human.gaf.gz' data/raw/go/current/goa_human.gaf.gz
dl 'https://current.geneontology.org/ontology/go-basic.obo' data/raw/go/current/go-basic.obo
dl 'https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt' data/raw/hgnc_complete_set.txt
dl 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz' data/raw/Homo_sapiens.gene_info.gz
dl 'https://raw.githubusercontent.com/jingcheng-du/Gene2vec/master/pre_trained_emb/gene2vec_dim_200_iter_9.txt' data/raw/gene2vec_dim_200_iter_9.txt
echo "=== $(date -u +%Y-%m-%dT%H:%M:%SZ) START gene2pubmed stream"
curl -sSL https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz | gunzip -c | awk -F'\t' 'NR==1 || $1=="9606"' > data/raw/gene2pubmed_human.tsv
echo "=== $(date -u +%Y-%m-%dT%H:%M:%SZ) END gene2pubmed pipestatus=${PIPESTATUS[*]}"
wc -l data/raw/gene2pubmed_human.tsv
echo ALL DONE
