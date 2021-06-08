#!/bin/bash
nextflow run /shared/bioinformatics/reference-database/main.nf \
-profile bio \
--source_database ensembl \
--release 103 \
--species homo_sapiens \
--assembly GRCh38 \
--data_ingestion \
--star_index \
--rsem_reference \
--fasta_url http://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
--gtf_url http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.gtf.gz
 
