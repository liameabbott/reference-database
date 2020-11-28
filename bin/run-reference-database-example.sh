#!/bin/bash
nextflow run gates-mri-bioinformatics/reference-database \
-user liameabbott \  # Github user name
-profile bio \  # execution profile
--source_database ensembl \
--release 101 \
--species homo_sapiens \
--assembly GRCh38 \
--data_ingestion \
--star_index \
--rsem_reference
