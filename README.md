# reference-database

This repository contains pipeline code to ingest and process raw source data from three main databases: Ensembl, GENCODE, and RefSeq.

Use this Nextflow pipeline and its configuration parameters to pull reference data from those source databases and create metadata files, such as interval files and alignment indices.

## Running the pipeline

To run the pipeline, you must have Docker and Nextflow installed on your machine. Then, the base command is simply:
```
nextflow run gates-mri-bioinformatics/reference-database -user <Github username>
```

To supply required and optional parameters to the pipeline, you may use either command-line arguments or a local configuration file.

With command-line arguments:
```
nextflow run gates-mri-bioinformatics/reference-database \
-user <Github username> \
--database_directory s3://gates-mri-bioinformatics/reference-database \
--source_database ensembl \
--release 101 \
--species homo_sapiens \
--assembly GRCh38 \
--fasta_url "ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" \
--gtf_url "ftp://ftp.ensembl.org/pub/release-101/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.101.gtf.gz"
```

The same command, with a configuration file `myparamters.config`:
```
nextflow run gates-mri-bioinformatics/reference-database -c myparameters.config -user <Github username>
```
where `myparamters.config` contains:
```
params {
  database_directory: "s3://gates-mri-bioinformatics/reference-database",
  source_database: "ensembl",
  release: 101,
  species: "homo_sapiens",
  assembly: "GRCh38",
  fasta_url: "ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
  gtf_url: "ftp://ftp.ensembl.org/pub/release-101/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.101.gtf.gz"
}
```

## Parameters

### Core
These parameters are shared across all pipeline execution components.

* `database_directory` (optional, default `"s3://gates-mri-bioinformatics/reference-database"`):
  * The local directory or S3 bucket of your database.
* `source_database` (required): 
  * One of `{"ensembl", "gencode", "refseq"}`. The database from which to pull raw sequence and annotation data.
* `release` (required): 
  * The release or version of the source database (e.g. `101`).
* `species` (required): 
  * The species of the organism (e.g. `"homo_sapiens"`, `"macaca_mulatta"`).
* `assembly` (required): 
  * The assembly build of the species' genome (e.g. `"GRCh38"`, `"Mmul_10"`).

### Reference data ingestion
If both `fasta_url` and `gtf_url` are not `null`, the pipeline will fetch the reference sequence and annotation from the external database and store them as well as other derived metadata files in `<database_directory>/genomes`. See details below.

* `fasta_url` (optional, default `null`):
  * The URL of the reference sequence to download (e.g. `"ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"`).
* `gtf_url` (optional, default `null`):
  * The URL of the reference annotation file to download (e.g. `"ftp://ftp.ensembl.org/pub/release-101/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.101.gtf.gz"`).

### STAR alignment index
If `star_read_lengths` is not null, the pipeline will create STAR indices during execution. See details below.

* `star_read_lengths` (optional, default `null`): 
  * A string of comma-separated integers representing the read lengths of STAR indices to create. For example, `--star_read_lengths 50,60,100` will create three STAR indices, with STAR parameter `sjdbOverhang` equal to `read_length - 1` (so `49`, `59`, and `99`, respectively).
* `star_genomeSAindexNbases` (optional, default `14`): 
  * This value is passed through to STAR's `genomeSAindexNbases` parameter.
  
## Pipeline components
  
The pipeline consists of separate components that may be run together in one run of the pipeline, or separately as needed (as long as the necessary reference files already exist in the database).
  
### Reference data ingestion
  
If both `fasta_url` and `gtf_url` are valid URLs, the reference data ingestion component will be run during execution.

The following reference and metadata files are published in the `<database_directory>`:

#### Reference sequence

* `<source_database>/release-<release>/<species>/<assembly>/fasta/<source_database>.release_<release>.<species>.<assembly>.fa.gz`
  * The reference genome sequence.
* `<source_database>/release-<release>/<species>/<assembly>/fasta/<source_database>.release_<release>.<species>.<assembly>.fa.fai`
  * FASTA index.
* `<source_database>/release-<release>/<species>/<assembly>/fasta/<source_database>.release_<release>.<species>.<assembly>.fa.gzi`
  * Block compression index.
* `<source_database>/release-<release>/<species>/<assembly>/fasta/<source_database>.release_<release>.<species>.<assembly>.dict`
  * Picard Tools reference sequence dictionary. 

#### Reference annotation

* `<source_database>/release-<release>/<species>/<assembly>/gtf/<source_database>.release_<release>.<species>.<assembly>.gtf.gz`
  * Reference annotation file.
* `<source_database>/release-<release>/<species>/<assembly>/gtf/<source_database>.release_<release>.<species>.<assembly>.reduced.gtf.gz`
  * Reduced GTF created with Drop-seq Tools' `ReduceGtf`.
* `<source_database>/release-<release>/<species>/<assembly>/gtf/<source_database>.release_<release>.<species>.<assembly>.refFlat`
  * refFlat file used by Picard Tools, created with `ConvertToRefFlat` tool in Drop-seq Tools package.
  
#### Reference genome regions in BED format

* `<source_database>/release-<release>/<species>/<assembly>/bed/<source_database>.release_<release>.<species>.<assembly>.bed.gz`
  * The main annotation file, in BED format.
* `<source_database>/release-<release>/<species>/<assembly>/bed/<source_database>.release_<release>.<species>.<assembly>.genes.bed.gz`
  * All genes in the reference genome.
* `<source_database>/release-<release>/<species>/<assembly>/bed/<source_database>.release_<release>.<species>.<assembly>.exons.bed.gz`
  * All exons in the reference genome.
* `<source_database>/release-<release>/<species>/<assembly>/bed/<source_database>.release_<release>.<species>.<assembly>.CDS.bed.gz`
  * All CDS regions in the reference genome.
* `<source_database>/release-<release>/<species>/<assembly>/bed/<source_database>.release_<release>.<species>.<assembly>.genes.rRNA.bed.gz`
  * All ribosomal RNA genes in the reference genome.
* `<source_database>/release-<release>/<species>/<assembly>/bed/<source_database>.release_<release>.<species>.<assembly>.genes.MT.bed.gz`
  * All mitochondrial genes in the reference genome.
* `<source_database>/release-<release>/<species>/<assembly>/bed/<source_database>.release_<release>.<species>.<assembly>.intronic.bed.gz`
  * All intronic regions in the reference genome.
* `<source_database>/release-<release>/<species>/<assembly>/bed/<source_database>.release_<release>.<species>.<assembly>.intergenic.bed.gz`
  * All intergenic regions in the reference genome.
  
#### Reference genome regions in interval_list format

* `<source_database>/release-<release>/<species>/<assembly>/interval-list/<source_database>.release_<release>.<species>.<assembly>.interval_list`
  * The main annotation file, in BED format.
* `<source_database>/release-<release>/<species>/<assembly>/interval-list/<source_database>.release_<release>.<species>.<assembly>.genes.interval_list`
  * All genes in the reference genome.
* `<source_database>/release-<release>/<species>/<assembly>/interval-list/<source_database>.release_<release>.<species>.<assembly>.exons.interval_list`
  * All exons in the reference genome.
* `<source_database>/release-<release>/<species>/<assembly>/interval-list/<source_database>.release_<release>.<species>.<assembly>.CDS.interval_list`
  * All CDS regions in the reference genome.
* `<source_database>/release-<release>/<species>/<assembly>/interval-list/<source_database>.release_<release>.<species>.<assembly>.genes.rRNA.interval_list`
  * All ribosomal RNA genes in the reference genome.
* `<source_database>/release-<release>/<species>/<assembly>/interval-list/<source_database>.release_<release>.<species>.<assembly>.genes.MT.interval_list`
  * All mitochondrial genes in the reference genome.
* `<source_database>/release-<release>/<species>/<assembly>/interval-list/<source_database>.release_<release>.<species>.<assembly>.intronic.interval_list`
  * All intronic regions in the reference genome.
* `<source_database>/release-<release>/<species>/<assembly>/interval-list/<source_database>.release_<release>.<species>.<assembly>.intergenic.interval_list`
  * All intergenic regions in the reference genome.
  
#### RSEM index files

* `<source_database>/release-<release>/<species>/<assembly>/rsem/<source_database>.release_<release>.<species>.<assembly>.chrlist`
* `<source_database>/release-<release>/<species>/<assembly>/rsem/<source_database>.release_<release>.<species>.<assembly>.grp`
* `<source_database>/release-<release>/<species>/<assembly>/rsem/<source_database>.release_<release>.<species>.<assembly>.idx.fa`
* `<source_database>/release-<release>/<species>/<assembly>/rsem/<source_database>.release_<release>.<species>.<assembly>.n2g.idx.fa`
* `<source_database>/release-<release>/<species>/<assembly>/rsem/<source_database>.release_<release>.<species>.<assembly>.seq`
* `<source_database>/release-<release>/<species>/<assembly>/rsem/<source_database>.release_<release>.<species>.<assembly>.ti`
* `<source_database>/release-<release>/<species>/<assembly>/rsem/<source_database>.release_<release>.<species>.<assembly>.transcripts.fa`

#### STAR alignment index
  
If `star_read_lengths` is a string of comma-separated integers (e.g. `"50,60,100"`), then STAR indices corresponding to those read lengths will be created during execution.
  
This component requires that either the reference data ingestion component is also run during the same execution, or that valid reference sequence and annotation files (as defined by the `source_database`, `release`, `species`, and `assembly` parameters) already exist in the database.

This component creates the file:

* `<source_database>/release-<release>/<species>/<assembly>/star/<source_database>.release_<release>.<species>.<assembly>.<read_length>_bp.star_idx/`
