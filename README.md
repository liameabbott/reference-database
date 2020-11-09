# reference-database

This repository contains pipeline code to ingest and process raw source data from three main databases: Ensembl, GENCODE, and RefSeq.

Use this Nextflow pipeline and its configuration parameters to pull reference data from those source databases and create metadata files, such as interval files and alignment indices.

### Running the pipeline

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

### Parameters

#### Core
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

#### Reference data ingestion
If both `fasta_url` and `gtf_url` are not `null`, the pipeline will fetch the reference sequence and annotation from the external database and store them as well as other derived metadata files in `<database_directory>/genomes`. See details below.

* `fasta_url` (optional, default `null`):
  * The URL of the reference sequence to download (e.g. `"ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"`).
* `gtf_url` (optional, default `null`):
  * The URL of the reference annotation file to download (e.g. `"ftp://ftp.ensembl.org/pub/release-101/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.101.gtf.gz"`).

#### STAR alignment index
If `star_read_lengths` is not null, the pipeline will create STAR indices during execution. See details below.

* `star_read_lengths` (optional, default `null`): 
  * A string of comma-separated integers representing the read lengths of STAR indices to create. For example, `--star_read_lengths 50,60,100` will create three STAR indices, with STAR parameter `sjdbOverhang` equal to `read_length - 1` (so `49`, `59`, and `99`, respectively).
* `star_genomeSAindexNbases` (optional, default `14`): 
  * This value is passed through to STAR's `genomeSAindexNbases` parameter.
  
### Pipeline components
  
The pipeline consists of separate components that may be run together in one run of the pipeline, or separately as needed (as long as the necessary reference files already exist in the database).
  
#### Reference data ingestion
  
If both `fasta_url` and `gtf_url` are valid URLs, the reference data ingestion component will be run during pipeline execution. 
  
#### STAR alignment index
  
If `star_read_lengths` is a string of comma-separated integers (e.g. `"50,60,100"`), then STAR indices corresponding to those read lengths will be created during pipeline execution.
  
This component requires that either the reference data ingestion component is also run during the same execution, or that valid reference sequence and annotation files (as defined by the `source_database`, `release`, `species`, and `assembly` parameters) already exist in the database.

