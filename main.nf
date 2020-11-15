/*
 *
 *
 */

// enable DSL 2 syntax
nextflow.enable.dsl=2

// core parameters
params.database_directory = null
params.source_database = null
params.release = null
params.species = null
params.assembly = null

// reference data ingestion parameters
params.fasta_url = null
params.gtf_url = null

// STAR alignment index parameters
params.star_read_lengths = null
params.star_genomeSAindexNbases = 14

// check that all required parameters are defined appropriately
if ( !params.database_directory ) {
    error "Missing 'database_directory' parameter."
}
if ( ! (['ensembl', 'refseq', 'gencode'].contains(params.source_database)) ) {
    error "'source_database' parameter must be one of {'ensembl', 'refseq', 'gencode'}."
}
if ( !params.release ) {
    error "Missing 'release' parameter."
}
if ( !params.species ) {
    error "Missing 'species' parameter."
}
if ( !params.assembly ) {
    error "Missing 'assembly' parameter."
}
if ( params.run_reference_data_ingestion ) {
    if ( !params.fasta_url ) {
        error "Missing 'fasta_url' parameter required by reference data ingestion."
    }
    if ( !params.gtf_url ) {
        error "Missing 'gtf_url' parameter required by reference data ingestion."
    }
}

params.genomes_directory = [
    params.database_directory,
    "genomes",
    params.source_database,
    "release_${params.release}",
    params.species,
    params.assembly
].join("/")

include {
    get_reference_fasta;
    extract_primary_assembly;
    normalize_fasta;
    index_fasta;
    create_sequence_dictionary;
    get_reference_gtf;
    reduce_gtf;
    gtf_to_refFlat;
    gtf_to_bed;
    extract_genes;
    extract_exons;
    extract_CDS;
    extract_rRNA_genes;
    extract_MT_genes;
    extract_intronic_regions;
    extract_intergenic_regions;
    bed_to_interval_list;
} from './genomes.nf'

include {
    rsem_prepare_reference;
    rsem_prepare_reference_with_star;
    generate_star_index
} from './indices.nf'

workflow {

    if ( params.run_reference_data_ingestion ) {

        // fetch reference genome sequence
        fasta = get_reference_fasta(params.fasta_url)

        // extract primary assembly if refseq  
        // (unsure how to do this for gencode or ensembl?)
        // can also just provide primary assembly URL
        if ( params.source_database == "refseq" ) {
            fasta = extract_primary_assembly(fasta[0])
        }

        // normalize and index fasta
        normalized_fasta = normalize_fasta(fasta[0])
        fasta_index = index_fasta(normalized_fasta[1])

        // create Picard Tools sequence dictionary
        dict = create_sequence_dictionary(fasta[0])

        // fetch reference genome annotation
        gtf = get_reference_gtf(params.gtf_url)

        // parse GTF into reduced and refFlat formats
        reduced_gtf = reduce_gtf(gtf[0], dict)[0]
        refflat = gtf_to_refFlat(gtf[0], dict)[0]

        // parse annotations into BED format
        bed = gtf_to_bed(gtf[0])[0]

        // extract regions into individual BED files
        genes_bed = extract_genes(bed)[0]
        exons_bed = extract_exons(bed)[0]
        cds_bed = extract_CDS(bed)[0]
        rrna_bed = extract_rRNA_genes(bed)[0]
        mt_bed = extract_MT_genes(bed)[0]
        intronic_bed = extract_intronic_regions(
            genes_bed, exons_bed)[0]
        intergenic_bed = extract_intergenic_regions(
            genes_bed, fasta_index[0])[0]

        // convert the BED files into interval_list format
        bed_to_interval_list(
            dict,
            genes_bed,
            exons_bed,
            cds_bed,
            rrna_bed,
            mt_bed,
            intronic_bed,
            intergenic_bed)

        // gzipped reference files for downstream work
        fasta = normalized_fasta[1]
        gtf = gtf[1]

    } else {
        // read fasta and gtf files from database 
        fasta = Channel.fromPath(
            "${params.genomes_directory}/fasta/reference.fa.gz")
        gtf = Channel.fromPath(
            "${params.genomes_directory}/gtf/reference.gtf.gz")
    }

    // create RSEM indices
    if ( params.run_rsem_prepare_reference ) {
        if ( params.rsem_with_star ) {
            star_sjdb_overhang = Channel.fromList(
                params.rsem_star_sjdb_overhang.toString().replaceAll("\\s", "").tokenize(","))
            rsem_idx = rsem_prepare_reference_with_star(
                fasta, gtf, star_sjdb_overhang)
        } else {
            rsem_idx = rsem_prepare_reference(
                fasta, gtf)
        }
    }

    // create STAR alignment indices
    if ( params.run_star_indexing ) {
        star_sjdb_overhang = Channel
            .fromList(params.star_sjdb_overhang.toString().replaceAll("\\s", "").tokenize(","))
        star_idx = generate_star_index(
            fasta, gtf, star_sjdb_overhang)
    }

}
