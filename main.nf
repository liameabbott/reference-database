/*
 *
 *
 */

// enable DSL 2 syntax
nextflow.enable.dsl=2

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
if ( params.data_ingestion ) {
    if ( !params.fasta_url ) {
        error "Missing 'fasta_url' parameter required by reference data ingestion."
    }
    if ( !params.gtf_url ) {
        error "Missing 'gtf_url' parameter required by reference data ingestion."
    }
}

// construct path to specific database directory
params.genomes_directory = [
    params.database_directory,
    "genomes",
    params.source_database,
    "release-${params.release}",
    params.species,
    params.assembly
].join("/")

// include all the necessary modules
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
    star_create_index;
    rsem_prepare_reference
} from './modules.nf'

workflow {

    // data ingestion routine
    if ( params.data_ingestion ) {

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
        gtf = get_reference_gtf(params.gtf_url, fasta_index[0])

        // parse GTF into reduced and refFlat formats
        reduced_gtf = reduce_gtf(gtf[0], dict)[0]
        refflat = gtf_to_refFlat(gtf[0], dict)[0]

        // parse annotations into BED format
        bed = gtf_to_bed(gtf[0])[0]

        // extract regions into individual BED files
        genes_bed = extract_genes(bed)[0]
        exons_bed = extract_exons(bed)[0]
        cds_bed = extract_CDS(bed)[0]
        utr_bed = extract_UTR(exons_bed, cds_bed)[0]
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
            utr_bed,
            rrna_bed,
            mt_bed,
            intronic_bed,
            intergenic_bed)

        // gzipped reference files for downstream work
        fasta = normalized_fasta[1]
        gtf = gtf[1]

    } else {
        // read files from database if not running data ingestion
        fasta = Channel.fromPath(
            "${params.genomes_directory}/fasta/reference.fa.gz")
        gtf = Channel.fromPath(
            "${params.genomes_directory}/gtf/reference.gtf.gz")
    }

    // STAR indexing routine
    if ( params.star_index ) {
        star_index = star_create_index(fasta, gtf)
    }

    // RSEM reference routine
    if ( params.rsem_reference ) {
        rsem_reference = rsem_prepare_reference(fasta, gtf)
    }

}
