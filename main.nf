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
if ( params.star_read_lengths ) {
    if ( !params.star_genomeSAindexNbases ) {
        error "Missing 'star_genomeSAindexNbases' parameter."
    }
}

params.basename = [
    params.source_database,
    "release_${params.release}",
    params.species,
    params.assembly
].join(".") 

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
    bed_to_interval_list
} from './genomes.nf'

process generate_star_index {
    cpus Runtime.runtime.availableProcessors()
 
    publishDir "${params.database_directory}/star-index/", \
        mode: "copy", overwrite: true //, \
        //saveAs: { filename -> "${basename}.${read_length}bp") }
    //publishDir "${params.database_directory}/star-index", \
    //    pattern: "reference.star_idx.cmd", \
    //    mode: "move", overwrite: true, \
    //    saveAs: { filename -> "${basename}.${read_length}bp.star_idx.cmd" }

    input:
    path(reference_fasta)
    path(reference_gtf)
    val(read_length)

    output:
    path("reference.star_idx.cmd")
    path("reference.star_idx/Genome")
    path("reference.star_idx/Log.out")
    path("reference.star_idx/SA")
    path("reference.star_idx/SAindex")
    path("reference.star_idx/chrLength.txt")
    path("reference.star_idx/chrName.txt")
    path("reference.star_idx/chrNameLength.txt")
    path("reference.star_idx/chrStart.txt")
    path("reference.star_idx/exonGeTrInfo.tab")
    path("reference.star_idx/exonInfo.tab")
    path("reference.star_idx/geneInfo.tab")
    path("reference.star_idx/genomeParameters.txt")
    path("reference.star_idx/sjdbInfo.txt")
    path("reference.star_idx/sjdbList.fromGTF.out.tab")
    path("reference.star_idx/sjdbList.out.tab")
    path("reference.star_idx/transcriptInfo.tab")

    """
    gunzip -c ${reference_fasta} > reference.fa
    gunzip -c ${reference_gtf} > reference.gtf
    cmd=\$(cat <<-EOF
    STAR \
        --runThreadN ${task.cpus} \
        --runMode genomeGenerate \
        --genomeDir reference.star_idx/ \
        --genomeFastaFiles reference.fa \
        --sjdbGTFfile reference.gtf \
        --sjdbOverhang \$((${read_length} - 1)) \
        --genomeSAindexNbases ${params.star_genomeSAindexNbases} \
        --limitGenomeGenerateRAM 200000000000
    EOF
    )
    eval \$cmd
    echo "reference.fa=${params.fasta_url}" >> reference.star_idx.cmd
    echo "reference.gtf=${params.gtf_url}" >> reference.star_idx.cmd
    echo \$cmd >> reference.star_idx.cmd 
    cp -r reference.star_idx/ ${basename}.tester.star_idx/
    """
}

def run_star_indexing(fasta, gtf, read_lengths) {
    star_read_lengths = Channel
        .fromList(read_lengths.toString().replaceAll("\\s", "").tokenize(","))
    star_index = generate_star_index(
        fasta, gtf, star_read_lengths)
}

workflow {

    if ( params.fasta_url && params.gtf_url ) {
        // fetch, normalize, index reference genome sequence
        fasta = get_reference_fasta(params.fasta_url)
        normalized_fasta = normalize_fasta(fasta[0])
        fasta_index = index_fasta(normalized_fasta)
        dict = create_sequence_dictionary(fasta[0])

        // fetch reference genome annotation, and parse into 
        // reduced and refFlat formats
        gtf = get_reference_gtf(params.gtf_url)[0]
        reduced_gtf = reduce_gtf(gtf, dict)[0]
        refflat = gtf_to_refFlat(gtf, dict)[0]

        // parse annotations into BED format and extract
        // regions in individual BED files
        bed = gtf_to_bed(gtf)[0]
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

        if ( params.star_read_lengths ) {
            //run_star_indexing(
            //    fasta, gtf, params.star_read_lengths)
        }
    
    } else if ( params.star_read_lengths ) {
        // read fasta and gtf files from database 
        fasta = Channel.fromPath(
            "${params.genomes_directory}/fasta/${params.basename}.fa.gz")
        gtf = Channel.fromPath(
            "${params.genomes_directory}/gtf/${params.basename}.gtf.gz")

        // create STAR indices with specified read lengths
        //run_star_indexing(
        //    fasta, gtf, params.star_read_lengths)
    }

}
