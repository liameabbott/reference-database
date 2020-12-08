
process get_reference_fasta {
    publishDir "${params.genomes_directory}/fasta", \
        pattern: "reference.fa.gz.url", \
        mode: "copy", overwrite: true

    input:
    val(fasta_url)

    output:
    path("reference.fa")
    path("reference.fa.gz.url")

    """
    wget -qO- ${fasta_url} | \
    gunzip -c > reference.fa
    
    printf "${fasta_url}\n" > reference.fa.gz.url
    """
}

process extract_primary_assembly {
    input:
    path(fasta)

    output:
    path("primary.fa")

    """
    awk '
        BEGIN { p=0; }
        \$0 ~ /^>/ { if(\$0 ~ /Primary Assembly/) p=1; else p=0; }
        { if(p) print \$0 }' ${fasta} > primary.fa
    """
}

process normalize_fasta {
    publishDir "${params.genomes_directory}/fasta", \
        pattern: "reference.fa.gz", \
        mode: "copy", overwrite: true

    input:
    path(fasta)

    output:
    path("normalized.fa")
    path("reference.fa.gz")

    """
    picard NormalizeFasta \
    --INPUT ${fasta} \
    --OUTPUT normalized.fa \
    --LINE_LENGTH 60 \
    --USE_JDK_DEFLATER true \
    --USE_JDK_INFLATER true
    
    bgzip -c normalized.fa > reference.fa.gz
    """
}

process index_fasta {
    publishDir "${params.genomes_directory}/fasta", \
        pattern: "reference.fa.gz.{fai,gzi}", \
        mode: "copy", overwrite: true

    input:
    path(fasta)

    output:
    path("reference.fa.gz.fai")
    path("reference.fa.gz.gzi")

    """
    samtools faidx ${fasta}
    """
}

process create_sequence_dictionary {
    publishDir "${params.genomes_directory}/fasta", \
        pattern: "reference.dict", \
        mode: "copy", overwrite: true

    input:
    path(fasta)

    output:
    path("reference.dict")

    """
    picard CreateSequenceDictionary \
    --REFERENCE ${fasta} \
    --OUTPUT reference.dict \
    --SPECIES ${params.species} \
    --GENOME_ASSEMBLY ${params.assembly} \
    --USE_JDK_DEFLATER true \
    --USE_JDK_INFLATER true
     """
}

process get_reference_gtf {
    publishDir "${params.genomes_directory}/gtf", \
        pattern: "reference.gtf.gz{.url,}", \
        mode: "copy", overwrite: true

    input:
    val(gtf_url)
    path(fasta_fai)

    output:
    path("reference.gtf")
    path("reference.gtf.gz")
    path("reference.gtf.gz.url")

    """
    wget -O - "${gtf_url}" | \
    gunzip -c | \
    awk -v FS=\$'\t' -v OFS=\$'\t' '
        (NR==FNR) { a[\$1]++; next; }
        (\$0 ~ /^#/) { print \$0; next; }
        (\$1 in a) { 
            gsub(/%/, "%%");
            gsub(/gene \\"/, "gene_name \\"");
            if (\$3=="gene") { 
                gsub(/transcript_id \\"\\"/, "");
            }
            printf \$0;
            if (\$0 !~ /gene_name/) { 
                printf " gene_name \\"NA\\";"; 
            }
            if ((\$3!="gene") && (\$0 !~ /transcript_name/)) {
                printf " transcript_name \\"NA\\";";
            }
            print "";
        }
    ' ${fasta_fai} - > reference.gtf
    gzip -c reference.gtf > reference.gtf.gz
    printf "${gtf_url}" > reference.gtf.gz.url
    """
}

process reduce_gtf {
    publishDir "${params.genomes_directory}/gtf", \
        pattern: "reference.reduced.gtf.gz", \
        mode: "copy", overwrite: true

    input:
    path(gtf)
    path(dict)

    output:
    path("reference.reduced.gtf")
    path("reference.reduced.gtf.gz")

    """
    ReduceGtf \
    SEQUENCE_DICTIONARY=${dict} \
    GTF=${gtf} \
    OUTPUT=reference.reduced.gtf \
    ENHANCE_GTF=true \
    USE_JDK_DEFLATER=true \
    USE_JDK_INFLATER=true

    gzip -c reference.reduced.gtf > reference.reduced.gtf.gz
    """
}

process gtf_to_refFlat {
    publishDir "${params.genomes_directory}/gtf", \
        pattern: "reference.refFlat", \
        mode: "copy", overwrite: true

    input:
    path(gtf)
    path(dict)

    output:
    path("reference.refFlat")

    """
    ConvertToRefFlat \
    ANNOTATIONS_FILE=${gtf} \
    SEQUENCE_DICTIONARY=${dict} \
    OUTPUT=reference.refFlat
    """
}

process gtf_to_bed {
    publishDir "${params.genomes_directory}/bed", \
        pattern: "reference.bed.gz", \
        mode: "copy", overwrite: true

    input:
    path(gtf)

    output:
    path("reference.bed")
    path("reference.bed.gz")

    """
    grep -v '^#' ${gtf} | \
    awk -v FS=\$'\t' -v OFS=\$'\t' '
        {
            print \$1,\$4-1,\$5,\$1":"\$4-1"-"\$5,\$6,\$7,".",".",".",".",".",".",\$3,\$9
        }' | \
    sort -k1V -k2n -k3n > reference.bed
    gzip -c reference.bed > reference.bed.gz
    """
}

process extract_genes {
    publishDir "${params.genomes_directory}/bed", \
        pattern: "reference.genic.bed.gz", \
        mode: "copy", overwrite: true

    input:
    path(bed)

    output:
    path("reference.genic.bed")
    path("reference.genic.bed.gz")

    """
    awk '\$13 == "gene"' ${bed} | \
    cut -f1-6 | \
    bedtools merge -i stdin > reference.genic.bed
    gzip -c reference.genic.bed > reference.genic.bed.gz
    """
}

process extract_exons {
    publishDir "${params.genomes_directory}/bed", \
        pattern: "reference.exonic.bed.gz", \
        mode: "copy", overwrite: true

    input:
    path(bed)

    output:
    path("reference.exonic.bed")
    path("reference.exonic.bed.gz")

    """
    awk '\$13 == "exon"' ${bed} | \
    cut -f1-6 | \
    bedtools merge -i stdin > reference.exonic.bed
    gzip -c reference.exonic.bed > reference.exonic.bed.gz
    """
}

process extract_CDS {
    publishDir "${params.genomes_directory}/bed", \
        pattern: "reference.CDS.bed.gz", \
        mode: "copy", overwrite: true

    input:
    path(bed)

    output:
    path("reference.CDS.bed")
    path("reference.CDS.bed.gz")

    """
    awk '\$13 == "CDS"' ${bed} | \
    cut -f1-6 | \
    bedtools merge -i stdin > reference.CDS.bed
    gzip -c reference.CDS.bed > reference.CDS.bed.gz
    """
}

process extract_UTR {
    publishDir "${params.genomes_directory}/bed", \
        pattern: "reference.UTR.bed.gz", \
        mode: "copy", overwrite: true

    input:
    path(exons_bed)
    path(CDS_bed)

    output:
    path("reference.UTR.bed")
    path("reference.UTR.bed.gz")

    """
    bedtools subtract -a ${exons_bed} -b ${CDS_bed} | \
    bedtools merge -i stdin > reference.UTR.bed
    gzip -c reference.UTR.bed > reference.UTR.bed.gz
    """
}

process extract_rRNA_genes {
    publishDir "${params.genomes_directory}/bed", \
        pattern: "reference.rRNA.bed.gz", \
        mode: "copy", overwrite: true

    input:
    path(bed)

    output:
    path("reference.rRNA.bed")
    path("reference.rRNA.bed.gz")

    """
    grep -E 'gene_(bio)?type "rRNA(_pseudogene)?"' ${bed} | \
    awk '\$13 == "gene"' | \
    cut -f1-6 | \
    bedtools merge -i stdin > reference.rRNA.bed
    gzip -c reference.rRNA.bed > reference.rRNA.bed.gz
    """
}

process extract_MT_genes {
    publishDir "${params.genomes_directory}/bed", \
        pattern: "reference.MT.bed.gz", \
        mode: "copy", overwrite: true
    
    input:
    path(bed)

    output:
    path("reference.MT.bed")
    path("reference.MT.bed.gz")

    """
    sed 's/^chr//' ${bed} | \
    grep -E '^(M|MT)\t' | \
    awk '\$13 == "gene"' | \
    cut -f1-6 | \
    bedtools merge -i stdin > reference.MT.bed
    gzip -c reference.MT.bed > reference.MT.bed.gz
    """
}

process extract_intronic_regions {
    publishDir "${params.genomes_directory}/bed", \
        pattern: "reference.intronic.bed.gz", \
        mode: "copy", overwrite: true

    input:
    path(genes_bed)
    path(exons_bed)

    output:
    path("reference.intronic.bed")
    path("reference.intronic.bed.gz")

    """
    bedtools subtract -a ${genes_bed} -b ${exons_bed} | \
    bedtools merge -i stdin > reference.intronic.bed
    gzip -c reference.intronic.bed > reference.intronic.bed.gz
    """
}

process extract_intergenic_regions {
    publishDir "${params.genomes_directory}/bed", \
        pattern: "reference.intergenic.bed.gz", \
        mode: "copy", overwrite: true

    input:
    path(genes_bed)
    path(fasta_fai)

    output:
    path("reference.intergenic.bed")
    path("reference.intergenic.bed.gz")

    """
    awk -v OFS=\$'\t' '{print \$1,\$2}' ${fasta_fai} | \
    sort -k1V > chr_sizes.tsv

    bedtools complement -i ${genes_bed} -g chr_sizes.tsv | \
    bedtools merge -i stdin > reference.intergenic.bed
    gzip -c reference.intergenic.bed > reference.intergenic.bed.gz
    """
}

process bed_to_interval_list {
    publishDir "${params.genomes_directory}/interval-list", \
        pattern: "*.interval_list", \
        mode: "copy", overwrite: true

    input:
    path(dict)
    path(genes_bed)
    path(exons_bed)
    path(CDS_bed)
    path(UTR_bed)
    path(rRNA_bed)
    path(MT_bed)
    path(intronic_bed)
    path(intergenic_bed)
    
    output:
    path("reference.genes.interval_list")
    path("reference.exons.interval_list")
    path("reference.CDS.interval_list")
    path("reference.UTR.interval_list")
    path("reference.genes.rRNA.interval_list")
    path("reference.genes.MT.interval_list")
    path("reference.intronic.interval_list")
    path("reference.intergenic.interval_list")

    """
    declare -A files=(
        [genes]=${genes_bed}
        [exons]=${exons_bed}
        [CDS]=${CDS_bed}
        [UTR]=${UTR_bed}
        [genes.rRNA]=${rRNA_bed}
        [genes.MT]=${MT_bed}
        [intronic]=${intronic_bed}
        [intergenic]=${intergenic_bed}
    )
    
    for name in "\${!files[@]}"; do
        picard BedToIntervalList \
        --INPUT \${files[\$name]} \
        --OUTPUT reference.\${name}.interval_list \
        --SEQUENCE_DICTIONARY ${dict} \
        --USE_JDK_DEFLATER true \
        --USE_JDK_INFLATER true
    done
    """
}

process star_create_index {
    publishDir "${params.genomes_directory}/star-index", \
        mode: "move", overwrite: true

    input:
    path(fasta)
    path(gtf)

    output:
    path("Genome")
    path("Log.out")
    path("SA")
    path("SAindex")
    path("chrLength.txt")
    path("chrName.txt")
    path("chrNameLength.txt")
    path("chrStart.txt")
    path("exonGeTrInfo.tab")
    path("exonInfo.tab")
    path("geneInfo.tab")
    path("genomeParameters.txt")
    path("sjdbInfo.txt")
    path("sjdbList.fromGTF.out.tab")
    path("sjdbList.out.tab")
    path("transcriptInfo.tab")

    """
    gunzip -c ${fasta} > reference.fa
    gunzip -c ${gtf} > reference.gtf
    STAR \
    --runThreadN ${task.cpus} \
    --runMode genomeGenerate \
    --genomeDir . \
    --genomeFastaFiles reference.fa \
    --sjdbGTFfile reference.gtf \
    --sjdbOverhang ${params.star_sjdbOverhang} \
    --genomeSAindexNbases ${params.star_genomeSAindexNbases} \
    --limitGenomeGenerateRAM 200000000000
    """
}

process rsem_prepare_reference {
    publishDir "${params.genomes_directory}/rsem-reference",
        mode: "move", overwrite: true

    input:
    path(fasta)
    path(gtf)

    output:
    path("reference.chrlist")
    path("reference.grp")
    path("reference.idx.fa")
    path("reference.n2g.idx.fa")
    path("reference.seq")
    path("reference.ti")
    path("reference.transcripts.fa")

    """
    gunzip -c ${fasta} > reference.fa
    gunzip -c ${gtf} > reference.gtf
    rsem-prepare-reference \
    --gtf reference.gtf \
    reference.fa reference
    """
}
