
genomes_directory = params.genomes_directory
basename = params.basename
picard_jar = "/working/software/picard.jar"

process get_reference_fasta {
    publishDir "${genomes_directory}/fasta", \
        pattern: "reference.fa.gz.url", \
        mode: "copy", overwrite: true, \
        saveAs: { filename -> "${basename}.fa.gz.url" }

    input:
    val(fasta_url)

    output:
    path("reference.fa.gz")
    path("reference.fa.gz.url")

    """
    wget -O - ${fasta_url} | \
    gunzip -c | \
    awk '
        BEGIN { p=0; }
        \$0 ~ /^>/ { if(\$0 ~ /(Primary Assembly)|(primary_assembly)/) p=1; else p=0; }
        { if(p) print \$0 }
    ' | \
    gzip -c > reference.fa.gz
    printf "${fasta_url}\n" > reference.fa.gz.url
    """
}

process extract_primary_assembly {
    input:
    path(fasta)

    output:
    path("primary.fa.gz")

    """
    gunzip -c ${fasta} | \
    awk '
        BEGIN { p=0; }
        \$0 ~ /^>/ { if(\$0 ~ /(Primary Assembly)|(primary_assembly)/) p=1; else p=0; }
        { if(p) print \$0 }' | \
    gzip -c > primary.fa.gz
    """
}

process normalize_fasta {
    publishDir "${genomes_directory}/fasta", \
        pattern: "normalized.fa.gz", \
        mode: "copy", overwrite: true, \
        saveAs: { filename -> "${basename}.fa.gz" }

    input:
    path(fasta)

    output:
    path("normalized.fa")
    path("normalized.fa.gz")

    """
    java -jar ${picard_jar} \
        NormalizeFasta \
        --INPUT ${fasta} \
        --OUTPUT normalized.fa \
        --LINE_LENGTH 60 \
        --USE_JDK_DEFLATER true \
        --USE_JDK_INFLATER true
    bgzip -c normalized.fa > normalized.fa.gz
    """
}

process index_fasta {
    publishDir "${genomes_directory}/fasta", \
        pattern: "normalized.fa.gz.fai", \
        mode: "copy", overwrite: true, \
        saveAs: { filename -> "${basename}.fa.gz.fai" }
    publishDir "${genomes_directory}/fasta", \
        pattern: "normalized.fa.gz.gzi", \
        mode: "copy", overwrite: true, \
        saveAs: { filename -> "${basename}.fa.gz.gzi" }

    input:
    path(fasta)

    output:
    path("normalized.fa.gz.fai")
    path("normalized.fa.gz.gzi")

    """
    samtools faidx \
    --fai-idx normalized.fa.gz.fai  \
    --gzi-idx normalized.fa.gz.gzi \
    ${fasta}
    """
}

process create_sequence_dictionary {
    publishDir "${genomes_directory}/fasta", \
        pattern: "reference.dict", \
        mode: "copy", overwrite: true, \
        saveAs: { filename -> "${basename}.dict" }

    input:
    path(fasta)

    output:
    path("reference.dict")

    """
    java -jar ${picard_jar} \
        CreateSequenceDictionary \
        --REFERENCE ${fasta} \
        --OUTPUT reference.dict \
        --SPECIES ${params.species} \
        --GENOME_ASSEMBLY ${params.assembly} \
        --USE_JDK_DEFLATER true \
        --USE_JDK_INFLATER true
     """
}

process get_reference_gtf {
    publishDir "${genomes_directory}/gtf", \
        pattern: "reference.gtf.gz", \
        mode: "copy", overwrite: true, \
        saveAs: { filename -> "${basename}.gtf.gz" }
    publishDir "${genomes_directory}/gtf", \
        pattern: "reference.gtf.gz.url", \
        mode: "copy", overwrite: true, \
        saveAs: { filename -> "${basename}.gtf.gz.url" }

    input:
    val(gtf_url)

    output:
    path("reference.gtf")
    path("reference.gtf.gz")
    path("reference.gtf.gz.url")

    """
    wget -O reference.gtf.gz "${gtf_url}"
    echo "${gtf_url}" > reference.gtf.gz.url
    gunzip -c reference.gtf.gz > reference.gtf
    """
}

process reduce_gtf {
    publishDir "${genomes_directory}/gtf", \
        pattern: "reduced.gtf.gz", \
        mode: "copy", overwrite: true, \
        saveAs: { filename -> "${basename}.reduced.gtf.gz" }

    input:
    path(gtf)
    path(dict)

    output:
    path("reduced.gtf")
    path("reduced.gtf.gz")

    """
    ReduceGtf \
    SEQUENCE_DICTIONARY=${dict} \
    GTF=${gtf} \
    OUTPUT=reduced.gtf \
    ENHANCE_GTF=true \
    USE_JDK_DEFLATER=true \
    USE_JDK_INFLATER=true

    gzip -c reduced.gtf > reduced.gtf.gz
    """
}

process gtf_to_refFlat {
    publishDir "${genomes_directory}/gtf", \
        pattern: "reference.refFlat", \
        mode: "copy", overwrite: true, \
        saveAs: { filename -> "${basename}.refFlat" }

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
    publishDir "${genomes_directory}/bed", \
        pattern: "reference.bed.gz", \
        mode: "copy", overwrite: true, \
        saveAs: { filename -> "${basename}.bed.gz" }

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
    publishDir "${genomes_directory}/bed", \
        pattern: "genes.bed.gz", \
        mode: "copy", overwrite: true, \
        saveAs: { filename -> "${basename}.genes.bed.gz" }

    input:
    path(bed)

    output:
    path("genes.bed")
    path("genes.bed.gz")

    """
    awk '\$13 == "gene"' ${bed} > genes.bed
    gzip -c genes.bed > genes.bed.gz
    """
}

process extract_exons {
    publishDir "${genomes_directory}/bed", \
        pattern: "exons.bed.gz", \
        mode: "copy", overwrite: true, \
        saveAs: { filename -> "${basename}.exons.bed.gz" }

    input:
    path(bed)

    output:
    path("exons.bed")
    path("exons.bed.gz")

    """
    awk '\$13 == "exon"' ${bed} > exons.bed
    gzip -c exons.bed > exons.bed.gz
    """
}

process extract_CDS {
    publishDir "${genomes_directory}/bed", \
        pattern: "CDS.bed.gz", \
        mode: "copy", overwrite: true, \
        saveAs: { filename -> "${basename}.CDS.bed.gz" }

    input:
    path(bed)

    output:
    path("CDS.bed")
    path("CDS.bed.gz")

    """
    awk '\$13 == "CDS"' ${bed} > CDS.bed
    gzip -c CDS.bed > CDS.bed.gz
    """
}

process extract_rRNA_genes {
    publishDir "${genomes_directory}/bed", \
        pattern: "genes.rRNA.bed.gz", \
        mode: "copy", overwrite: true, \
        saveAs: { filename -> "${basename}.genes.rRNA.bed.gz" }

    input:
    path(bed)

    output:
    path("genes.rRNA.bed")
    path("genes.rRNA.bed.gz")

    """
    grep -E 'gene_(bio)?type "rRNA(_pseudogene)?"' ${bed} | \
    awk '\$13 == "gene"' > genes.rRNA.bed
    gzip -c genes.rRNA.bed > genes.rRNA.bed.gz
    """
}

process extract_MT_genes {
    publishDir "${genomes_directory}/bed", \
        pattern: "genes.MT.bed.gz", \
        mode: "copy", overwrite: true, \
        saveAs: { filename -> "${basename}.genes.MT.bed.gz" }
    
    input:
    path(bed)

    output:
    path("genes.MT.bed")
    path("genes.MT.bed.gz")

    """
    sed 's/^chr//' ${bed} | \
    grep -E '^(M|MT)\t' | \
    awk '\$13 == "gene"' > genes.MT.bed
    gzip -c genes.MT.bed > genes.MT.bed.gz
    """
}

process extract_intronic_regions {
    publishDir "${genomes_directory}/bed", \
        pattern: "intronic.bed.gz", \
        mode: "copy", overwrite: true, \
        saveAs: { filename -> "${basename}.intronic.bed.gz" }

    input:
    path(genes_bed)
    path(exons_bed)

    output:
    path("intronic.bed")
    path("intronic.bed.gz")

    """
    bedtools merge -i ${exons_bed} -s | \
    bedtools subtract -a ${genes_bed} -b stdin -s > intronic.bed
    gzip -c intronic.bed > intronic.bed.gz
    """
}

process extract_intergenic_regions {
    publishDir "${genomes_directory}/bed", \
        pattern: "intergenic.bed.gz", \
        mode: "copy", overwrite: true, \
        saveAs: { filename -> "${basename}.intergenic.bed.gz" }

    input:
    path(genes_bed)
    path(fasta_fai)

    output:
    path("intergenic.bed")
    path("intergenic.bed.gz")

    """
    awk -v OFS=\$'\t' '{print \$1,\$2}' ${fasta_fai} | \
    sort -k1n -k2n > chr_sizes.tsv

    bedtools complement -i ${genes_bed} -g chr_sizes.tsv > intergenic.bed
    gzip -c intergenic.bed > intergenic.bed.gz
    """
}

process bed_to_interval_list {
    publishDir "${genomes_directory}/interval-list", \
        pattern: "*.interval_list", \
        mode: "copy", overwrite: true, \
        saveAs: { filename -> "${basename}.${filename}" }

    input:
    path(dict)
    path(genes_bed)
    path(exons_bed)
    path(CDS_bed)
    path(rRNA_bed)
    path(MT_bed)
    path(intronic_bed)
    path(intergenic_bed)
    
    output:
    path("genes.interval_list")
    path("exons.interval_list")
    path("CDS.interval_list")
    path("genes.rRNA.interval_list")
    path("genes.MT.interval_list")
    path("intronic.interval_list")
    path("intergenic.interval_list")

    """
    java -jar ${picard_jar} \
        BedToIntervalList \
        --INPUT ${genes_bed} \
        --OUTPUT genes.interval_list \
        --SEQUENCE_DICTIONARY ${dict} \
        --USE_JDK_DEFLATER true \
        --USE_JDK_INFLATER true

    java -jar ${picard_jar} \
        BedToIntervalList \
        --INPUT ${exons_bed} \
        --OUTPUT exons.interval_list \
        --SEQUENCE_DICTIONARY ${dict} \
        --USE_JDK_DEFLATER true \
        --USE_JDK_INFLATER true

    java -jar ${picard_jar} \
        BedToIntervalList \
        --INPUT ${CDS_bed} \
        --OUTPUT CDS.interval_list \
        --SEQUENCE_DICTIONARY ${dict} \
        --USE_JDK_DEFLATER true \
        --USE_JDK_INFLATER true

    java -jar ${picard_jar} \
        BedToIntervalList \
        --INPUT ${rRNA_bed} \
        --OUTPUT genes.rRNA.interval_list \
        --SEQUENCE_DICTIONARY ${dict} \
        --USE_JDK_DEFLATER true \
        --USE_JDK_INFLATER true

    java -jar ${picard_jar} \
        BedToIntervalList \
        --INPUT ${MT_bed} \
        --OUTPUT genes.MT.interval_list \
        --SEQUENCE_DICTIONARY ${dict} \
        --USE_JDK_DEFLATER true \
        --USE_JDK_INFLATER true

    java -jar ${picard_jar} \
        BedToIntervalList \
        --INPUT ${intronic_bed} \
        --OUTPUT intronic.interval_list \
        --SEQUENCE_DICTIONARY ${dict} \
        --USE_JDK_DEFLATER true \
        --USE_JDK_INFLATER true

    java -jar ${picard_jar} \
        BedToIntervalList \
        --INPUT ${intergenic_bed} \
        --OUTPUT intergenic.interval_list \
        --SEQUENCE_DICTIONARY ${dict} \
        --USE_JDK_DEFLATER true \
        --USE_JDK_INFLATER true
    """
}

process generate_star_index {
    cpus Runtime.runtime.availableProcessors()

    publishDir "${genomes_directory}/star", \
        pattern: "reference.star_idx.cmd", \
        mode: "copy", overwrite: true, \
        saveAs: { filename -> "${basename}.star_idx.cmd" }
    publishDir "${genomes_directory}/star/${basename}.${read_length}_bp.star_idx", \
        pattern: "reference.star_idx/Genome", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/star/${basename}.${read_length}_bp.star_idx", \
        pattern: "reference.star_idx/Log.out", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/star/${basename}.${read_length}_bp.star_idx", \
        pattern: "reference.star_idx/SA", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/star/${basename}.${read_length}_bp.star_idx", \
        pattern: "reference.star_idx/SAindex", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/star/${basename}.${read_length}_bp.star_idx", \
        pattern: "reference.star_idx/chrLength.txt", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/star/${basename}.${read_length}_bp.star_idx", \
        pattern: "reference.star_idx/chrName.txt", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/star/${basename}.${read_length}_bp.star_idx", \
        pattern: "reference.star_idx/chrNameLength.txt", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/star/${basename}.${read_length}_bp.star_idx", \
        pattern: "reference.star_idx/chrStart.txt", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/star/${basename}.${read_length}_bp.star_idx", \
        pattern: "reference.star_idx/exonGeTrInfo.tab", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/star/${basename}.${read_length}_bp.star_idx", \
        pattern: "reference.star_idx/exonInfo.tab", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/star/${basename}.${read_length}_bp.star_idx", \
        pattern: "reference.star_idx/geneInfo.tab", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/star/${basename}.${read_length}_bp.star_idx", \
        pattern: "reference.star_idx/genomeParameters.txt", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/star/${basename}.${read_length}_bp.star_idx", \
        pattern: "reference.star_idx/sjdbInfo.txt", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/star/${basename}.${read_length}_bp.star_idx", \
        pattern: "reference.star_idx/sjdbList.fromGTF.out.tab", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/star/${basename}.${read_length}_bp.star_idx", \
        pattern: "reference.star_idx/sjdbList.out.tab", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/star/${basename}.${read_length}_bp.star_idx", \
        pattern: "reference.star_idx/transcriptInfo.tab", \
        mode: "copy", overwrite: true

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
    echo \$cmd >> reference.star_idx.cmd
    """
}

process rsem_prepare_reference {
    publishDir "${genomes_directory}/rsem", \
        mode: "copy", overwrite: true, \
        pattern: "rsem.chrlist",
        saveAs: { filename -> "${basename}.chrlist" }
    publishDir "${genomes_directory}/rsem", \
        mode: "copy", overwrite: true, \
        pattern: "rsem.grp",
        saveAs: { filename -> "${basename}.grp" }
    publishDir "${genomes_directory}/rsem", \
        mode: "copy", overwrite: true, \
        pattern: "rsem.idx.fa",
        saveAs: { filename -> "${basename}.idx.fa" }
    publishDir "${genomes_directory}/rsem", \
        mode: "copy", overwrite: true, \
        pattern: "rsem.n2g.idx.fa",
        saveAs: { filename -> "${basename}.n2g.idx.fa" }
    publishDir "${genomes_directory}/rsem", \
        mode: "copy", overwrite: true, \
        pattern: "rsem.seq",
        saveAs: { filename -> "${basename}.seq" }
    publishDir "${genomes_directory}/rsem", \
        mode: "copy", overwrite: true, \
        pattern: "rsem.ti",
        saveAs: { filename -> "${basename}.ti" }
    publishDir "${genomes_directory}/rsem", \
        mode: "copy", overwrite: true, \
        pattern: "rsem.transcripts.fa",
        saveAs: { filename -> "${basename}.transcripts.fa" }

    input:
    path(fasta)
    path(gtf)

    output:
    path("rsem.chrlist")
    path("rsem.grp")
    path("rsem.idx.fa")
    path("rsem.n2g.idx.fa")
    path("rsem.seq")
    path("rsem.ti")
    path("rsem.transcripts.fa")

    """
    gunzip -c ${fasta} > reference.fa
    gunzip -c ${gtf} > reference.gtf
    rsem-prepare-reference --gtf reference.gtf reference.fasta rsem
    """
}


