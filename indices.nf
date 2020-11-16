
genomes_directory = params.genomes_directory


process rsem_prepare_reference {
    publishDir "${genomes_directory}/index/rsem", \
        mode: "copy", overwrite: true, \
        pattern: "reference.chrlist"
    publishDir "${genomes_directory}/index/rsem", \
        mode: "copy", overwrite: true, \
        pattern: "reference.grp"
    publishDir "${genomes_directory}/index/rsem", \
        mode: "copy", overwrite: true, \
        pattern: "reference.idx.fa"
    publishDir "${genomes_directory}/index/rsem", \
        mode: "copy", overwrite: true, \
        pattern: "reference.n2g.idx.fa"
    publishDir "${genomes_directory}/index/rsem", \
        mode: "copy", overwrite: true, \
        pattern: "reference.seq"
    publishDir "${genomes_directory}/index/rsem", \
        mode: "copy", overwrite: true, \
        pattern: "reference.ti"
    publishDir "${genomes_directory}/index/rsem", \
        mode: "copy", overwrite: true, \
        pattern: "reference.transcripts.fa"

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

process rsem_prepare_reference_with_star {
    cpus Runtime.runtime.availableProcessors()

    publishDir "${genomes_directory}/index/rsem-with-star-overhang-${overhang}", \
        mode: "copy", overwrite: true, \
        pattern: "reference.chrlist"
    publishDir "${genomes_directory}/index/rsem-with-star-overhang-${overhang}", \
        mode: "copy", overwrite: true, \
        pattern: "reference.grp"
    publishDir "${genomes_directory}/index/rsem-with-star-overhang-${overhang}", \
        mode: "copy", overwrite: true, \
        pattern: "reference.idx.fa"
    publishDir "${genomes_directory}/index/rsem-with-star-overhang-${overhang}", \
        mode: "copy", overwrite: true, \
        pattern: "reference.n2g.idx.fa"
    publishDir "${genomes_directory}/index/rsem-with-star-overhang-${overhang}", \
        mode: "copy", overwrite: true, \
        pattern: "reference.seq"
    publishDir "${genomes_directory}/index/rsem-with-star-overhang-${overhang}", \
        mode: "copy", overwrite: true, \
        pattern: "reference.ti"
    publishDir "${genomes_directory}/index/rsem-with-star-overhang-${overhang}", \
        mode: "copy", overwrite: true, \
        pattern: "reference.transcripts.fa"
        publishDir "${genomes_directory}/index/rsem-with-star-overhang-${overhang}", \
        pattern: "Genome", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/rsem-with-star-overhang-${overhang}", \
        pattern: "Log.out", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/rsem-with-star-overhang-${overhang}", \
        pattern: "SA", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/rsem-with-star-overhang-${overhang}", \
        pattern: "SAindex", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/rsem-with-star-overhang-${overhang}", \
        pattern: "chrLength.txt", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/rsem-with-star-overhang-${overhang}", \
        pattern: "chrName.txt", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/rsem-with-star-overhang-${overhang}", \
        pattern: "chrNameLength.txt", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/rsem-with-star-overhang-${overhang}", \
        pattern: "chrStart.txt", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/rsem-with-star-overhang-${overhang}", \
        pattern: "exonGeTrInfo.tab", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/rsem-with-star-overhang-${overhang}", \
        pattern: "exonInfo.tab", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/rsem-with-star-overhang-${overhang}", \
        pattern: "geneInfo.tab", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/rsem-with-star-overhang-${overhang}", \
        pattern: "genomeParameters.txt", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/rsem-with-star-overhang-${overhang}", \
        pattern: "sjdbInfo.txt", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/rsem-with-star-overhang-${overhang}", \
        pattern: "sjdbList.fromGTF.out.tab", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/rsem-with-star-overhang-${overhang}", \
        pattern: "sjdbList.out.tab", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/rsem-with-star-overhang-${overhang}", \
        pattern: "transcriptInfo.tab", \
        mode: "copy", overwrite: true

    input:
    path(fasta)
    path(gtf)
    val(overhang)

    output:
    path("reference.chrlist")
    path("reference.grp")
    path("reference.idx.fa")
    path("reference.n2g.idx.fa")
    path("reference.seq")
    path("reference.ti")
    path("reference.transcripts.fa")
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
    rsem-prepare-reference \
        --gtf reference.gtf \
        --star -p ${task.cpus} \
        --star-sjdboverhang ${overhang} \
        reference.fa reference
    """
}


process generate_star_index {
    cpus Runtime.runtime.availableProcessors()

    publishDir "${genomes_directory}/index/star-overhang-${overhang}", \
        pattern: "Genome", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/star-overhang-${overhang}", \
        pattern: "Log.out", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/star-overhang-${overhang}", \
        pattern: "SA", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/star-overhang-${overhang}", \
        pattern: "SAindex", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/star-overhang-${overhang}", \
        pattern: "chrLength.txt", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/star-overhang-${overhang}", \
        pattern: "chrName.txt", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/star-overhang-${overhang}", \
        pattern: "chrNameLength.txt", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/star-overhang-${overhang}", \
        pattern: "chrStart.txt", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/star-overhang-${overhang}", \
        pattern: "exonGeTrInfo.tab", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/star-overhang-${overhang}", \
        pattern: "exonInfo.tab", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/star-overhang-${overhang}", \
        pattern: "geneInfo.tab", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/star-overhang-${overhang}", \
        pattern: "genomeParameters.txt", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/star-overhang-${overhang}", \
        pattern: "sjdbInfo.txt", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/star-overhang-${overhang}", \
        pattern: "sjdbList.fromGTF.out.tab", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/star-overhang-${overhang}", \
        pattern: "sjdbList.out.tab", \
        mode: "copy", overwrite: true
    publishDir "${genomes_directory}/index/star-overhang-${overhang}", \
        pattern: "transcriptInfo.tab", \
        mode: "copy", overwrite: true

    input:
    path(fasta)
    path(gtf)
    val(overhang)

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
    --sjdbOverhang ${overhang} \
    --genomeSAindexNbases ${params.star_genomeSAindexNbases} \
    --limitGenomeGenerateRAM 200000000000
    """
}