process DVPOLISH_PBMM2_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    // Note: the versions here need to match the versions used in pbmm2/index
    conda 'bioconda::pbmm2=1.13.1'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbmm2:1.13.1--h9ee0642_0' :
        'biocontainers/pbmm2:1.13.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("*.bam", arity: '1..*'), optional: true, emit: bam
    tuple val(meta), path("*.bai", arity: '1..*'), optional: true, emit: bai
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''

    def out_name_part1 = reference.name.endsWith(".gz") ? reference.getBaseName(2) : reference.baseName
    def out_name_part2 = reads.name.endsWith(".gz") ?  reads.getBaseName(2) : reads.baseName

    """
    pbmm2 align --sort \\
        $args \\
        -j $task.cpus \\
        "$reference" \\
        "$reads" \\
        ${out_name_part1}_${out_name_part2}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(pbmm2 --version 2>&1 | head -n 1)
    END_VERSIONS
    """
}
