nextflow.enable.dsl = 2

include { ART_ILLUMINA } from 'modules/nf-core/art_illumina'
include { RNASEQ } from 'modules/nf-core/rnaseq'

params {
    genome_fasta = './genome/methylorubrum.fna' // Reference genome
    genome_gtf = './genome/methylorubrum.gff' // Annotation file
    read_length = 150  // Read length for simulated reads
    coverage = 10  // Coverage depth
    output_dir = './simulated_reads/'  // Output directory for simulated reads
    rnaseq_outdir = './rnaseq_results/' // Output directory for RNA-seq analysis
}

process SimulateReads {
    tag "Simulating reads"
    
    input:
    path genome_fasta

    output:
    path "${params.output_dir}methylorubrum_R1.fq", emit: read1
    path "${params.output_dir}methylorubrum_R2.fq", emit: read2

    script:
    """
    mkdir -p ${params.output_dir}
    art_illumina -ss HS25 -i ${genome_fasta} \\
      -l ${params.read_length} -f ${params.coverage} -m 200 -s 10 \\
      -o ${params.output_dir}methylorubrum

    mv ${params.output_dir}methylorubrum1.fq ${params.output_dir}methylorubrum_R1.fq
    mv ${params.output_dir}methylorubrum2.fq ${params.output_dir}methylorubrum_R2.fq
    """
}

workflow {
    simulated_reads = SimulateReads(params.genome_fasta)

    rnaseq_input = simulated_reads.collect().map { files -> tuple('SimulatedSample', files[0], files[1]) }

    rnaseq_input | RNASEQ \
        --fasta params.genome_fasta \
        --gtf params.genome_gtf \
        --outdir params.rnaseq_outdir
}
