process rsem_prepare_reference {
    input:
    path fasta

    output:
    path "rsem_ref*" , emit: ref

    script:
    """
    rsem-prepare-reference --bowtie2 ${fasta} rsem_ref
    """
}

process rsem_calculate_expression {

  publishDir "$params.outdir/rsem", mode: 'copy'

  input:
  tuple val(meta), path(reads)
  path(ref)

  output:
  tuple val(meta), path("*.results"), emit: gene_expr
  tuple val(meta), path("*.bam"), emit: bams  
  tuple val(meta), path("*.stat"), emit: stats

  script:

  def args = task.ext.args ?: ''

  if (meta.single_end) {
  """
    rsem-calculate-expression \\
        --bowtie2 \\
        --num-threads ${task.cpus} \\
         ${reads.join(',')} \\
         rsem_ref ${meta.sample}
  """
    } else {
   read1=[]
   read2=[]
   reads.eachWithIndex{ v, ix -> ( ix & 1 ? read2 : read1 ) << v }

  """
    rsem-calculate-expression \\
        --paired-end --bowtie2 \\
        --num-threads ${task.cpus} \\
         ${read1.join(',')} ${read2.join(',')} \\
         rsem_ref ${meta.sample}
  """
    }

}


