nextflow.enable.dsl=2

include { multiqc_fastp;multiqc_rsem } from './modules/multiqc.nf'
include { fastp } from './modules/fastp.nf'
include { combine_refs;combine_maps } from './modules/fasta.nf'
include { rsem_prepare_reference;rsem_prepare_reference_t2gmap;rsem_calculate_expression } from './modules/rsem.nf'

params.refa=null
params.refb=null
params.refa_map=null
params.refb_map=null
params.outdir=null

if(!params.outdir){
  log.error "No outdir provided. Provide one with --outdir myoutdir"
  exit 1
}

workflow preprocess {
  take:
    fastqin

  main:
    fastp(fastqin) 
    fastp.out.json | collect | multiqc_fastp

  emit:
    fastp.out.reads
}

workflow {
  refa_fasta = Channel.fromPath(file(params.refa, checkIfExists:true)) | collect
  
  if ( params.refa_map ){
    refa_map = Channel.fromPath(file(params.refa_map,checkIfExists:true)) | collect
  }

  if ( params.refb ){
    refb_fasta = Channel.fromPath(file(params.refb, checkIfExists:true)) | collect
    ref_fasta = combine_refs(refa_fasta,refb_fasta) | collect

    if ( params.refb_map ){
      refb_map = Channel.fromPath(file(params.refb_map,checkIfExists:true)) | collect
    }

  } else {
    ref_fasta = refa_fasta
  }

  if ( params.refa_map && params.refb_map ){
    println "Will combine maps"
    ref_map = combine_maps(refa_map,refb_map) | collect
    rsem_ref = rsem_prepare_reference_t2gmap(ref_fasta,ref_map) | collect
  } else {
    if ( params.refb && ((params.refa_map && !params.refb_map) || (!params.refa_map && params.refb_map)) ){
      log.error "Two refs provided but only one transcript2gene map. If you provide one map you must provide both"
      exit 1
    }
    if ( params.refa_map && !params.refa_map){
      println "No transcript2gene map provided. Proceeding without map"
      rsem_ref = rsem_prepare_reference(ref_fasta) | collect
    } else {
        println "Found transcript2gene map. Will attempt to use this map"
        rsem_ref = rsem_prepare_reference_t2gmap(ref_fasta,refa_map) | collect
    }
  }
// Preprocess data
  ch_input_sample = extract_csv(file(params.samples, checkIfExists: true))

  ch_trimmed_samples = ch_input_sample | preprocess

  ch_grouped_samples = ch_trimmed_samples.map { meta,reads ->
    meta = [sample:meta.sample,single_end:meta.single_end]
    tuple(meta,reads)
  }.groupTuple().map { meta,pairList ->
      flatreads=[]
      if ( !meta.single_end ){
        pairList.each { v -> 
          flatreads << v[0]
          flatreads << v[1]
        }
      } else {
        flatreads << pairList[0]
      }
    tuple(meta,flatreads)
  }

  rsem_calculate_expression(ch_grouped_samples,rsem_ref)

  rsem_calculate_expression.out.stats.map{ meta,stats -> stats } | collect | multiqc_rsem

}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def longestCommonPrefix(a,b){
  int minLength = Math.min(a.length(),b.length())

  prefix=a.substring(0, minLength);
  for (int i = 0; i < minLength; i++) {
        if (a.charAt(i) != b.charAt(i)) {
            prefix = a.substring(0, i);
        }
    }
  if ( prefix.endsWith(".")){
    prefix=prefix.substring(0,prefix.length()-1)
  }
  return prefix
}

def resolve_path(pathstring){
  if(pathstring =~ /^\//){
    pathstring
  } else {
    "${params.base_path}/${pathstring}"
  }
}

def extract_csv(csv_file) {
    Channel.from(csv_file).splitCsv(header: true)
    .map{ row -> 
      def meta = [:]
      meta.sample = row.sample

      def fastq_1     = file(resolve_path(row.fastq_1), checkIfExists: true)
      def fastq_2     = row.fastq_2 ? file(resolve_path(row.fastq_2), checkIfExists: true) : null
      meta.single_end = row.fastq_2 ? false : true

      if (!meta.single_end){
        prefix=longestCommonPrefix(fastq_1.getName(),fastq_2.getName())
        meta.prefix = prefix
        reads = [fastq_1,fastq_2]
      } else {
        reads = [fastq_1]
        meta.prefix = fastq_1.getBaseName()         
      }

      reads.removeAll([null])

      [meta,reads]
    }
}
