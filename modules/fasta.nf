process combine_refs {
  input:
    path(refa)
    path(refb)

  output:
    path("combined.fasta"), emit:ref

  script:
  """
  cat ${refa} | sed 's/>/>a_/' > tmp.a
  cat ${refb} | sed 's/>/>b_/' > tmp.b
  cat tmp.a tmp.b > combined.fasta
  """
}

process combine_maps {
  input:
    path(refa_map)
    path(refb_map)

  output:
    path("combined.gene2trans.map"), emit:map

  script:
  """
  cat ${refa_map} | awk -F '\t' '{ln=sprintf("a_%s:a_%s",\$1,\$2);print ln}' | sed 's/:/\t/'> tmp.a
  cat ${refb_map} | awk -F '\t' '{ln=sprintf("b_%s:b_%s",\$1,\$2);print ln}' | sed 's/:/\t/'> tmp.b
  cat tmp.a tmp.b > combined.gene2trans.map  
  """

}

