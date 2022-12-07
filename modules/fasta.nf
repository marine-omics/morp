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