strip_mdpi = function(in_rmd = "doc/peakcharacterization_manuscript.Rmd", out_rmd = "doc/peakcharacterization_nostyle.Rmd"){
  rmd_doc = readLines(in_rmd)

  template_line = grepl("metabolites-template", rmd_doc)
  rmd_notemplate = rmd_doc[!template_line]

  has_style = grepl(":::", rmd_notemplate)
  rmd_nostyle = rmd_notemplate[!has_style]

  cat(rmd_nostyle, file = out_rmd, sep = "\n")
  out_rmd
}

strip_headers = function(in_rmd = "doc/peakcharacterization_manuscript.Rmd",
                         out_rmd = "doc/peakcharacterization_mdpi.Rmd"){

  rmd_doc = readLines(in_rmd)

  has_header = grepl("^##+ \\w+", rmd_doc)
  doc_out = rmd_doc
  doc_out[has_header] = gsub("^##+ ", "", rmd_doc[has_header])
  cat(doc_out, file = out_rmd, sep = "\n")
  out_rmd
}
