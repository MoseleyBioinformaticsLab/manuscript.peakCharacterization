strip_mdpi_render = function(in_rmd = "doc/peakcharacterization_manuscript.Rmd", out_rmd = "doc/peakcharacterization_nostyle.Rmd", out_file = here::here("doc/peakcharacterization_nostyle.docx"), output_format = "word_document"){
  rmd_doc = readLines(in_rmd)

  template_line = grepl("metabolites-template", rmd_doc)
  rmd_notemplate = rmd_doc[!template_line]

  has_style = grepl(":::", rmd_notemplate)
  rmd_nostyle = rmd_notemplate[!has_style]

  cat(rmd_nostyle, file = out_rmd, sep = "\n")
  rmarkdown::render(out_rmd,
                    output_format = output_format,
                    output_file = out_file,
                    knit_root_dir = getwd(),
                    quiet = TRUE)
  beepr::beep(4)
  out_rmd
}

strip_headers_render = function(in_rmd = "doc/peakcharacterization_manuscript.Rmd",
                         out_rmd = "doc/peakcharacterization_mdpi.Rmd"){

  rmd_doc = readLines(in_rmd)

  has_header = grepl("^##+ \\w+", rmd_doc)
  has_bullet = grepl("^\\* ", rmd_doc)
  doc_out = rmd_doc
  doc_out[has_header] = gsub("^##+ ", "", rmd_doc[has_header])
  doc_out[has_bullet] = gsub("^\\* ", "", rmd_doc[has_bullet])
  cat(doc_out, file = out_rmd, sep = "\n")
  rmarkdown::render(out_rmd, knit_root_dir = getwd(), quiet = TRUE)
  beepr::beep(4)
  out_rmd
}
