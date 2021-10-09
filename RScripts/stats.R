stats <- function(path = path2coverage, outfile_pdf = coverage_out, stats_td,
      stats_gd, path_output, protocol){
  #' Statistics
  #'
  #' @description Information about Coverage and Readdepth
  #'
  #' @param path string. Path to statistics data
  #' @param outfile_pdf string. Filename for coverage plot
  #' @param stats_td numerical. Total reads for tumor
  #' @param stats_gd numerical. Total reads for germline
  #' @param path_output string. Output directory
  #' @param protocol strin. Protocol of analyses
  #' 
  #' @return list of
  #' @return cover vector. Mean coverage 
  #' @return avreads vector. Total reads for germline and tumor
  #'
  #' @details Statistical numbers are extracted from alignment's statistics.
  cover <- coverage_plot(path = path, outfilePDF = coverage_out, protocol = protocol)
  cover_exons <- coverage_exon(path = path, protocol = protocol)
  if (protocol == "somaticGermline" | protocol == "somatic"){
    avreads <- reads(stats_td, stats_gd)
  } else {
    avreads <- treads(stats_td)
  }
  qc_check <- quality_check(path = path, nsamples = cover$labs, protocol = protocol)
  
  return(list(cover = cover, cover_exons = cover_exons, avreads = avreads, qc_check = qc_check))
}
