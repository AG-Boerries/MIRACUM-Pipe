stats <- function(path, outfile_pdf, stats_td, stats_gd = NULL,
                  protocol = "Tumor_Normal"){
  #' Statistics
  #'
  #' @description Information about Coverage and Readdepth
  #'
  #' @param path string. Path to statistics data
  #' @param outfile_pdf string. Filename for coverage plot
  #' @param stats_td numerical. Total reads for tumor
  #' @param stats_gd numerical. Total reads for germline
  #' @param protocol string. Protocol of analysis
  #' 
  #' @return list of
  #' @return cover vector. Mean coverage 
  #' @return avreads vector. Total reads for germline and tumor
  #'
  #' @details Statistical numbers are extracted from alignment's statistics.
  cover <- coverage_plot(path = path, outfilePDF = coverage_out, protocol = protocol)
  if (protocol == "Tumor_Normal"){
    avreads <- reads(stats_td, stats_gd)
  } else {
    avreads <- treads(stats_td)
  }
  return(list(cover = cover, avreads = avreads))
}