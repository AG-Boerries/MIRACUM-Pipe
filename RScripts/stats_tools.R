#################
# Coverage Plot #
coverage_plot <- function(path, outfilePDF, protocol){
  #' Coverage Plot
  #'
  #' @description Coverage Plot
  #'
  #' @param path string. Path to data
  #' @param outfilePDF string. Name of output file
  #' @param protocol strin. Protocl of analyses
  #'
  #' @return list of
  #' @return cov numerical. Mean coverage
  #' @return labs vector of strings. IDs for coverage files
  #' @return files vector of strings. Files that contain coverage
  #' 
  #' 
  #' @details A plot is generated to describe the coverage over the sequenced
  #' @details sample.
  require(RColorBrewer)
  # Assumes you've already run coverageBed -hist, and grep'd '^all'. E.g. something like:
  # find *.bam | parallel 'bedtools -abam {} -b capture.bed -hist | grep ^all > {}.all.txt'
  # Get a list of the bedtools output files you'd like to read in
  print(files <- list.files(path = path, pattern="all.txt$"))
  # Optional, create short sample names from the filenames. 
  # For example, in this experiment, my sample filenames might look like this:
  # prefixToTrash-01.pe.on.pos.dedup.realigned.recalibrated.bam
  # prefixToTrash-02.pe.on.pos.dedup.realigned.recalibrated.bam
  # prefixToTrash-03.pe.on.pos.dedup.realigned.recalibrated.bam
  # This regular expression leaves me with "samp01", "samp02", and "samp03" in the legend.
  
  print(labs <- gsub("_coverage.all.txt", "", files, perl = TRUE))
  files <- paste(path, files, sep = '')
  
  # Create lists to hold coverage and cumulative coverage for each alignment,
  # and read the data into these lists.
  cov <- list()
  cov_cumul <- list()
  for (i in 1:length(files)) {
    cov[[i]] <- read.table(files[i])
    cov_cumul[[i]] <- 1 - cumsum(cov[[i]][, 5])
  }
  
  # Pick some colors
  # Ugly:
  # cols <- 1:length(cov)
  # Prettier:
  # ?colorRampPalette
  # display.brewer.all()
  l <- length(cov)
  if(l <= 2){l <- 3}
  cols <- brewer.pal(l, "Dark2")
  
  # Save the graph to a file
  pdf(file = outfilePDF, height = 12, width = 12, pointsize=20)
  
  # Create plot area, but do not plot anything. Add gridlines and axis labels.
 if (protocol != "panelTumor"){
    plot(cov[[1]][2:401, 2], cov_cumul[[1]][1:400], type='n', xlab="Depth",
       ylab = expression("Fraction of capture target bases " >= " depth"),
       ylim = c(0,1.0), main = "Target Region Coverage")
    abline(v = 20, col = "gray60")
    abline(v = 50, col = "gray60")
    abline(v = 80, col = "gray60")
    abline(v = 100, col = "gray60")
    abline(h = 0.50, col = "gray60")
    abline(h = 0.90, col = "gray60")
    axis(1, at = c(20, 50, 80), labels = c(20, 50, 80))
    axis(2, at = c(0.90), labels = c(0.90))
    axis(2, at = c(0.50), labels = c(0.50))
  } else {
    plot(cov[[1]][2:10001, 2], cov_cumul[[1]][1:10000], type='n', xlab="Depth",
         ylab = expression("Fraction of capture target bases " >= " depth"),
         ylim = c(0,1.0), main = "Target Region Coverage")
    abline(v = 2000, col = "gray60")
    abline(v = 5000, col = "gray60")
    abline(v = 8000, col = "gray60")
    abline(v = 10000, col = "gray60")
    abline(h = 0.50, col = "gray60")
    abline(h = 0.90, col = "gray60")
    axis(1, at = c(2000, 5000, 8000), labels = c(2000, 5000, 8000))
    axis(2, at = c(0.90), labels = c(0.90))
    axis(2, at = c(0.50), labels = c(0.50))
  }
  
  # Actually plot the data for each of the alignments (stored in the lists).
if (protocol != "panelTumor"){
    for (i in 1:length(cov)) points(cov[[i]][2:401, 2], cov_cumul[[i]][1:400],
                                  type = "l", lwd = 3, col = cols[i])
  } else {
    for (i in 1:length(cov)) points(cov[[i]][2:10001, 2], cov_cumul[[i]][1:10000],
                                    type = "l", lwd = 3, col = cols[i])
  }
  
  # Add a legend using the nice sample labeles rather than the full filenames.
  legend("topright", legend = labs, col = cols, lty = 1, lwd = 4)
  
  dev.off()
  #return(cov)
  
  
  # Mean Coverage
  for (i in 1:length(files)) {
    print(paste('Mean Coverage', labs[i], ':',
                sum(cov[[i]][,2] * cov[[i]][, 5]), sep = " "))
  }
  return(list(cov = cov,labs = labs, files = files))
}

reads <- function(tfile, gfile){
  #' Reads
  #'
  #' @description Extract the mean readcount
  #'
  #' @param tfile string. Name of germline input file
  #' @param gfile string. Name of germline input file
  #'
  #' @return list of
  #' @return nRT numerical. Total number of reads in tumor data 
  #' @return nRG numerical. Total number of reads in germline data
  #' 
  #' 
  #' @details The statistics files contain information about properly paired
  #' @details reads. The information is extracted for the report.
  treads <- read.table(file = tfile, sep = "\t", skip = 7, nrows = 31, fill = TRUE)
  id <- which (as.character(treads$V2) == "reads properly paired:")
  treads <- as.character(treads$V3[id])
  treads <- as.numeric(treads)/1000000
  ntreads <- round(treads)
  
  greads <- read.table(file = gfile, sep = "\t", skip = 7, nrows = 31, fill = TRUE)
  id <- which (as.character(greads$V2) == "reads properly paired:")
  greads <- as.character(greads$V3[id])
  greads <- as.numeric(greads)/1000000
  ngreads <- round(greads)
  
  return(list(nRT = ntreads, nRG = ngreads))
}

treads <- function(tfile){
  treads <- read.table(file = tfile, sep = "\t", skip = 7, nrows = 31, fill = TRUE)
  id <- which (as.character(treads$V2) == "reads properly paired:")
  treads <- as.character(treads$V3[id])
  treads <- as.numeric(treads)/1000000
  ntreads <- round(treads)

return(list(nRT = ntreads, nRG = NULL))
}