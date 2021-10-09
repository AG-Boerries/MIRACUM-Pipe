#################
# Coverage Plot #
coverage_plot <- function(path, outfilePDF, protocol) {
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
  # Calculate Percentage of Targeted Bases with more reads than cutoff (WES = 8, Panel = 20)
  if (protocol != "panelTumor") {
    min_cov <- 8
    mit_cov <- 40
  } else {
    min_cov <- 20
    mit_cov <- 100
  }
  perc <- list()
  for (i in 1:length(files)) {
    perc_mc <- 1 - sum(cov[[i]][which(cov[[i]][, 2] < min_cov), 5])
    perc_min <- 1 - sum(cov[[i]][which(cov[[i]][, 2] < mit_cov), 5])
    perc[[i]] <- c(perc_mc, perc_min)
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
  return(list(cov = cov, perc = perc, labs = labs, files = files))
}

reads <- function(tfile, gfile) {
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
  treads_tab <- read.table(file = tfile, sep = "\t", skip = 7, nrows = 38, fill = TRUE)
  id <- which (as.character(treads_tab$V2) == "reads properly paired:")
  treads <- as.character(treads_tab$V3[id])
  treads <- as.numeric(treads)/1000000
  ntreads <- round(treads)
  tin_size <- as.character(treads_tab$V3[which(as.character(treads_tab$V2) == "insert size average:")])
  tin_sd <- as.character(treads_tab$V3[which(as.character(treads_tab$V2) == "insert size standard deviation:")])

  greads_tab <- read.table(file = gfile, sep = "\t", skip = 7, nrows = 38, fill = TRUE)
  id <- which (as.character(greads_tab$V2) == "reads properly paired:")
  greads <- as.character(greads_tab$V3[id])
  greads <- as.numeric(greads)/1000000
  ngreads <- round(greads)
  gin_size <- as.character(greads_tab$V3[which(as.character(greads_tab$V2) == "insert size average:")])
  gin_sd <- as.character(greads_tab$V3[which(as.character(greads_tab$V2) == "insert size standard deviation:")])
  
  
  return(list(nRT = ntreads, tin = tin_size, tin_sd = tin_sd,
              nRG = ngreads, gin = gin_size, gin_sd = gin_sd))
}

treads <- function(tfile) {
  treads_tab <- read.table(file = tfile, sep = "\t", skip = 7, nrows = 31, fill = TRUE)
  id <- which (as.character(treads_tab$V2) == "reads properly paired:")
  treads <- as.character(treads_tab$V3[id])
  treads <- as.numeric(treads)/1000000
  ntreads <- round(treads)
  tin_size <- as.character(treads_tab$V3[which(as.character(treads_tab$V2) == "insert size average:")])
  tin_sd <- as.character(treads_tab$V3[which(as.character(treads_tab$V2) == "insert size standard deviation:")])
  
  return(list(nRT = ntreads, tin = tin_size, tin_sd = tin_sd,
              nRG = NULL, gin = NULL, gin_sd = NULL))
}

coverage_exon <- function(path, protocol = protocol){
  #' Coverage Exons
  #'
  #' @description Coverage Exons
  #'
  #' @param path string. Path to data
  #' @param protocol string. Tumor_Normal oder Tumor_Only
  #'
  #' @return list of
  #' @return cov numerical. Mean coverage
  #' @return labs vector of strings. IDs for coverage files
  #' @return perc vector. Vector containing percentage of targeted exons with at least
  #' @return 8 or 40 reads for Tumor_Normal, 20 or 100 reads for Tumor_Only.
  #'
  #'

  # Get a list of the bedtools output files you'd like to read in
  print(files <- list.files(path = path, pattern = "exons.txt$"))
  print(labs <- gsub("_coverage.exons.txt", "", files, perl = TRUE))
  files <- paste(path, files, sep = "/")

  # Create lists to hold coverage and cumulative coverage for each alignment,
  # and read the data into these lists.
  cov <- list()
  cov_cumul <- list()
  for (i in 1:length(files)) {
    cov[[i]] <- read.table(files[i])
    cov_cumul[[i]] <- 1 - cumsum(cov[[i]][, 5])
  }
  # Calculate Percentage of Targeted Bases with more reads than cutoff (WES = 8, Panel = 20)
  if (protocol != "panelTumor") {
    min_cov <- 8
    mit_cov <- 40
  } else {
    min_cov <- 20
    mit_cov <- 100
  }
  perc <- list()
  for (i in 1:length(files)) {
    perc_mc <- 1 - sum(cov[[i]][which(cov[[i]][, 2] < min_cov), 5])
    perc_min <- 1 - sum(cov[[i]][which(cov[[i]][, 2] < mit_cov), 5])
    perc[[i]] <- c(perc_mc, perc_min)
  }
return(list(cov = cov, labs = labs, perc = perc))
}

quality_check <- function(path, nsamples, protocol) {
  gc_content <- rep(0, times = length(nsamples))
  mean_QC <- rep(0, times = length(nsamples))
  for (i in 1:length(nsamples)){
    if (protocol %in% c("somatic", "somaticGermline")){
      filename <- paste0(path, "/", nsamples[i], "_output.sort.filtered.rmdup.realigned.fixed.recal_fastqc/fastqc_data.txt")
    }
    if (protocol == "tumorOnly") {
      filename <- paste0(path, "/", nsamples[i], "_output.sort.rmdup.realigned.fixed.recal_fastqc/fastqc_data.txt")
    }
    if (protocol == "panelTumor") {
      filename <- paste0(path, "/", nsamples[i], "_output.sort.realigned.fixed.recal_fastqc/fastqc_data.txt")
    }
    sectionEndings <- which(str_detect(readLines(file(filename, "r", blocking = F)), ">>END_MODULE") == TRUE)
    fastq_data <- read.table(file = filename, skip = 2 , nrows = 7, sep = "\t")
    gc_content[i] <- as.character(fastq_data[which(fastq_data$V1 == "%GC"), 2])
    fastq_data <- read.table(file = filename, skip = 12 , nrows = sectionEndings[2]-14, sep = "\t")
    readlength <- convert_readlength(fastq_data$V1)
    sum_QC <- sum(readlength * fastq_data$V2)
    mean_QC[i] <- sum_QC/sum(readlength)
  }
  return(list(labs = nsamples, gc_content = gc_content, mean_QC = mean_QC)) 
}

convert_readlength <- function(x) {
  split <- str_split(x, "-", simplify = T)
  start <- as.numeric(split[,1])
  end <- as.numeric(split[,2])
  end[is.na(end)] = start[is.na(end)]
  return((end-start)+1)
}