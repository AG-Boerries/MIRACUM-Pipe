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
  # Calculate Percentage of Targeted Bases with more reads than cutoff (WES = 8, Panel = 20)
  if (protocol %in% c("somatic", "somaticGermline", "tumorOnly")) {
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

reads <- function(tfile, gfile = NULL){
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
  treads_tab <- readSamtoolsStats(tfile, "SN")$SN
  ntreads <- round(as.numeric(treads_tab$value[which(as.character(treads_tab$description) == "reads properly paired:")])/1000000)
  tin_size <- as.character(treads_tab$value[which(as.character(treads_tab$description) == "insert size average:")])
  tin_sd <- as.character(treads_tab$value[which(as.character(treads_tab$description) == "insert size standard deviation:")])

  ngreads <- NULL
  gin_size <- NULL
  gin_sd <- NULL
  
  if (!is.null(gfile)) {
    greads_tab <- readSamtoolsStats(gfile, "SN")$SN
    ngreads <- round(as.numeric(treads_tab$value[which(as.character(treads_tab$description) == "reads properly paired:")])/1000000)
    gin_size <- as.character(treads_tab$value[which(as.character(treads_tab$description) == "insert size average:")])
    gin_sd <- as.character(treads_tab$value[which(as.character(treads_tab$description) == "insert size standard deviation:")])
  }
  
  return(list(nRT = ntreads, tin = tin_size, tin_sd = tin_sd,
              nRG = ngreads, gin = gin_size, gin_sd = gin_sd))
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
  if (protocol %in% c("somatic", "somaticGermline", "tumorOnly")) {
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

quality_check <- function(path, nsamples, protocol){
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

#' Parse samtools stat output
#' Source: https://github.com/mcjmigdal/sumsamstats/blob/master/R/parse.R
#'
#' readSamtoolsStats parses output of samtools stat making it easy to
#' work with it in R.
#'
#' @param file the name of the file which the data are to be read from. File must be
#' of type \code{character}, providing path to input file. If it does not provide an
#' absolute path(s), the file name is relative to the current working directory, \code{getwd()}.
#'
#' @param section A character vector of patterns to extract from the report. Each section of samtools stats report is uniquley
#' labeled, for example Summary Numbers are labeled as SN. By default all sections are read.
#'
#' @return list of data frames holding data from different parts of \code{samtools stat} output.
#' Sections names are: SN (summary numbers), FFQ (first fragment qualities), LFQ (last fragment qualities),
#' GCF (GC Content of first fragments), GCL (GC content of last fragments), GCC (ACGT content per cycle),
#' IS (insert size), RL (read lengths), ID (indel distribution), IC (indels per cycle), COV (coverage distribution),
#' GCD (GC-depth)
#'
#' @examples
#' readSamtoolsStats(file = system.file('extdata', 'sample1', package = 'sumsamstats'))
#'
#' @export
readSamtoolsStats <- function(file, section = c("SN", "FFQ", "LFQ", "GCF", "GCL", "GCC", "IS", "RL", "ID", "IC", "COV", "GCD")) {
    if (!is.character(file)) {
        stop("File argument must be of class character!\n")
    }
    if (!file.exists(file)) {
        stop(paste0("File '", file, "' doesn't exists!\n"))
    }
    if (all(!section %in% c("SN", "FFQ", "LFQ", "GCF", "GCL", "GCC", "IS", "RL", "ID", "IC", "COV", "GCD"))) {
        stop(paste0("Parsing '", section, "' is not supported!"))
    }

    grepTable <- list(
      SN = list(
        columns = c(2, 3),
        col.names = c("description", "value")
        ),
      FFQ = list(
        columns = 2:43,
        col.names = c("cycle", paste0("Qual", 1:41))
        ),
      LFQ = list(
        columns = 2:43,
        col.names = c("cycle", paste0("Qual", 1:41))
        ),
      GCF = list(
        columns = c(2, 3),
        col.names = c("GC", "count")
        ),
      GCL = list(
        columns = c(2, 3),
        col.names = c("GC", "count")
        ),
      GCC = list(
        columns = 2:8,
        col.names = c("cycle", "A", "C", "G", "T", "N", "O")
        ),
      IS = list(
        columns = c(2, 3),
        col.names = c("insert_size", "pairs_total")
        ),
      RL = list(
        columns = c(2, 3),
        col.names = c("read_length", "count")
        ),
      ID = list(
        columns = c(2, 3, 4),
        col.names = c("length", "number_of_insertions", "number_of_deletions")
        ),
      IC = list(
        columns = c(2, 3, 4, 5, 6),
        col.names = c("cycle", "number_of_insertions_fwd", "number_of_insertions_rwd",
                      "number_of_deletions_fwd", "number_of_deletions_rwd")
        ),
      COV = list(
        columns = c(3, 4),
        col.names = c("coverage", "bases")
        ),
      GCD = list(
        columns = c(2, 3, 4, 5, 6, 7, 8),
        col.names = c("GC", "unique_sequence_percentiles", "10th", "25th", "50th", "75th", "90th")
        )
      )

    stats <- list()
    inputFile <- readLines(file)
    for (i in 1:length(section)) {
        stats[[section[i]]] <- .grepData(inputFile, section[i], columns = grepTable[[section[i]]]$columns, col.names = grepTable[[section[i]]]$col.names)
    }
    return(stats)
}

#' Helper function for parsing output of samtools stats
#' Source: https://github.com/mcjmigdal/sumsamstats/blob/master/R/parse.R
#'
#' As readSamtoolsStats parses output of samtools stat this function is used to grep relevant
#' parts of the output.
#'
#' @param input A character vector containing output of samtools stats, each line as one element.
#'
#' @param columns A integer vector containing numbers of columns to extract from specified section of samtools
#' stats report.
#'
#' @param section Pattern to extract from the report. Each section of samtools stats report is uniquley
#' labeled, for example Summary Numbers are labeled as SN.
#'
#' @param col.names A character vector containing names of extracted columns. Optional.
#'
#' @return data frame holding data from selected section of \code{samtools stat} output.
#'
#' @examples
#' file <- system.file('extdata', 'sample1', package = 'sumsamstats')
#' inputFile <- readLines(file)
#' .grepData(input=inputFile, section='SN', columns=c(2,3), col.names=c('description', 'value'))
#'
.grepData <- function(input, section, columns, col.names = "") {
  section <- paste("^", section, sep = "")
  handle <- strsplit(input[grep(pattern = section, input)], split = "\t")
  if (length(handle) == 0) {
    stop("Could not find pattern '", section, "'!")
  }
  handle <- data.frame(vapply(columns, function(i) {
    vapply(handle, `[`, i, FUN.VALUE = character(1))
  }, FUN.VALUE = character(length(handle))), stringsAsFactors = FALSE)
  if (length(col.names > 0)) {
    colnames(handle) <- col.names
  }
  return(handle)
}