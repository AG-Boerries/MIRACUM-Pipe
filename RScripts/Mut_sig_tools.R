#' Wrapper to compute confidence intervals for a cohort
#' 
#' Wrapper function around \link{confIntExp}, which is applied to every
#' signature/sample pair in a cohort. The extracted upper and lower bounds of
#' the confidence intervals are added to the input data which is reordered and
#' melted in order to prepare for visualization with ggplot2.
#'
#' @param in_catalogue_df
#'  Input numerical data frame of the mutational catalog of the cohort to be
#'  analyzed.
#' @param in_sig_df
#'  Numerical data frame of the signatures used for analysis.
#' @param in_exposures_df
#'  Input numerical data frame of the exposures computed for the cohort to be
#'  analyzed.
#' @param in_sigLevel
#'  Significance level, parameter passed to \link{confIntExp}.
#' @param in_delta
#'  Inflation parameter for the alternative model, parameter passed on to
#'  \link{confIntExp}
#' @param in_pdf
#'  Probability distribution function, parameter passed on to \link{confIntExp},
#'  if NULL assumed to be normal distribution.
#'
#' @return A melted data frame.
#' @import magrittr
#' @import dplyr
#' @export
#'
#' @examples
#'  library(BSgenome.Hsapiens.UCSC.hg19)
#'  data(lymphoma_test)
#'  data(lymphoma_cohort_LCD_results)
#'  data(sigs)
#'  word_length <- 3
#'  temp_list <- create_mutation_catalogue_from_df(
#'    lymphoma_test_df,this_seqnames.field = "CHROM",
#'    this_start.field = "POS",this_end.field = "POS",
#'    this_PID.field = "PID",this_subgroup.field = "SUBGROUP",
#'    this_refGenome = BSgenome.Hsapiens.UCSC.hg19,
#'    this_wordLength = word_length)
#'  lymphoma_catalogue_df <- temp_list$matrix
#'  lymphoma_PIDs <- colnames(lymphoma_catalogue_df)
#'  data("lymphoma_cohort_LCD_results")
#'  lymphoma_exposures_df <-
#'    lymphoma_Nature2013_COSMIC_cutoff_exposures_df[, lymphoma_PIDs]
#'  lymphoma_sigs <- rownames(lymphoma_exposures_df)
#'  lymphoma_sig_df <- AlexCosmicValid_sig_df[, lymphoma_sigs]
#'  lymphoma_complete_df <- variateExp(in_catalogue_df = lymphoma_catalogue_df,
#'                                     in_sig_df = lymphoma_sig_df,
#'                                     in_exposures_df = lymphoma_exposures_df,
#'                                     in_sigLevel = 0.025, in_delta = 0.4)
#'  head(lymphoma_complete_df)
#'  lymphoma_complete_df$sample <- 
#'    factor(lymphoma_complete_df$sample, 
#'           levels = colnames(lymphoma_exposures_df)[
#'             order(colSums(lymphoma_exposures_df), decreasing = TRUE)])
#'  sig_colour_vector <- c("black", AlexCosmicValid_sigInd_df$colour)
#'  names(sig_colour_vector) <- 
#'    c("total", as.character(AlexCosmicValid_sigInd_df$sig))
#'  ggplot(data = lymphoma_complete_df,
#'         aes(x = sample, y = exposure, fill = sig)) +
#'    geom_bar(stat = "identity") +
#'    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
#'    facet_wrap(~sig, nrow = nrow(lymphoma_exposures_df) + 1) +
#'    theme_grey() +
#'    theme(panel.border = element_rect(fill = NA, colour = "black"),
#'          strip.background = element_rect(colour = "black"),
#'          legend.position = "none") +
#'    scale_fill_manual(values = sig_colour_vector)
#'  
variateExp <- function(in_catalogue_df, in_sig_df, in_exposures_df,
                       in_sigLevel = 0.05, in_delta = 0.4, in_pdf = NULL){
  require(magrittr)
  require(dplyr)
  temp_list <- 
    lapply(seq_len(ncol(in_exposures_df)), function(z){
      temp_list <-
        lapply(seq_len(ncol(in_sig_df)), 
               function(y){
                 confIntExp(in_ind = y, in_sigLevel = in_sigLevel,
                            in_delta = in_delta,
                            in_catalogue_vector = in_catalogue_df[, z],
                            in_sig_df = in_sig_df,
                            in_exposure_vector = in_exposures_df[, z], 
                            in_pdf = in_pdf)
               })
      names(temp_list) <- names(in_sig_df)
      return(list(upper = unlist(lapply(temp_list, function(x) x$upper)),
                  lower = unlist(lapply(temp_list, function(x) x$lower))))
    })
  names(temp_list) <- names(in_exposures_df)
  upper_df <- do.call(cbind, lapply(temp_list, function(x) x$upper))
  lower_df <- do.call(cbind, lapply(temp_list, function(x) x$lower))
  melt_df <- Reduce(function(x, y) left_join(x, y, by = c("Var1", "Var2")), 
                    lapply(list(exposures = in_exposures_df,
                                lower = lower_df, 
                                upper = upper_df),
                           function(x) melt(as.matrix(x))))
  names(melt_df) <- c("sig", "sample", "exposure", "relLower", "relUpper")
  total_df <- data.frame(sig = "total", 
                         sample = names(in_exposures_df),
                         exposure = colSums(in_exposures_df),
                         relLower = 1, relUpper = 1)
  melt_df <- rbind(melt_df, total_df)
  melt_df %<>% mutate(lower = exposure * relLower,
                      upper = exposure * relUpper)
  melt_df$sig <- factor(melt_df$sig, levels = c("total",
                                                names(in_sig_df)))
  return(melt_df)
}

#' Compute confidence intervals
#' 
#' Compute confidence intervals using the (log-)likelihood ratio test,
#' primarily for one input sample.
#'
#' @param in_ind
#'  Index of the input signature to be variated.
#' @param in_sigLevel
#'  Significance leve (one-sided)
#' @param in_delta
#'  Inflation parameter for the alternative model.
#' @param in_exposure_vector
#'  Exposure vector computed for the input sample.
#' @param ... Input parameters passed on to \link{variateExpSingle}.
#'
#' @return
#'  A list with entries
#'  \itemize{
#'    \item \code{upper}: Upper bound of the confidence interval
#'    \item \code{lower}: Lower bound of the confidence interval
#'  }
#' @importFrom pracma newtonsys
#' 
#' @export
#'
#' @examples
#'  library(BSgenome.Hsapiens.UCSC.hg19)
#'  data(lymphoma_test)
#'  data(lymphoma_cohort_LCD_results)
#'  data(sigs)
#'  word_length <- 3
#'  temp_list <- create_mutation_catalogue_from_df(
#'    lymphoma_test_df,this_seqnames.field = "CHROM",
#'    this_start.field = "POS",this_end.field = "POS",
#'    this_PID.field = "PID",this_subgroup.field = "SUBGROUP",
#'    this_refGenome = BSgenome.Hsapiens.UCSC.hg19,
#'    this_wordLength = word_length)
#'  lymphoma_catalogue_df <- temp_list$matrix
#'  lymphoma_PIDs <- colnames(lymphoma_catalogue_df)
#'  data("lymphoma_cohort_LCD_results")
#'  lymphoma_exposures_df <-
#'    lymphoma_Nature2013_COSMIC_cutoff_exposures_df[, lymphoma_PIDs]
#'  lymphoma_sigs <- rownames(lymphoma_exposures_df)
#'  lymphoma_sig_df <- AlexCosmicValid_sig_df[, lymphoma_sigs]
#'  confIntExp(in_ind = 1, in_sigLevel = 0.05, in_delta = 0.4,
#'             in_exposure_vector = lymphoma_exposures_df[, 1],
#'             in_catalogue_vector = lymphoma_catalogue_df[, 1],
#'             in_sig_df = lymphoma_sig_df)
#'             
confIntExp <- function(in_ind = 1, in_sigLevel = 0.05, in_delta = 1,
                       in_exposure_vector = NULL, in_verbose = FALSE, ...){
  require(pracma)
  stopifnot(!is.null(in_exposure_vector))
  if(in_exposure_vector[in_ind] == 0){
    upper = 0
    lower = 0
  } else {
    current_fun <- function(x){
      variateExpSingle(in_ind = in_ind, in_factor = x, 
                       in_exposure_vector = in_exposure_vector,
                       verbose = in_verbose, ...) - in_sigLevel
    }
    start_delta <- in_delta
    upper <- 1
    while (upper <= 1 & start_delta < 10) {
      delta <- start_delta
      #cat("upper: start_delta: ", start_delta, "\n")
      while (upper <= 1 & delta > 1e-4) {
        #cat("\tupper: delta: ", delta, "\n")
        tryCatch({
          upper <- newtonsys(current_fun, 1 + delta)$zero
          if(upper <= 1) delta <- delta / 2
        },
        error = function(e){
          delta <<- delta / 2
          warning("\n\tupper; delta: ", delta)
        })
      }
      start_delta <- start_delta * 2
    }
    if(upper < 1){
      upper <- 1
      if(in_verbose) cat("No upper bound for CI could be found.\n")
    }
    start_delta <- in_delta
    lower <- 1
    while (lower >= 1 & start_delta < 10) {
      delta <- start_delta
      #cat("lower: start_delta: ", start_delta, "\n")
      while (lower >=1 & delta > 1e-4) {
        #cat("\tlower: delta: ", delta, "\n")
        tryCatch({
          lower <- newtonsys(current_fun, 1 - delta)$zero
          if(lower >= 1) delta <- delta / 2
        },
        error = function(e){
          delta <<- delta / 2
          warning("\n\tlower; delta: ", delta)
        })
      }
      start_delta <- start_delta * 2
    }
    if(lower > 1){
      lower <- 1
      if(in_verbose) cat("No lower bound for CI could be found.\n")
    }
  }
  return(list(upper = upper,
              lower = lower))
}

#' Wrapper for the likelihood ratio test
#' 
#' Application of the likelihood ratio test to mutational signatures, primarily
#' for one single sample.
#'
#' @param in_catalogue_vector
#'  Mutational catalog of the input sample.
#' @param in_sig_df
#'  Data frame encoding the signatures used for the analysis.
#' @param in_exposure_vector
#'  Exposure vector computed for the input sample.
#' @param in_ind
#'  Index specifying which signature among \code{in_sig_df} is to be tested.
#' @param in_factor
#'  Deviation factor of the altered alternative model.
#' @param in_pdf
#'  Probability distibution function, parameter passed on to
#'  \link{logLikelihood} and later to \link{computeLogLik}
#' @param verbose 
#'
#' @return
#' @export
#'
#' @examples
#'  library(BSgenome.Hsapiens.UCSC.hg19)
#'  data(lymphoma_test)
#'  data(lymphoma_cohort_LCD_results)
#'  data(sigs)
#'  word_length <- 3
#'  temp_list <- create_mutation_catalogue_from_df(
#'    lymphoma_test_df,this_seqnames.field = "CHROM",
#'    this_start.field = "POS",this_end.field = "POS",
#'    this_PID.field = "PID",this_subgroup.field = "SUBGROUP",
#'    this_refGenome = BSgenome.Hsapiens.UCSC.hg19,
#'    this_wordLength = word_length)
#'  lymphoma_catalogue_df <- temp_list$matrix
#'  lymphoma_PIDs <- colnames(lymphoma_catalogue_df)
#'  data("lymphoma_cohort_LCD_results")
#'  lymphoma_exposures_df <- 
#'    lymphoma_Nature2013_COSMIC_cutoff_exposures_df[, lymphoma_PIDs]
#'  lymphoma_sigs <- rownames(lymphoma_exposures_df)
#'  lymphoma_sig_df <- AlexCosmicValid_sig_df[, lymphoma_sigs]
#'  variateExpSingle(
#'    in_ind = 1,
#'    in_factor = 1.5, 
#'    in_catalogue_vector = lymphoma_catalogue_df[, 1],
#'    in_sig_df = lymphoma_sig_df,
#'    in_exposure_vector = lymphoma_exposures_df[, 1])
#'    
variateExpSingle <- function(in_catalogue_vector, in_sig_df,
                             in_exposure_vector, in_ind,
                             in_factor = 1, in_pdf = NULL, verbose = FALSE){
  ini_lsfit <- lsfit(as.matrix(in_sig_df), in_catalogue_vector, 
                     intercept = FALSE)
  delta_catalogue_vector <- 
    in_sig_df[, in_ind] * in_exposure_vector[in_ind]
  my_catalogue_vector <- in_catalogue_vector -
    in_factor * delta_catalogue_vector
  my_lsfit <- lsfit(as.matrix(in_sig_df[, -in_ind]), my_catalogue_vector, 
                    intercept = FALSE)
  return(logLikelihood(ini_lsfit$residuals, my_lsfit$residuals,
                       in_pdf = in_pdf, verbose = verbose)$p.value)
}

#' Compute a loglikelihood ratio test
#'
#' Compute a likelihood ratio test based on the loglikelihoods of the residuals
#' of two different models of the same data.
#'
#' @param in_1
#'  Residuals of model 1 of the input data.
#' @param in_2 
#'  Residuals of model 2 of the input data.
#' @param df_1
#'  Degrees of freedom of the input model 1. If either \code{df_1} or
#'  \code{df_2} is NULL, the difference between the degrees of freedom of the
#'  two models is assumed to be 1.
#' @param df_2 
#'  Degrees of freedom of the input model 2. If either \code{df_1} or
#'  \code{df_2} is NULL, the difference between the degrees of freedom of the
#'  two models is assumed to be 1.
#' @param in_pdf 
#'  Probability distribution function, passed on to \link{computeLogLik}, if
#'  NULL a normal distribution is used.
#' @param verbose 
#'
#' @return
#'  A list with entries
#'  \itemize{
#'    \item \code{statistic}:
#'     The test statistic
#'    \item \code{delta_df}:
#'     The difference in degrees of freedom between input model 1 and 2
#'    \item \code{p.value}:
#'     p value of the statistical test.
#'  }
#' @export
#'
#' @examples
#'  library(BSgenome.Hsapiens.UCSC.hg19)
#'  data(lymphoma_test)
#'  data(sigs)
#'  data(cutoffs)
#'  word_length <- 3
#'  temp_list <- create_mutation_catalogue_from_df(
#'    lymphoma_test_df,this_seqnames.field = "CHROM",
#'    this_start.field = "POS",this_end.field = "POS",
#'    this_PID.field = "PID",this_subgroup.field = "SUBGROUP",
#'    this_refGenome = BSgenome.Hsapiens.UCSC.hg19,
#'    this_wordLength = word_length)
#'  lymphoma_catalogue_df <- temp_list$matrix
#'  lymphoma_PIDs <- colnames(lymphoma_catalogue_df)
#'  current_sig_df <- AlexCosmicValid_sig_df
#'  current_sigInd_df <- AlexCosmicValid_sigInd_df
#'  current_cutoff_vector <- cutoffCosmicValid_rel_df[6, ]
#'  iniLCDList <- LCD_complex_cutoff(
#'    in_mutation_catalogue_df = lymphoma_catalogue_df[, 1, drop = F],
#'    in_signatures_df = current_sig_df,
#'    in_cutoff_vector = current_cutoff_vector,
#'    in_method = "relative", in_rescale = TRUE,
#'    in_sig_ind_df = current_sigInd_df)
#'  current_sig_df <- AlexCosmicValid_sig_df[, -9]
#'  current_sigInd_df <- AlexCosmicValid_sigInd_df[-9,]
#'  current_cutoff_vector <- cutoffCosmicValid_rel_df[6, -9]
#'  redLCDList <- LCD_complex_cutoff(
#'    in_mutation_catalogue_df = lymphoma_catalogue_df[, 1, drop = F],
#'    in_signatures_df = current_sig_df,
#'    in_cutoff_vector = current_cutoff_vector,
#'    in_method = "relative", in_rescale = TRUE,
#'    in_sig_ind_df = current_sigInd_df)
#'  logLikelihood(iniLCDList, redLCDList)
#'  
logLikelihood <- function(in_1, in_2,
                          df_1 = NULL, df_2 = NULL,
                          in_pdf = NULL, verbose = FALSE){
  if(is.null(df_1) & is.null(df_2)){
    delta_df <- 1
  } else {
    delta_df <- abs(df_1 - df_2)
  }
  if(inherits(in_1, "data.frame") & inherits(in_2, "data.frame")){
    my_res_1 <- as.numeric(as.matrix(in_1))
    my_res_2 <- as.numeric(as.matrix(in_2))
  } else if(inherits(in_1, "numeric") & inherits(in_2, "numeric")){
    my_res_1 <- in_1
    my_res_2 <- in_2
  } else if(inherits(in_1, "list") & inherits(in_2, "list")){
    my_res_1 <- as.numeric(as.matrix(in_1$residual_catalogue))
    my_res_2 <- as.numeric(as.matrix(in_2$residual_catalogue))
  } else {
    stop("Unvalid input data type.")
  }
  STATISTIC <-
    2 * (computeLogLik(my_res_1, in_pdf = in_pdf, verbose = verbose) - 
           computeLogLik(my_res_2, in_pdf = in_pdf, verbose = verbose))
  PVAL <- pchisq(STATISTIC, delta_df, lower.tail = FALSE)
  return(list(statistic = STATISTIC, delta_df = delta_df, p.value = PVAL))
}

#' Compute the loglikelihood
#'
#' @param in_vector
#'  Numeric vector of input values of which the loglikelihood is computed.
#' @param in_pdf
#'  Probability distribution function, if NULL a normal distribution is used.
#' @param verbose 
#'
#' @return A numeric value (sum of the logarithms of the likelihoods of the
#'  input vector)
#' @export
#'
#' @examples
#' 
computeLogLik <- function(in_vector, in_pdf = NULL, verbose = FALSE){
  if(verbose) cat("class(in_pdf): ", class(in_pdf), "\n")
  if(inherits(in_pdf, "function")){
    if(verbose) cat("Selection 1.\n")
    my_pdf <- function(x) in_pdf(x)
  } else if(inherits(in_pdf, "character")){
    if(in_pdf == "distrib"){
      if(verbose) cat("Selection 2.\n")
      fn <- approxfun(density(in_vector))
      my_pdf <- function(x) fn(x)
    } else if(in_pdf %in% c("zero", "0", "zeroNorm")){
      if(verbose) cat("Selection 3.\n")
      my_sd <- sd(in_vector)
      my_pdf <- function(x) dnorm(x, mean = 0, sd = my_sd)
    }
  } else {
    if(verbose) cat("Selection 4.\n")
    my_mean <- mean(in_vector)
    my_sd <- sd(in_vector)
    my_pdf <- function(x) dnorm(x, mean = my_mean, sd = my_sd)
  }
  log_likelihood <- log2(my_pdf(in_vector))
  return(sum(log_likelihood))
}

#### Plot?!

#' Plot exposures including confidence intervals
#'
#' Plot the exposures to extracted signatures including confidence intervals
#' computed e.g. by \link{variateExp}.
#'
#' @param in_complete_df
#'  Melted numeric input data frame e.g. as computed by \link{variateExp}
#' @param in_subgroups_df
#'  Data frame containing meta information on subgroup attribution of the
#'  samples in the cohort of interest.
#' @param in_sigInd_df
#'  Data frame with meta information on the signatures used in the analysis.
#'
#' @return The function doesn't return any value but plots instead.
#' @export
#'
#' @examples
plotExposuresConfidence <- function(in_complete_df,
                                    in_subgroups_df,
                                    in_sigInd_df){
  sig_colour_vector <- c("black", in_sigInd_df$colour)
  names(sig_colour_vector) <- c("total", as.character(in_sigInd_df$sig))
  in_complete_df$sample <- 
    factor(in_complete_df$sample, 
           levels = in_subgroups_df$PID[order(in_subgroups_df$index)])
  in_complete_df$subgroup <- 
    in_subgroups_df$subgroup[match(in_complete_df$sample, in_subgroups_df$PID)]
  exposure_plot <- in_complete_df[which(in_complete_df$sig != "total"),] %>%
    # exposure_plot <- in_complete_df %>%
    ggplot(aes(x = sample, y = exposure, fill = sig)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
    facet_grid(sig ~ subgroup, scales = "free_x", space = "free_x",
               switch = "x") +
    theme_grey() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          panel.border = element_rect(fill = NA, colour = "black"),
          strip.background = element_rect(colour = "black"),
          legend.position = "none") +
    scale_fill_manual(values = sig_colour_vector)
  subgroup_aggregate_df <- aggregate(col ~ subgroup, data = in_subgroups_df, 
                                     FUN = head, 1)
  index_aggregate_df <- aggregate(index ~ subgroup, data = in_subgroups_df, 
                                  FUN = mean)
  subgroup_colour_vector <-
    subgroup_aggregate_df$col[order(index_aggregate_df$index)]
  exposure_g <- ggplot_gtable(ggplot_build(exposure_plot))
  stripr <- which(grepl('strip-b', exposure_g$layout$name))
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', exposure_g$grobs[[i]]$grobs[[1]]$childrenOrder))
    exposure_g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- 
      subgroup_colour_vector[k]
    k <- k+1
  }
  total_p <- in_complete_df[which(in_complete_df$sig == "total"),] %>%
    ggplot(aes(x = sample, y = exposure, fill = sig)) +
    geom_bar(stat = "identity", fill = "black") +
    facet_grid(sig ~ subgroup, scales = "free_x", space = "free_x",
               switch = "x") +
    theme_grey() +
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black"),
          strip.background = element_rect(colour = "black"),
          strip.text.x = element_blank(),
          legend.position = "none")
  total_g <- ggplotGrob(total_p)
  all_g <- rbind(total_g, exposure_g)
  # grid.draw(exposure_g)
  grid.draw(all_g)
}

