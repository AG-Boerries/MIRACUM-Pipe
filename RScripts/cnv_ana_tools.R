####################################################
# Assigning Significance to Copy Number Variations #
# ControlFreec
assess_significance <- function(ratio_file, cnvs_file, outfile){
  #' Assess Significance
  #'
  #' @description Print Number of Somatic Mutations
  #'
  #' @param ratio_file string. Filename for ratio
  #' @param cnvs_file string. Filename for CNVs
  #' @param outfile string. Filename for output
  #'
  #' @details Adds p-value of Wilcoxon Rank Sum Test and p-value of Kolmogorov
  #' @details Smirnov to table of copy number variations and writes results in
  #' @details outfile.
  dataTable <-read.table(ratio_file, header = TRUE);
  ratio <-data.frame(dataTable)
  dataTable <- read.table(cnvs_file, header = FALSE)
  cnvs <- data.frame(dataTable)
  ratio$Ratio[which(ratio$Ratio == - 1)] = NA
  cnvs.bed = GRanges(cnvs[, 1], IRanges(cnvs[, 2], cnvs[, 3]))  
  ratio.bed = GRanges(ratio$Chromosome, IRanges(ratio$Start, ratio$Start),
                      score=ratio$Ratio)
  overlaps <- subsetByOverlaps(ratio.bed, cnvs.bed)
  normals <- setdiff(ratio.bed, cnvs.bed)
  normals <- subsetByOverlaps(ratio.bed, normals)
  numberOfCol=length(cnvs)
  for (i in c(1:length(cnvs[,1]))) {
    values <- score(subsetByOverlaps(ratio.bed, cnvs.bed[i]))
    
    W <- function(values, normals){
      resultw <- try(wilcox.test(values,score(normals)), silent = TRUE)
    if (class(resultw) == "try-error") {
      return(list("statistic" = NA,"parameter" = NA,"p.value" = NA,
                  "null.value" = NA, "alternative" = NA, "method" = NA,
                  "data.name" = NA))
    } else {
        resultw
    }}
    
    KS <- function(values,normals){
      resultks <- try(ks.test(values,score(normals)), silent = TRUE)
    if (class(resultks) == "try-error") {
      return(list("statistic" = NA, "p.value" = NA, "alternative" = NA,
                  "method" = NA, "data.name" = NA))
    } else {
        resultks
    }
    }
    cnvs[i,numberOfCol + 1] = W(values,normals)$p.value
    cnvs[i,numberOfCol + 2] = KS(values,normals)$p.value
  }
  if (numberOfCol == 5) {
    names(cnvs) = c("chr", "start", "end", "copy number", "status",
                    "WilcoxonRankSumTestPvalue", "KolmogorovSmirnovPvalue")  
  }
  if (numberOfCol == 7) {
    names(cnvs) = c("chr", "start", "end", "copy number", "status", "genotype",
                    "uncertainty", "WilcoxonRankSumTestPvalue",
                    "KolmogorovSmirnovPvalue")  
  }
  if (numberOfCol == 9) {
    names(cnvs) = c("chr", "start", "end", "copy number", "status", "genotype",
                    "uncertainty", "somatic/germline", "precentageOfGermline",
                    "WilcoxonRankSumTestPvalue", "KolmogorovSmirnovPvalue")  
  }
  write.table(cnvs, file = outfile, sep = "\t", quote = F, row.names = F)
}

make_cnv_graph <- function(ratio_file, ploidity = '2', outfile_plot,
                         outfile_ideogram){
  #' Make CNV Graph
  #'
  #' @description Plot CNVs on Genome
  #'
  #' @param ratio_file string. Filename of ratio file
  #' @param ploidity numerical. Ploidity (default: 2)
  #' @param outfile_plot string. Filename of output file
  #' @param outfile_ideogram string. Filename of ideogram 
  #'
  #' @details This function produces two plots that visualizes the CNVs on the
  #' @details genome.
  require(gtrellis)
  require(circlize)
  require(ComplexHeatmap)

  dataTable <-read.table(ratio_file, header=TRUE);

  ratio <-data.frame(dataTable)
  ploidy <- type.convert(ploidity)

  maxLevelToPlot <- 3
  for (i in c(1:length(ratio$Ratio))) {
    if (ratio$Ratio[i] > maxLevelToPlot) {
      ratio$Ratio[i] = maxLevelToPlot;
    }
  }

  ratio_new <- ratio
  ratio_new$Value <- ratio_new$Ratio*ploidy
  ratio_new$Chromosome <- paste('chr', ratio_new$Chromosome, sep = '')
  
  pdf(file = outfile_plot, width = 15, height = 15, pointsize = 20)
  gtrellis_layout(n_track = 3, nrow = 5, byrow = F, equal_width = FALSE,
                  track_axis = c(FALSE, TRUE, FALSE), 
                  track_height = unit.c(2 * grobHeight(textGrob("chr1")), 
                                        unit(1, "null"),
                                        grobHeight(textGrob("chr1"))),  
                  track_ylim = c(0, 1, 0, max(ratio_new$Value), 0, 1),
                  track_ylab = c("", "normalized copy number profile", ""))
  add_track(panel_fun = function(gr) {
    chr = get_cell_meta_data("name")  
    grid.rect(gp = gpar(fill = "#EEEEEE"))
    grid.text(chr)
  })
  add_track(ratio_new, track = 2, panel_fun = function(ratio_new) {
    x = ratio_new$Start
    y = ratio_new$Value
    grid.points(x,y,pch = ".", gp = gpar(col = colors()[88]))
  })
  tt <- which(ratio_new$CopyNumber>ploidy)
  tmp <- ratio_new[tt,]
  add_track(tmp, track = 2, panel_fun = function(tmp) {
    x = tmp$Start
    y = tmp$Value
    grid.points(x,y,pch = ".", gp = gpar(col = colors()[136]))
  })
  tt <- which(ratio_new$Ratio==maxLevelToPlot & ratio_new$CopyNumber>ploidy)
  tmp <- ratio_new[tt,]
  add_track(tmp, track = 2, panel_fun = function(tmp) {
    x = tmp$Start
    y = tmp$Value
    grid.points(x,y,pch = ".", gp = gpar(col = colors()[136], cex = 4))
  })  
  tt <- which(ratio_new$CopyNumber<ploidy & ratio_new$CopyNumber != -1)
  tmp <- ratio_new[tt,]
  add_track(tmp, track = 2, panel_fun = function(tmp) {
    x = tmp$Start
    y = tmp$Value
    grid.points(x,y,pch = ".", gp = gpar(col = colors()[461]))
  })
  cytoband_df = circlize::read.cytoband()$df
  add_track(cytoband_df, panel_fun = function(gr) {
    cytoband_chr = gr
    grid.rect(cytoband_chr[[2]], unit(0, "npc"),
              width = cytoband_chr[[3]] - cytoband_chr[[2]],
              height = unit(1, "npc"), default.units = "native", hjust = 0,
              vjust = 0,
              gp = gpar(fill = circlize::cytoband.col(cytoband_chr[[5]])))
    grid.rect(min(cytoband_chr[[2]]), unit(0, "npc"),
              width = max(cytoband_chr[[3]]) - min(cytoband_chr[[2]]),
              height = unit(1, "npc"),
              default.units = "native", hjust = 0, vjust = 0,
              gp = gpar(fill = "transparent"))
  })
  dev.off()
  
  col_fun = colorRamp2(c(1, 2, 3), c(colors()[461], colors()[88],
                                     colors()[136]))
  cm = ColorMapping(col_fun = col_fun)
  lgd = color_mapping_legend(cm, plot = FALSE, title = "Value")
  
  pdf(file = outfile_ideogram, width = 12, height = 12, pointsize = 20)
  gtrellis_layout(n_track = 1, ncol = 1, track_axis = FALSE,
                  xpadding = c(0.1, 0), gap = unit(4, "mm"), border = FALSE,
                  asist_ticks = FALSE, add_ideogram_track = TRUE, 
                  ideogram_track_height = unit(3, "mm"))
  
  tt <- which(ratio_new$CopyNumber>ploidy
              | (ratio_new$CopyNumber<ploidy & ratio_new$CopyNumber != -1))
  tmp <- ratio_new[tt,]
  
  add_track(track = 1, tmp, panel_fun = function(gr) {
    grid.rect(gr$Start, unit(0.2, "npc"), unit(1, "mm"), unit(0.8, "npc"),
              hjust = 0, vjust = 0, default.units = "native",
              gp = gpar(fill = col_fun(gr$CopyNumber), col = NA))
  })
  add_track(track = 2, clip = FALSE, panel_fun = function(gr) {
    chr = get_cell_meta_data("name")
    if(chr == "chrY") {
      grid.lines(get_cell_meta_data("xlim"), unit(c(0, 0), "npc"), 
                 default.units = "native")
    }
    grid.text(chr, x = 0, y = 0, just = c("left", "bottom"))
  })
  dev.off()
}

del_dup_query <- function(da.fr){
  require(RMySQL)
  require(doMC)
  chroms <- unlist(da.fr[[1]])
  ucscchroms <- unlist(lapply(chroms, function(i) paste0("chr",i)))
  start_list <- unlist(da.fr[[2]])
  stop_list <- unlist(da.fr[[3]])
  out <- mclapply(1:length(ucscchroms), function(i){
    con<- dbConnect(RMySQL::MySQL(),
                    host="genome-euro-mysql.soe.ucsc.edu",
                    user = "genome",
                    password = '',
                    port = 3306,
                    dbname = "hg19"
    )
    chrom <-ucscchroms[i]
    start <-start_list[i]
    stop <- stop_list[i]
    tmp <- paste0("Select name2 FROM refGene WHERE chrom='",chrom,"'AND ((cdsStart BETWEEN'",start,"'AND'",stop,"')OR (cdsEnd BETWEEN'",start,"'AND'",stop,"'))")
    res <- dbSendQuery(con, tmp)
    da.fr <- dbFetch(res, n =-1)
    da.fr <- unique(da.fr)
    dbClearResult(res)
    da.fr <- sort(unlist(as.list(da.fr)))
    names(da.fr) <- NULL
    dbDisconnect(con)
    return(da.fr)
  }, mc.cores = detectCores())

  name_vector <- unlist(lapply(1:length(ucscchroms), function(i) paste0(ucscchroms[i],":",start_list[i],"-",stop_list[i])))
  names(out) <- name_vector
  return(out)
}


cnv_annotation <- function(cnv_pvalue_txt, outfile, outfile_onco, outfile_tumorsuppressors, dbfile, path_data, path_script){
  #' CNV Annotation
  #'
  #' @description Annotate the CNV Regions
  #'
  #' @param cnv_pvalue_txt string. Filename of CNV regions file
  #' @param outfile string. Filename of CNV list
  #' @param outfile_onco string. Filename of Oncogene list
  #' @param outfile_tumorsuppressors string. Filename of Tumorsuppressor list
  #' @param path_data string. Directory of the databases
  #' @param path_script string. Directory of required scripts
  #'
  #' @return returns list of
  #' @return CNVsAnnotated dataframe. Table of genes with CNVs
  #' @return CNVOncogenes dataframe. Table of oncogenes with CNVs
  #' @return CNVTumorSuppressors dataframe. Table of tumorsuppressor genes with
  #' @return CNVs
  #'
  #' @details The CNV regions are correlated to the affected genes. A table
  #' @details of these genes is generated. Each gene has also a column with the
  #' @details total copy number. The table is divided in two subtables
  #' @details containing oncogenes and tumorsuppressor genes. 
  x <- read.table(cnv_pvalue_txt, header = T, sep = "\t")
  x <- data.frame(x)

  db <- read.delim(dbfile, header = T, sep = "\t", colClasses = "character")
  db <- data.frame(db)
  ts <- which(db$Is.Tumor.Suppressor.Gene == "Yes")
  ts <- db[ts, ]
  og <- which(db$Is.Oncogene == "Yes")
  onc <- db[og, ]
  
  x$genes <- "."
  x$TumorSuppressor <- "."
  x$Oncogene <- "."
  x$Length <- "."
  not.significant <- c()
  
if(x$WilcoxonRankSumTestPvalue[i] < 0.05
       & x$KolmogorovSmirnovPvalue[i] < 0.05
       & !is.na(x$WilcoxonRankSumTestPvalue[i])
       & !is.na(x$KolmogorovSmirnovPvalue[i])) {
      location <- list(x$chr[i], x$start[i], x$end[i])
      
      query <- del_dup_query(location)
      query2 <- unlist(query)
      x$genes[i] <- paste(query2, collapse = ", ")
      
      # Test for Tumor Suppressors
      id.ts <- which(ts$Hugo.Symbol %in% query2)
      if(sum(id.ts) > 0) {
        x$TumorSuppressor[i] <- paste(ts$Hugo.Symbol[id.ts], collapse = ",")
      }
      id.onc <- which(onc$Hugo.Symbol %in% query2)
      if(sum(id.onc) > 0) {
        x$Oncogene[i] <- paste(onc$Hugo.Symbol[id.onc], collapse = ",")
      }
      x$Length[i] <- x$end[i] - x$start[i]
    } else {
      not.significant <- c(not.significant, i)
    }
  }
  
  if (length(not.significant) > 0){
    x1 <- x[-not.significant, ]
  } else {
    x1 <- x
  }
  
  write.xlsx(x1, outfile, row.names = F)
  onco <- x1[x1$Oncogene != ".", c("chr", "copy.number", "status", "Oncogene")]
  write.xlsx(onco, outfile_onco, row.names = F)
  tumorSupp <- x1[x1$TumorSuppressor != ".", c("chr", "copy.number", "status", "TumorSuppressor")]
  write.xlsx(tumorSupp, outfile_tumorsuppressors, row.names = F)
  
  return(list(CNVsAnnotated = x1, CNVOncogenes = onco, CNVTumorSuppressors = tumorSupp))
}



cnv_processing <- function(cnv_file, targets,
                           outfile_loss, outfile_gain, go.bp,
                           path_data, path_script){
  #' CNV Processing
  #'
  #' @description Copy Number Analysis
  #'
  #' @param cnv_file string. Filename of genes with a CNV
  #' @param targets string. Filename of targeted genes
  #' @param outfile_loss string. Filename of genes with loss
  #' @param outfile_gain string. Filename of genes with gain
  #' @param go.bp dataframe. GO Database
  #' @param path_data string. Directory of the databases
  #' @param path_script string. Directory of required scripts
  #'
  #' @return returns list of
  #' @return loss_go dataframe. Table of pathways affected by genes with losses
  #' @return gain_go dataframe. Table of pathways affected by genes with gains
  #'
  #' @details The genes with CNVs are playing a role in some pathways. The
  #' @details dataset of GO is used to find these pathways. They are stored
  #' @details in the output files.
  cnv <- read.xlsx(cnv_file)
  
  cnv.genes <- c()
  cnv$genes <- as.character(cnv$genes)
  cnv$status <- as.character(cnv$status)
  cnv$chr <- as.character(cnv$chr)
# LOSS 
  for (i in 1:length(cnv$status)) {
    if(cnv$status[i] == "loss") {
      s <- strsplit(cnv$genes[i], split = ",", fixed = TRUE)
      cnv.genes <- c(cnv.genes, s[[1]])
    }
  }
  cnv.genes <- unique(cnv.genes)
  cnv.genes <- cnv.genes[!is.na(cnv.genes)]
  prep <- prep_pwa(targets, cnv.genes)
  loss_go <- get_terms(go.bp, outfile_loss, prep$de_genes, prep$universe)
  ########
  # GAIN #
  cnv.genes <- c()
  for(i in 1:length(cnv$status)) {
    if(cnv$copy.number[i] > 3) {
      s <- strsplit(cnv$genes[i], split = ",", fixed = TRUE)
      cnv.genes <- c(cnv.genes, s[[1]])
    }
  }
  cnv.genes <- unique(cnv.genes)
  cnv.genes <- cnv.genes[!is.na(cnv.genes)]
  prep <- prep_pwa(targets, cnv.genes)
  gain_go <- get_terms(go.bp, outfile_gain, prep$de_genes, prep$universe)
  return(list(loss_go = loss_go, gain_go = gain_go))
}

cnv_table <- function(cnvs){
  #' CNV Table
  #'
  #' @description Write CNV Table
  #'
  #' @param cnvs dataframe. Table of regions with CNVs
  #'
  #' @return output dataframe. Table of genes with CNVs
  #'
  #' @details Regions with CNVs are correlated to the corresponding genes.
  output <- data.frame()
  for(i in 1:dim(cnvs)[1]){
    genesSplitList <- strsplit(cnvs$genes, split = ',')
    entryLength <- lapply(genesSplitList, length)
    chr <- rep(cnvs$chr[i],entryLength[[i]])
    start <- rep(cnvs$start[i],entryLength[[i]])
    end <- rep(cnvs$end[i],entryLength[[i]])
    genes <- genesSplitList[[i]]
    copyNumber <- rep(cnvs$copy.number[i],entryLength[[i]])
    status <- rep(cnvs$status[i],entryLength[[i]])
    tmp <- data.frame(Gene = genes, Chr = chr, Start = start, End = end, CopyNumber = copyNumber, Status = status)
    output <- rbind(output, tmp)
  }
  return(output)
}

# check cnvs for dna damage response
cnv_dna_damage<- function(input, outfile_dna_damage, db){
  #' CNV DNA Damage Response
  #'
  #' @description Get DNA Damage Response genes with CNVs
  #'
  #' @param input dataframe. Table of genes with CNVs
  #' @param outfile_dna_damage string. Filename for output file
  #' @param db dataframe. Table of DNA Damage Response genes
  #'
  #' @return l dataframe. Table of DNA Damage Reponse genes affected by CNVs
  #'
  #' @details Genes with CNVs could play a role in the DNA Damage Response
  #' @details pathway. This function extracts these genes and stores them
  #' @details in the output file.
  dnad <- read.delim(db, header = FALSE, sep = "\t", colClasses = "character")
  idx <- which (input$Gene %in% dnad$V1)
  l <-  list()
  gene <- input$Gene[idx]
  copynumber <- input$CopyNumber[idx]
  status <- input$Status[idx]
  l <- cbind(gene, copynumber, status)
  write.xlsx(l, file = outfile_dna_damage, row.names = FALSE, col.names = TRUE)
  return(l)
}
