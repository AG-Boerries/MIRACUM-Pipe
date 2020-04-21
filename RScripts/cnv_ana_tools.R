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

make_cnv_ideo_sig <- function(ratio_file, outfile_ideogram, protocol){
  #' Make CNV Graph for significant ones
  #'
  #' @description Plot significant CNVs on Genome
  #'
  #' @param ratio_file string. Filename of ratio file with p-values for CNVs
  #' @param outfile_ideogram string. Filename of ideogram
  #'
  #' @details This function produces a plot containing an ideogram that
  #' @details visualizes the significant CNVs on the human genome.
  require(gtrellis)
  require(circlize)
  require(ComplexHeatmap)
  if (protocol != "panelTumor"){
    dataTable <- read.table(ratio_file, header = FALSE, skip = 1)
    ratio <- data.frame(dataTable)
  } else {
    ratio <- as.data.frame(ratio_file)
  }
  colnames(ratio) <- c("chr", "start", "end", "CopyNumber", "status", "WRST", "KS")
  ratio_new <- ratio
  ratio_new$chr <- paste('chr', ratio_new$chr, sep = '')

  col_fun = colorRamp2(c(0, 1, 3, 4, 5, 6), c("darkblue", "mediumblue",
                                              "red", "firebrick3", "red4",
                                              "darkred"))
  cm = ColorMapping(col_fun = col_fun)
  lgd = color_mapping_legend(cm, plot = FALSE, title = "total Copy Number")

  ratio_new$CopyNumber[which(as.numeric(ratio_new$CopyNumber) > 6)] <- 6

  pdf(file = outfile_ideogram, width = 12, height = 12, pointsize = 20)
  gtrellis_layout(n_track = 1, ncol = 1, track_axis = FALSE,
                  xpadding = c(0.1, 0), gap = unit(4, "mm"), border = FALSE,
                  asist_ticks = FALSE, add_ideogram_track = TRUE,
                  ideogram_track_height = unit(3, "mm"),legend = lgd)

  tt <- which(as.numeric(ratio_new$CopyNumber) > 2
              | (as.numeric(ratio_new$CopyNumber) < 2 & as.numeric(ratio_new$CopyNumber) != -1))
  tmp <- ratio_new[tt,]

  add_track(track = 1, tmp, panel_fun = function(gr) {
    grid.rect(x = gr$start, y = unit(0.2, "npc"), width = (gr$end - gr$start),
              height = unit(0.8, "npc"), hjust = 0, vjust = 0,
              default.units = "native", gp = gpar(fill = col_fun(gr$CopyNumber),
                                                  col = NA))
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
#  te <- try(read.delim(file = sef_snp, sep = "\t", header = FALSE, comment.char = "#"), silent = TRUE)
#    if (inherits(te, 'try-error')) {
#  # Annotaion via MySQL database
#  for(i in 1:nrow(x)) {
#    cat("Processing CNV#", i, "\n")
#
#    if(x$WilcoxonRankSumTestPvalue[i] < 0.05 & x$KolmogorovSmirnovPvalue[i] < 0.05 & !is.na(x$WilcoxonRankSumTestPvalue[i]) & !is.na(x$KolmogorovSmirnovPvalue[i])) {
#      location <- list(x$chr[i], x$start[i], x$end[i])
#      query <- del_dup_query(location)
#      query2 <- unlist(query)
#      x$genes[i] <- paste(query2, collapse = ", ")
#      
#      # Test for Tumor Suppressors
#      id.ts <- which(ts$Hugo.Symbol %in% query2)
#      if(sum(id.ts) > 0) {
#        x$TumorSuppressor[i] <- paste(ts$Hugo.Symbol[id.ts], collapse = ",")
#      }
#      id.onc <- which(onc$Hugo.Symbol %in% query2)
#      if(sum(id.onc) > 0) {
#        x$Oncogene[i] <- paste(onc$Hugo.Symbol[id.onc], collapse = ",")
#      }
#      x$Length[i] <- x$end[i] - x$start[i]
#    } else {
#      not.significant <- c(not.significant, i)
#    }
#  }

# Try annotation first with SQL server if this fails try with biomaRt

  ensembl=useMart("ensembl", dataset="hsapiens_gene_ensembl")
  
  for(i in 1:nrow(x)) {
    cat("Processing CNV#", i, "\n")
    
    if(x$WilcoxonRankSumTestPvalue[i] < 0.05 & x$KolmogorovSmirnovPvalue[i] < 0.05 & !is.na(x$WilcoxonRankSumTestPvalue[i]) & !is.na(x$KolmogorovSmirnovPvalue[i])) {
      location <- list(x$chr[i], x$start[i], x$end[i])
      # try USCS SQL server first
      query <- try(del_dup_query(location), silent = TRUE)
      if (inherits(query, 'try-error')){
        # try biomaRt next
        query <- try(getBM(c('hgnc_symbol'), filters = c('chromosome_name', 'start', 'end'), values = location, mart = ensembl), silent = TRUE)
        if (inherits(query, 'try-error')){
          error("Check your internet connection. Connection to USCS SQL and biomaRt server failed!")
        }
      }
      query2 <- unlist(query)
      x$genes[i] <- paste(query2, collapse = ",")
      
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
  #' @details The table is written to "CNV_Table.xlsx".
  output <- data.frame()
  genesSplitList <- strsplit(as.character(cnvs$genes), split = ',')
  entryLength <- lapply(genesSplitList, length)
  for(i in 1:dim(cnvs)[1]){
    if (cnvs$genes[i] != "") {
      chr <- rep(cnvs$chr[i], entryLength[[i]])
      start <- rep(cnvs$start[i], entryLength[[i]])
      end <- rep(cnvs$end[i], entryLength[[i]])
      genes <- genesSplitList[[i]]
      genes[2:length(genes)] <- substr(genes[2:length(genes)], start = 2,
                                       stop = nchar(as.character(genes[2:length(genes)])))
      copyNumber <- rep(cnvs$copy.number[i], entryLength[[i]])
      status <- rep(cnvs$status[i], entryLength[[i]])
      tmp <- data.frame(Gene = genes, Chr = chr, Start = start, End = end, CopyNumber = copyNumber, Status = status)
      output <- rbind(output, tmp)
    }
  }
  #write.xlsx(output, file = "CNV_Table.xlsx", keepNA = FALSE, rowNames = FALSE, firstRow = TRUE)
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

cnvs2cbioportal <- function(cnvs, id, outfile_cbioportal){
  #' Export annotated CNVs to a txt file
  #' 
  #' @describtion Export annotated CNVs to a text file for subsequent import in cbioportal
  #' 
  #' @param cnvs dataframe. Table of annotated CNVs
  #' @param id string. Sample ID
  #' @param outfile_cbioportal string. Name of the outputfile
  require(tidyr)
  require(stringi)
  require(org.Hs.eg.db)
  
  cnvs.sub <- subset(cnvs, select = c("copy.number", "genes"))
  cnvs.sub.extended <- data.frame(separate_rows(cnvs.sub,genes,sep=","))
  # remove empty genes
  cnvs.sub.extended <- cnvs.sub.extended[!stri_isempty(cnvs.sub.extended$genes),]
  cnvs.sub.extended$Entrez <- unlist(lapply(mget(cnvs.sub.extended$genes, org.Hs.egSYMBOL2EG, ifnotfound = NA), function(x) x[1]))
  cnvs.out <- data.frame(Hugo_Symbol = cnvs.sub.extended$genes, Entrez_Gene_Id = cnvs.sub.extended$Entrez, Sample_ID = cnvs.sub.extended$copy.number)
  
  # how to deal with multiple copy number variations for a single gene?
  # keep only the maximal CNA
  # cnvs.out <- cnvs.out[!duplicated(cnvs.out$Hugo_Symbol),]
  cnvs.out <- as.data.table(cnvs.out)
  cnvs.out <- cnvs.out[cnvs.out[, .I[which.max(Sample_ID)], by=Hugo_Symbol]$V1]
  
  # set sample specific column name
  colnames(cnvs.out)[3] <- paste(id,"TD",sep = "_")
  # convert total copy numbers from Control-FREEC to allowed values
  # allowed values -2, -1, 0, 1, 2, NA
  cnvs.out$SLGFSK_TD[cnvs.out$SLGFSK_TD == 0] <- -2
  cnvs.out$SLGFSK_TD[cnvs.out$SLGFSK_TD == 1] <- -1
  cnvs.out$SLGFSK_TD[cnvs.out$SLGFSK_TD == 2] <- 0
  cnvs.out$SLGFSK_TD[cnvs.out$SLGFSK_TD == 3] <- 1
  cnvs.out$SLGFSK_TD[cnvs.out$SLGFSK_TD > 3] <- 2
  
  write.table(x = cnvs.out, file = outfile_cbioportal, quote = F, sep = "\t", col.names = T, row.names = F)
}

# TODO
cnv_panel <- function(input_file, outfile, outfile_ts_og, outfile_ideogram, path_data, sureselect, targets_txt, protocol) {
  #' CNV Analysis for Tumor-Only Data
  #'
  #' @description Get Tables for CNVs
  #'
  #' @param input_file string. Filename for Data to be analysed
  #' @param outfile string. Filename for output file
  #' @param outfile_ts_og string. Filename for TSG and OG table
  #' @param outfile_ideogram string. Filename for Plot
  #' @param path_data string. Path to databases
  #' @param mode string. Capture Kit
  #'
  #' @return loss_gain table. Table containing all losses and gains
  #' @return tsg_og table. Table containing CNVs in TSG and OG
  #'
  #' @details CNVs detected in the Pipeline are analysed according to their
  #' @details status. That means if they are losses or gains. Further CNVs
  #' @details for tumorsuppressor genes (TSG) and oncogenes (OG) are listed
  #' @details in another table. If the capture kit is covering most of the
  #' @details genome, there is a plot of the CNVS created.
  cnvs <- read.table(file = input_file, header = TRUE)
  cnvs$CN <- 2
  cnvs$CN <- (2^cnvs$log2)
  cnvs$CN <- round(cnvs$CN)
  cnvs$length <- 0
  cnvs$length <- (cnvs$end - cnvs$start)
  #if (mode == "TruSight_Amplicon"){
  #  genelist <- c("ABL1", "AKT1", "ALK", "APC", "ATM", "BRAF", "CDH1", "CDKN2A",
  #                "CSF1R", "CTNNB1", "EGFR", "ERBB2", "ERBB4", "FBXW7", "FGFR1",
  #                "FGFR2", "FGFR3", "FLT3", "GNA11", "GNAQ", "GNAS", "HNF1A",
  #                "HRAS", "IDH1", "JAK2", "JAK3", "KDR", "KIT", "KRAS", "MET",
  #                "MLH1", "MPL", "NOTCH1", "NPM1", "NRAS", "PDGFRA", "PIK3CA",
  #                "PTEN", "PTPN11", "RB1", "RET", "SMAD4", "SMARCB1", "SMO",
  #                "SRC", "STK11", "TP53", "VHL")
  #} else if (mode == "TruSight_Tumor"){
  #  genelist <- c("AKT1", "BRIP1", "AKT2", "BTK", "AKT3", "CARD11", "ALK",
  #                "CCND1", "APC", "CCND2", "AR", "CCNE1", "ARID1A", "CD79A",
  #                "ATM", "CD79B", "ATR", "CDH1", "BAP1", "CDK12", "BARD1",
  #                "CDK4", "BCL2", "CDK6", "BCL6", "CDKN2A", "BRAF", "CEBPA",
  #                "BRCA1", "CHEK1", "BRCA2", "CHEK2", "CREBBP", "CSF1R",
  #                "CTNNB1", "DDR2", "DNMT3A", "EGFR", "EP300", "ERBB2",
  #                "ERBB3", "ERBB4", "ERCC1", "ERCC2", "ERG", "ESR1", "EZH2",
  #                "FAM175A", "FANCI", "FGFR2", "FANCL", "FGFR3", "FBXW7",
  #                "FGFR4", "FGF1", "FLT1", "FGF2", "FLT3", "FGF3", "FOXL2",
  #                "FGF4", "GEN1", "FGF5", "GNA11", "FGF6", "GNAQ", "FGF7",
  #                "GNAS", "FGF8", "HNF1A", "FGF9", "HRAS", "FGF10", "IDH1",
  #                "FGF14", "IDH2", "FGF23", "INPP4B", "FGFR1", "JAK2", "JAK3",
  #                "MSH3", "KDR", "MSH6", "KIT", "MTOR", "KMT2A", "MLL",
  #                "MUTYH", "KRAS", "MYC", "MAP2K1", "MYCL1", "MAP2K2", "MYCN",
  #                "MCL1", "MYD88", "MDM2", "NBN", "MDM4", "NF1", "MET",
  #                "NOTCH1", "MLH1", "NOTCH2", "MLLT3", "NOTCH3", "MPL", "NPM1",
  #                "MRE11A", "NRAS", "MSH2", "NRG1", "PALB2", "PDGFRA",
  #                "PDGFRB", "PIK3CA", "RAD51D", "TSC1", "RAD54L", "TSC2", "RB1",
  #                "VHL", "RET", "XRCC2", "PIK3CB", "RICTOR", "PIK3CD", "ROS1",
  #                "PIK3CG", "RPS6KB1", "PIK3R1", "SLX4", "PMS2", "SMAD4",
  #                "PPP2R2A", "SMARCB1", "PTCH1", "SMO", "PTEN", "SRC", "PTPN11",
  #                "STK11", "RAD51", "TERT", "RAD51B", "TET2", "RAD51C", "TP53")
  #}
  if (protocol == "panelTumor"){
    geneList <- read.table(file = targets_txt, stringsAsFactors = FALSE)$V1
  }

  cnvs$genename<- lapply(strsplit(x = as.character(cnvs$gene), split = ".",
                                  fixed = T), function(x) x[1])
  cnvs$genename<- lapply(strsplit(x = as.character(cnvs$genename), split = "_",
                                  fixed = T), function(x) x[1])
  i <- which (cnvs$genename %in% genelist)
  cnvs$genename[-i] <- gsub("\\d+$", "", cnvs$genename[-i])
  
  # Find all tumorsuppressors and oncogenes with cnv
  db <- read.delim(paste(path_data, "cancerGeneList.tsv", sep = "/"),
                   header = T, sep = "\t", colClasses = "character")
  cnvs$is_tsg <- 0
  ts <- which(db$OncoKB.TSG == "Yes")
  idt <- which (cnvs$genename %in% db$Hugo.Symbol[ts])
  cnvs_ts <- cnvs[idt, c(1:10)]
  idt <- which( cnvs_ts$CN != 2)
  cnvs_ts <- cnvs_ts[idt, ]

  cnvs$is_og <- 0
  og <- which(db$OncoKB.Oncogene == "Yes")
  ido <- which (cnvs$genename %in% db$Hugo.Symbol[og])
  cnvs_og <- cnvs[ido, c(1:10)]
  ido <- which( cnvs_og$CN != 2)
  cnvs_og <- cnvs_og[ido, ]

  # Find all gains and losses and print them in xlsx table
  idg <- which(cnvs$CN > 2)
  idl <- which(cnvs$CN < 2)
  id <- c(idg, idl)
  cnvs_l <- cnvs[idl, ]
  cnvs_g <- cnvs[idg, ]

  lst1 <- list("CN-Losses" = cnvs_l[, c("chromosome", "start", "end",
                                        "genename", "CN", "length", "is_tsg",
                                        "is_og")],
               "CN-Gains" = cnvs_g[, c("chromosome", "start", "end", "genename",
                                       "CN", "length", "is_tsg", "is_og")])
  write.xlsx(lst1, file = outfile)


  lst2 <- list("CN-TSG" = cnvs_og, "CN-OG" = cnvs_ts)
  write.xlsx(lst2, file = outfile_ts_og)

  if (mode == "TruSight_Tumor"){
    make_cnv_ideo_sig(ratio_file = cnvs, outfile_ideogram = outfile_ideogram, protocol = "Tumor_Only")
  }

  return(list(loss_gain = lst1, tsg_og = lst2))
}
