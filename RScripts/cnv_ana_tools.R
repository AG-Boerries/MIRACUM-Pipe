####################################################
# Assigning Significance to Copy Number Variations #
# ControlFreec
assess_significance <- function(ratio_file, cnvs_file){
  #' Assess Significance
  #'
  #' @description Print Number of Somatic Mutations
  #'
  #' @param ratio_file string. Filename for ratio
  #' @param cnvs_file string. Filename for CNVs
  #'
  #' @details Adds p-value of Wilcoxon Rank Sum Test and p-value of Kolmogorov
  #' @details Smirnov to table of copy number variations.
  dataTable <-read.table(ratio_file, header = TRUE);
  ratio <-data.frame(dataTable)
  dataTable <- read.table(cnvs_file, header = FALSE)
  cnvs <- data.frame(dataTable)
  ratio$Ratio[which(ratio$Ratio == - 1)] = NA
  cnvs.bed = GRanges(cnvs[, 1], IRanges(cnvs[, 2], cnvs[, 3]))  
  ratio.bed = GRanges(ratio$Chromosome, IRanges(ratio$Start, ratio$Start),
                      score=ratio$Ratio)
  overlaps <- subsetByOverlaps(ratio.bed, cnvs.bed)
  normals <- GenomicRanges::setdiff(ratio.bed, cnvs.bed)
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
  return(cnvs)
}

make_cnv_graph <- function(ratio_file, ploidity = '2', outfile_plot){
  #' Make CNV Graph
  #'
  #' @description Plot CNVs on Genome
  #'
  #' @param ratio_file string. Filename of ratio file
  #' @param ploidity numerical. Ploidity (default: 2)
  #' @param outfile_plot string. Filename of output file
  #'
  #' @details This function produces a plot that visualizes all the CNVs on the
  #' @details genome. Therby it is not important, if the CNVs are significant.
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
  cytoband_df = circlize::read.cytoband(species = "hg19")$df
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
  # if (protocol != "panelTumor"){
  #   dataTable <- read.table(ratio_file, header = FALSE, skip = 1)
  #   ratio <- data.frame(dataTable)
  # } else {
  #   ratio <- as.data.frame(ratio_file)
  # }
  ratio <- as.data.frame(ratio_file)

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

del_dup_query <- function(da.fr, server){
  require(RMySQL)
  require(doMC)
  time1 <- Sys.time()
  chroms <- unlist(da.fr[[1]])
  ucscchroms <- unlist(lapply(chroms, function(i) paste0("chr",i)))
  start_list <- unlist(da.fr[[2]])
  stop_list <- unlist(da.fr[[3]])
  out <- mclapply(1:length(ucscchroms), function(i){
    con<- dbConnect(RMySQL::MySQL(),
                    host = server,
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

cnv_annotation <- function(cnv_pvalue_txt, dbfile,
                           path_data, path_script, ucsc_server){
  #' CNV Annotation
  #'
  #' @description Annotate the CNV Regions
  #'
  #' @param cnv_pvalue_txt string. Filename of CNV regions file
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
  x <- as.data.frame(cnv_pvalue_txt)
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

# Try annotation first with SQL server if this fails try with biomaRt

  ensembl <- NULL
  
  for(i in 1:nrow(x)) {
    cat("Processing CNV#", i, "\n")
    
    if (x$WilcoxonRankSumTestPvalue[i] < 0.05
       & x$KolmogorovSmirnovPvalue[i] < 0.05
       & !is.na(x$WilcoxonRankSumTestPvalue[i])
       & !is.na(x$KolmogorovSmirnovPvalue[i])) {
      location <- list(x$chr[i], x$start[i], x$end[i])
      # try USCS SQL server first
      query <- try(del_dup_query(location, ucsc_server), silent = TRUE)
      if (inherits(query, 'try-error')){
        # try biomaRt next
        if (is.null(ensembl)) {
          ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl", host="grch37.ensembl.org")
        }
        query <- try(getBM(c('hgnc_symbol'), filters = c('chromosome_name', 'start', 'end'), values = location, mart = ensembl), silent = TRUE)
        if (inherits(query, 'try-error')){
          error("Check your internet connection. Connection to UCSC SQL and biomaRt server failed!")
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

  onco <- x1[x1$Oncogene != ".", c("chr", "copy.number", "status", "Oncogene")]
  tumorSupp <- x1[x1$TumorSuppressor != ".", c("chr", "copy.number", "status", "TumorSuppressor")]
  
  return(list(CNVsAnnotated = x1, CNVOncogenes = onco, CNVTumorSuppressors = tumorSupp))
}


cnv_processing <- function(cnv_file, targets, db,
                           path_data, path_script){
  #' CNV Processing
  #'
  #' @description Copy Number Analysis
  #'
  #' @param cnv_file string. Filename of genes with a CNV
  #' @param targets string. Filename of targeted genes
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
  cnv <- as.data.frame(cnv_file)
  
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
  loss_go <- get_terms(dataset = db, mut.entrez = prep$de_genes, t2.entrez = prep$universe, outfile = NULL)
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
  gain_go <- get_terms(dataset = db, mut.entrez = prep$de_genes, t2.entrez = prep$universe, outfile = NULL)
  return(list(loss_go = loss_go, gain_go = gain_go))
}

cnv_table <- function(cnvs, ampl_genes){
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
      # if(length(genes) > 1) {
      #   genes[2:length(genes)] <- substr(genes[2:length(genes)], start = 2,
      #                                    stop = nchar(as.character(genes[2:length(genes)])))
      # }
      copyNumber <- rep(cnvs$copy.number[i], entryLength[[i]])
      status <- rep(cnvs$status[i], entryLength[[i]])
      tmp <- data.frame(Gene = genes, Chr = chr, Start = start, End = end, CopyNumber = copyNumber, Status = status)
      output <- rbind(output, tmp)
    }
  }
  if (!is.null(ampl_genes)){
    id <- which(output$Gene %in% ampl_genes)
    output <- output[id, ]
  }
  #write.xlsx(output, file = "CNV_Table.xlsx", keepNA = FALSE, rowNames = FALSE,firstRow = TRUE)
  return(output)
}

# check cnvs for dna damage response
cnv_pathways <- function(input, db) {
  #' CNV DNA Damage Response
  #'
  #' @description Get DNA Damage Response genes with CNVs
  #'
  #' @param input dataframe. Table of genes with CNVs
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
  gene <- as.character(input$Gene[idx])
  copynumber <- input$CopyNumber[idx]
  status <- as.character(input$Status[idx])
  l <- cbind(gene, copynumber, status)
  return(l)
}

cnvs2cbioportal <- function(cnvs, id, outfile_cbioportal, gender, ampl_genes = NULL){
  #print(gender)
  #print(ampl_genes)
  #' Export annotated CNVs to a txt file
  #'
  #' @describtion Export annotated CNVs to a text file for subsequent import in cbioportal
  #'
  #' @param cnvs dataframe. Table of annotated CNVs
  #' @param id string. Sample ID
  #' @param outfile_cbioportal string. Name of the outputfile
  require(tidyr)
  require(stringi)

  # due to the conversion to -2:2 (GISTIC) sex chromosomes by a male individuum need adjustments because these are only present as a single copy
  if (gender == "XY") {
    ids <- which(cnvs$Chr %in% c("X", "Y") & cnvs$CopyNumber >= 1)
    cnvs[ids, "CopyNumber"] <- cnvs[ids, "CopyNumber"] + 1
  }

  cnvs.sub <- subset(cnvs, select = c("CopyNumber", "Gene"))
  cnvs.sub.extended <- na.omit(cnvs.sub)
  # remove empty genes
  cnvs.sub.extended <- cnvs.sub.extended[!stri_isempty(cnvs.sub.extended$Gene),]
  cnvs.sub.extended$Entrez <- unlist(lapply(mget(as.character(cnvs.sub.extended$Gene), org.Hs.egSYMBOL2EG, ifnotfound = NA), function(x) x[1]))
  cnvs.out <- data.frame(Hugo_Symbol = cnvs.sub.extended$Gene, Entrez_Gene_Id = cnvs.sub.extended$Entrez, Sample_ID = cnvs.sub.extended$CopyNumber)
 
  # how to deal with multiple copy number variations for a single gene?
  # keep only the maximal CNA
  # cnvs.out <- cnvs.out[!duplicated(cnvs.out$Hugo_Symbol),]
  cnvs.out <- data.table::as.data.table(cnvs.out)
  cnvs.out <- cnvs.out[cnvs.out[, .I[which.max(Sample_ID)], by=Hugo_Symbol]$V1]

  # for TSO500 report all analyzed genes; all genes not found in cnvs.out are assumed to have no change in CN and get therefore a 2 (diploid)
  if(!is.null(ampl_genes)) {
    ampl_genes.entrez <- unlist(lapply(mget(as.character(ampl_genes), org.Hs.egSYMBOL2EG, ifnotfound = NA), function(x) x[1]))
    cnvs.complete <- data.frame(Hugo_Symbol = ampl_genes, Entrez_Gene_Id = ampl_genes.entrez, Sample_ID = 2)

    cnvs.complete[na.omit(match(cnvs.out$Hugo_Symbol, cnvs.complete$Hugo_Symbol)), "Sample_ID"] <- cnvs.out$Sample_ID
    cnvs.out <- cnvs.complete
  }

  # set sample specific column name
  colnames(cnvs.out)[3] <- paste(id,"TD",sep = "_")
  # convert total copy numbers from Control-FREEC to allowed values
  # allowed values -2, -1, 0, 1, 2, NA
  cnvs.out[,3][cnvs.out[,3] == 0] <- -2
  cnvs.out[,3][cnvs.out[,3] == 1] <- -1
  cnvs.out[,3][cnvs.out[,3] == 2] <- 0
  cnvs.out[,3][cnvs.out[,3] <= 6 & cnvs.out[,3] >= 3 ] <- 1
  cnvs.out[,3][cnvs.out[,3] > 6] <- 2

  write.table(x = cnvs.out, file = outfile_cbioportal, quote = F, sep = "\t", col.names = T, row.names = F)
}


get_type <- function(Oncogenes, Tumorsuppressor, CNVsAnnotated, sureselect_type) {
  require(EnsDb.Hsapiens.v75)
  require(GenomicRanges)
  
  if(sureselect_type != "TSO500"){
    genes_onc <- Oncogenes$Oncogene
    genes_onc <- unique(unlist(strsplit(genes_onc, split = ",")))
    genes_tsg <- Tumorsuppressor$TumorSuppressor
    genes_tsg <- unique(unlist(strsplit(genes_tsg, split = ",")))
    
    edb <- EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
    genes_edb <- genes(edb)
    
    gene_loci_tsg <- genes_edb[which(genes_edb$symbol %in% genes_tsg), ]
    gene_loci_tsg <- gene_loci_tsg[which(gene_loci_tsg$gene_biotype == "protein_coding"), ]
    gene_loci_tsg <- gene_loci_tsg[which(seqnames(gene_loci_tsg) %in% c(1:22, "X"))]
    
    gene_loci_onc <- genes_edb[which(genes_edb$symbol %in% genes_onc), ]
    gene_loci_onc <- gene_loci_onc[which(gene_loci_onc$gene_biotype == "protein_coding"), ]
    gene_loci_onc <- gene_loci_onc[which(seqnames(gene_loci_onc) %in% c(1:22, "X"))]
    
    cnvs <- data.frame(chr = CNVsAnnotated$chr,
                      start = CNVsAnnotated$start,
                      end = CNVsAnnotated$end,
                      cn = CNVsAnnotated$copy.number)
    gene_loci_onc <- as.data.frame(gene_loci_onc)[, c(1:3, 7)]
    gene_loci_tsg <- as.data.frame(gene_loci_tsg)[, c(1:3, 7)]
    
    get_cn <- function(x){
      cn <- cnvs$cn[which(x$start >= cnvs$start
                          & x$end <= cnvs$end
                          & x$seqnames == cnvs$chr)]
      return(cn)
    }
    gene_loci_onc$cn <- "-"
    gene_loci_onc$Type <- "focal"
    for (i in 1:dim(gene_loci_onc)[1]) {
      if(length(get_cn(gene_loci_onc[i, ])) > 0) {
        gene_loci_onc$cn[i] <- get_cn(gene_loci_onc[i, ])
        gene_loci_onc$Type[i] <- "complete"
      } else {
        gene_loci_onc$cn[i] <- paste0(Oncogenes$copy.number[grep(gene_loci_onc$gene_name[i],
                                                                x = Oncogenes$Oncogene)],
                                      sep = "", collapse = ", ")
      }
    }
    
    gene_loci_tsg$cn <- "-"
    gene_loci_tsg$Type <- "focal"
    for (i in 1:dim(gene_loci_tsg)[1]) {
      if(length(get_cn(gene_loci_tsg[i, ])) > 0) {
        gene_loci_tsg$cn[i] <- get_cn(gene_loci_tsg[i, ])
        gene_loci_tsg$Type[i] <- "complete"
      } else {
        gene_loci_tsg$cn[i] <- paste0(Tumorsuppressor$copy.number[grep(gene_loci_tsg$gene_name[i],
                                                                      x = Tumorsuppressor$TumorSuppressor)],
                                      sep = "", collapse = ", ")
      }
    }

    return(list(gene_loci_onc = gene_loci_onc, gene_loci_tsg = gene_loci_tsg))

  } else if(sureselect_type == "TSO500"){
    genes_onc <- as.character(Oncogenes$Gene)
    genes_onc <- unique(unlist(strsplit(genes_onc, split = ",")))
    
    edb <- EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
    genes_edb <- genes(edb)
    
    gene_loci <- genes_edb[which(genes_edb$symbol %in% genes_onc), ]
    gene_loci <- gene_loci[which(gene_loci$gene_biotype == "protein_coding"), ]
    gene_loci <- gene_loci[which(seqnames(gene_loci) %in% c(1:22, "X", "Y"))]
    
    cnvs <- data.frame(chr = Oncogenes$Chr,
                       start = Oncogenes$Start,
                       end = Oncogenes$End,
                       cn = Oncogenes$CopyNumber)
    gene_loci <- as.data.frame(gene_loci)[, c(1:3, 7)]
    
    get_cn <- function(x){
      cn <- cnvs$cn[which(x$start >= cnvs$start
                          & x$end <= cnvs$end
                          & x$seqnames == cnvs$chr)]
      return(cn)
    }
    gene_loci$cn <- "-"
    gene_loci$Type <- "focal"
    for (i in 1:dim(gene_loci)[1]) {
      if(length(get_cn(gene_loci[i, ])) > 0) {
        gene_loci$cn[i] <- get_cn(gene_loci[i, ])
        gene_loci$Type[i] <- "complete"
      } else {
        gene_loci$cn[i] <- paste0(Oncogenes$CopyNumber[grep(gene_loci$gene_name[i],
                                                                x = Oncogenes$Gene)],
                                      sep = "", collapse = ", ")
      }
    }
    
    if(any(unlist(lapply(strsplit(gene_loci$cn, split = ", "), length)) > 1)){
      ids <- which(unlist(lapply(strsplit(gene_loci$cn, split = ", "), length)) > 1)
      singleCNs <- strsplit(gene_loci$cn, split = ", ")
      gene_loci$cn <- unlist(lapply(strsplit(gene_loci$cn, split = ", "), function(c) c[1]))
      gene_loci2 <- gene_loci
      for (i in 1:length(ids)){
        dupRow <- gene_loci[ids[i],]
        dupRow$cn <- singleCNs[[ids[i]]][2]
        gene_loci2 <- rbind(gene_loci2, dupRow)
      }
      gene_loci <- gene_loci2
    }
    
    return(list(gene_loci = gene_loci))
  }

}
hrd_extr <- function(hrd_file) {
  hrd <- read.delim(hrd_file, sep = " ")
  sum <- hrd$HRD.sum
  score <- paste0(sum, " (", hrd$HRD, "|", hrd$Telomeric.AI, "|", hrd$LST, ")")
  return(list(sum = sum, score = score))
}

purity_extr <- function(purity_file) {
  alt_solution <- read.delim(purity_file, header = TRUE)
  puri <- alt_solution[1, 1]
  ploi <- alt_solution[1, 2]
  return(list(purity = puri, ploidy = ploi))
}

freec2seg <- function(cnvs_file, cpn_file, id, outfile_seg) {
  #' Export annotated CNVs to a seg file
  #'
  #' @describtion Export annotated CNVs to a seg file compliant to the GISTIC syntax
  #'
  #' @param cnvs_file string. Filename for CNVs
  #' @param cpn_file string. Filename for CPN
  #' @param id string. Sample ID
  #' @param outfile_seg string. Name of the outputfile

  bam <- read.table(cnvs_file, header = F)
  cpn <- read.table(cpn_file, header = F)
  bam$V4 <- log2(bam$V4) - 1
  bam$V5 <- paste(id, "TD", sep = "_")
  bam$V6 <- 0
  bam <- bam[c(5,1:3,6,4)]
  colnames(bam) <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")

  for (i in 1:nrow(bam)) {
    entry <- bam[i,]
    bam[i,]$num.mark <- sum(cpn[cpn$V1 == entry$chrom & (cpn$V2 >= entry$loc.end & cpn$V2 <= entry$loc.end | cpn$V3 >= entry$loc.start & cpn$V3 <= entry$loc.end),]$V4)
  }

  bam$seg.mean[!is.finite(bam$seg.mean)] <- -1

  write.table(bam, outfile_seg, sep = "\t", quote = F, row.names = F)
}
