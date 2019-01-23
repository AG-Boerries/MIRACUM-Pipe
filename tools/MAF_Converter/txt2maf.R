txt2maf <- function(input, Center = 'Freiburg', refBuild = 'hg19', idCol = NULL, id = NULL, sep = "\t", Mutation_Status = c("Somatic", "Germline","LoH")[1]){
	
  require("data.table")
  require("org.Hs.eg.db")
  require("gtools")
  
	ann <- read.table(input, header = T, sep = sep, stringsAsFactors = F)
	
	essential.col = c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "ExonicFunc.refGene",
	                  "AAChange.SnpEff", "CChange.SnpEff", "Transcript.SnpEff", "Ensembl.SnpEff", "avsnp150")
	
	 for(i in 1:length(essential.col)){
	   colId = suppressWarnings(grep(pattern = paste0("^",essential.col[i], "$"),
		 x = colnames(ann), ignore.case = TRUE))
		 if (length(colId) == 1) {
		   colnames(ann)[colId] = essential.col[i]
		 }
	 }
	 if(length(essential.col[!essential.col %in% colnames(ann)]) > 0){
	   message("Available fields:")
		 print(colnames(ann))
		 message(paste0("Missing required field in input file: "))
		 print(essential.col[!essential.col %in% colnames(ann)])
		 stop()
	 }
	 if(is.null(idCol) & is.null(id)) {
	   error('Provide either the column containing the Tumor Sample Barcode or the Tumor Sample Barcode as string!')
		 stop()
	 }
	 if(is.null(idCol) & !is.null(id)) {
	   ann$Tumor_Sample_Barcode <- as.character(id)
	 }
	 
	 if(!is.null(idCol)) {
	   colnames(ann)[which(colnames(ann) == idCol)] = "Tumor_Sample_Barcode"
	 }
     
	 if(is.null(Center)) {
       Center = NA
	 }
	
	 ann$uid = paste("uid", 1:nrow(ann), sep = "")
   ann.mand = c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "ExonicFunc.refGene",
                "AAChange.SnpEff", "CChange.SnpEff", "Transcript.SnpEff", "Ensembl.SnpEff", "Tumor_Sample_Barcode", "avsnp150", "uid")
				  
  ann.opt = colnames(ann)[!colnames(ann) %in% ann.mand]
  ann.opt = c(ann.opt, "uid")
  ann.opt = ann[, ann.opt]
  ann = ann[, ann.mand]
  ann$ExonicFunc.refGene = gsub(pattern = " SNV", replacement = "", x = ann$ExonicFunc.refGene)
  funcSpl = strsplit(x = as.character(ann$ExonicFunc.refGene), split = ";", fixed = TRUE)
  funcSpl = sapply(funcSpl, function(l) {l[length(l)]})
	ann$ExonicFunc.refGene = funcSpl
     
	funcRef = strsplit(x = as.character(ann$Func.refGene), split = ";", fixed = TRUE)
  funcRef = sapply(funcRef, function(l) {l[length(l)]})
  ann$Func.refGene = funcRef
	 
  ann$ExonicFunc.refGene = ifelse(test = ann$Func.refGene == "intronic", yes = "Intron", no = ann$ExonicFunc.refGene)
  ann$ExonicFunc.refGene = ifelse(test = ann$Func.refGene == "intergenic", yes = "IGR", no = ann$ExonicFunc.refGene)
  ann$ExonicFunc.refGene = ifelse(test = ann$Func.refGene == "downstream", yes = "3'Flank", no = ann$ExonicFunc.refGene)
  ann$ExonicFunc.refGene = ifelse(test = ann$Func.refGene == "upstream", yes = "5'Flank", no = ann$ExonicFunc.refGene)
  ann$ExonicFunc.refGene = ifelse(test = ann$Func.refGene == "splicing", yes = "Splice_Site", no = ann$ExonicFunc.refGene)
  ann$ExonicFunc.refGene = ifelse(test = ann$Func.refGene == "UTR3", yes = "3'UTR", no = ann$ExonicFunc.refGene)
  ann$ExonicFunc.refGene = ifelse(test = ann$Func.refGene == "UTR5", yes = "5'UTR", no = ann$ExonicFunc.refGene)
  ann$ExonicFunc.refGene = ifelse(test = ann$Func.refGene %in% c("ncRNA_exonic", "ncRNA_intronic", "ncRNA_UTR3", "ncRNA_UTR5", "ncRNA"), yes = "RNA", no = ann$ExonicFunc.refGene)
  ann.lvls = c("synonymous", "nonsynonymous", "stopgain", "stoploss", "frameshift insertion", "frameshift deletion", "nonframeshift insertion", "nonframeshift deletion",
               "Intron", "IGR", "Splice_Site", "3'UTR", "3'Flank", "5'UTR", "5'Flank", "unknown", "UNKNOWN", "RNA")
  ann.lbls = c("Silent", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del", "In_Frame_Ins", "In_Frame_Del",
               "Intron", "IGR", "Splice_Site", "3'UTR", "3'Flank", "5'UTR", "5'Flank", "UNKNOWN", "UNKNOWN", "RNA")
  names(ann.lbls) = ann.lvls
	ann$ExonicFunc.refGene = as.character(ann.lbls[as.character(ann$ExonicFunc.refGene)])
     
	ann.del = ann[ann$Alt %in% "-",]
	ann = ann[!ann$Alt %in% "-",]
	
	if(nrow(ann.del) > 0){
	  ann.del$var.type = "DEL"
    } else {
      ann.del$var.type = rep(NA, times = dim(ann.del)[1])
    }
	 
  ann.ins = ann[ann$Ref %in% "-",]
  ann = ann[!ann$Ref %in% "-",]
  if(nrow(ann.ins) > 0){
    ann.ins$var.type = "INS"
  } else {
      ann.ins$var.type =  rep(NA, times = dim(ann.ins)[1])
  }

	 
  if(nrow(ann) > 0){
    ann$var.type = "SNP"
  } else {
      ann$var.type =  rep(NA, times = dim(ann)[1])
  }
	 
  ann = rbind(ann, ann.del, ann.ins)
  ann.splice = ann[ann$ExonicFunc.refGene == "Splice_Site",]
  
  if(nrow(ann.splice) > 0){
    ann = ann[ann$ExonicFunc.refGene != "Splice_Site",]
    ann.splice$Gene.refGene = sapply(strsplit(x = as.character(ann.splice$Gene.refGene), split = "(", fixed = TRUE), "[[", 1)
    ann = rbind(ann, ann.splice)
  }
	 
	symbol = unlist(lapply(strsplit(ann$Gene.refGene, split = ";"), function(x) {x[1]}))
    idx <- which(!is.na(symbol))
    entrez <- rep(NA, times = length(symbol))
	entrez = mget(as.character(symbol), org.Hs.egSYMBOL2EG, ifnotfound = NA)
	 
	aa = unlist(lapply(strsplit(x = as.character(ann$AAChange.SnpEff), split = ";", fixed = T), function(x) x[1]))
	aa_short = c("H", "Q", "P", "R", "L", "D", "E", "A", "G", "V", "Y", "S", "C", "W", "F", "N", "K", "T", "I", "M", "fs", "X")
	aa_long = c("His", "Gln", "Pro", "Arg", "Leu", "Asp", "Glu", "Ala", "Gly", "Val", "Tyr", "Ser", "Cys", "Trp", "Phe", "Asn", "Lys", "Thr", "Ile", "Met", "fs", "X")
	names(aa_short) <- aa_long
	aa = gsub(aa, pattern = '*',replacement = 'X',fixed = T)
	aa.num = as.numeric(gsub("[^\\d]+", "", aa, perl=TRUE))
	aa = unlist(lapply(strsplit(aa , split = '.', fixed = T), function(s) s[2]))
	aa.split = strsplit(aa, split = "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])", perl=T)
	aa.split = lapply(aa.split, function(c) {aa_short[c]})
	aa.short = do.call(rbind,aa.split)
    if(length(which(is.na(aa.num))) != length(aa.num)) {
        aa.short[, 2] = aa.num
        proteinChange = paste0("p.", aa.short[, 1], aa.short[, 2], aa.short[, 3])
    } else {
        aa.short =  NA
        proteinChange = NA
    }
	proteinChange[proteinChange == "p.NANANA"] = ""
	Transcript_Id = ann$Transcript.SnpEff
	Transcript_Id[is.na(Transcript_Id)] <- ""
	TxChange = unlist(lapply(strsplit(x = as.character(ann$CChange.SnpEff), split = ";", fixed = T), function(x) x[1]))
	ensembl = ann$Ensembl.SnpEff
	ensembl[is.na(ensembl)] <- ""
	 
  ann.maf = data.table::data.table(Hugo_Symbol = ann$Gene.refGene,
                                   Entrez_Gene_Id = as.character(entrez),
                                   Center = Center,
                                   NCBI_Build = refBuild,
                                   Chromosome = ann$Chr,
                                   Start_Position = ann$Start,
                                   End_Position = ann$End,
                                   Strand = "+",
                                   Variant_Classification = ann$ExonicFunc.refGene,
                                   Variant_Type = ann$var.type,
                                   Reference_Allele = ann$Ref,
                                   Tumor_Seq_Allele1 = ann$Ref,
                                   Tumor_Seq_Allele2 = ann$Alt,
                                   dbSNP_RS = ann$avsnp150,
                                   Tumor_Sample_Barcode = ann$Tumor_Sample_Barcode,
                                   Mutation_Status = Mutation_Status,
                                   HGSVp_Short = proteinChange,
                                   Amino_Acid_Change = proteinChange,
                                   TxChange = TxChange,
                                   Transcript_Id = Transcript_Id,
                                   ENSEMBL_Gene_Id = ensembl,
                                   uid = ann$uid)
									  
  ann.maf = merge(ann.maf, ann.opt, by = "uid")
  ann.maf = ann.maf[, `:=`(uid, NULL)]
  
  ann.maf <- ann.maf[order(ann.maf[,"Chromosome"], ann.maf[,"Start_Position"]),]
  ann.maf <- ann.maf[mixedorder(ann.maf$Chromosome),]
  
  return(ann.maf)
}
