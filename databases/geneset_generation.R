library(GSA)
gmt <- GSA.read.gmt('h.all.v7.0.entrez.gmt')
genesets <- gmt$genesets
names <- data.frame(Names = gmt$geneset.names, Descriptions = gmt$geneset.descriptions)
names(genesets) <- names$Names
hallmarksOfCancer <- genesets
save(hallmarksOfCancer, file = "hallmarksOfCancer_GeneSets.RData")