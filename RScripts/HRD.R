args <- commandArgs()
ID <- args[6]
path <- args[7]

library(scarHRD)
library(sequenza)
input_file <- paste0(path, "/", ID, ".small.seqz.gz")
# Determine the HRD-Score
result <- scar_score(input_file, reference = "grch37", seqz = TRUE)
write.table(x = result, file = paste0(path, "/", ID, "_HRD.txt"))
# CNV analysis of SEQUENZA to determine tumor purity and ploidy
extr <- sequenza.extract(input_file, verbose = FALSE)
cp <- sequenza.fit(extr)
sequenza.results(
    sequenza.extract = extr, cp.table = cp, sample.id = ID,
    out.dir = paste0(path, "/", ID, "_sequenza")
    )
