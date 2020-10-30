#create markers_db dataframe 
library(GenomicFeatures)
#import Miscanthus gff3 file 
Msi_TxDb <-makeTxDbFromGFF("Msinensis_497_v7.1.gene.gff3",
                           format = "gff3",
                           organism = "Miscanthus sinensis",
                           dataSource = "Phytozome 12")
Msi_TxDb_gene_only <- genes(Msi_TxDb)



candidate_genes <- function(input) {
markers <- data.frame(Name = markers$SNP,
                      Chromosome = as.integer(substring(input$SNP , 4, 5)),
                      Position = as.integer(sapply(strsplit(as.character(input$SNP ), "-"), function(x) x[2])))
remove <- which(is.na(markers$Chromosome) == TRUE)
markers <- markers[-c(remove),]

#ensure your chromosome names match how they are documented in the gff file
markers$Chromosome <- paste("Chr", formatC(markers$Chromosome, flag="0", width=2), sep = "")

#make a granges object
markers_granges <- GRanges(Rle(markers$Chromosome, rep(1, length(markers$Chromosome))),
                           IRanges(start = markers$Position, end = markers$Position,
                                   # range 500bp up and downstream of SNP
                                   names = markers$SNP))
markers_granges <- GRanges(Rle(markers$Chromosome, rep(1, length(markers$Chromosome))),
        IRanges(start = markers$Position - 1000, end = markers$Position + 1000, # range 1000bp up and downstream of SNP
                names = markers$Name))


#markers_granges$SNP <- names(markers_granges)

#find genes in the granges
markers_genes <- transcriptsByOverlaps(Msi_TxDb, markers_granges, columns = c("tx_id", "tx_name"))
markers_genes <- as.data.frame(markers_genes)
}

# markers_db <-candidate_genes(GWAS$SNP)
gene_distance <- distanceToNearest(markers_granges, Msi_TxDb_gene_only)
distance_col <- as.data.frame(t(rbind(markers_granges,mcols(gene_distance)$distance)))
markers_granges$SNP <- names(markers_granges)
distance_col<- as.data.frame(t(rbind(markers_granges$SNP,gene_distance@elementMetadata@listData[["distance"]])))
#polyRAD_2 = full_genes
full_genes <- as.data.frame(cbind(markers$Name, distance_col))






#distance_polyRAD <- distanceToNearest()
write.csv(polyRAD_2,file =  "polyRAD_distances_full.csv")
#names(full_genes) <- c("V1", "V2")
full_genes$V1 <- as.numeric(full_genes$V1)
in_gene <- length(which(full_genes$V1 == 0))
#polyRAD_2$V2 <- polyRAD_2$V2 -1
orange_5 <- length(which(full_genes$V1 < 5000))- in_gene
under_30 <- length(which(full_genes$V1 < 30000)) - sum(orange_5, in_gene)
over_30 <- length(which(full_genes$V2 > 30000))
