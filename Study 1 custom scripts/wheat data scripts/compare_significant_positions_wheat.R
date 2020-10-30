#wheat significant SNPs
sig_SNP_file <- read.csv(file = "wheat_significant_marker_names_08012020.csv", header = TRUE, stringsAsFactors = FALSE)

library(GenomicFeatures)
species_list <-loadTaxonomyDb()
which(species_list$species == "aesitvum")

Msi_TxDb <-makeTxDbFromGFF("Taestivum_296_v2.2.gene.gff3",
                           format = "gff3",
                           organism = "Triticum aestivum",
                           dataSource = "Phytozome 12")
Msi_TxDb_gene_only <- genes(Msi_TxDb)

#candidate_genes <- function(sig_genes) {
sig_SNP_file$CHR.TASSEL <- as.character(substring(sig_SNP_file$TASSEL , 11, 13))
sig_SNP_file$CHR.TASSEL <- as.character(gsub("_", "", sig_SNP_file$CHR.TASSEL))
sig_SNP_file$Position.TASSEL <- substring(sig_SNP_file$TASSEL , 19, 28)
sig_SNP_file$Position.TASSEL <- gsub("f_", "", sig_SNP_file$Position.TASSEL)
sig_SNP_file$Position.TASSEL <- gsub("_", "", sig_SNP_file$Position.TASSEL)
markers <- data.frame(Name = sig_SNP_file$TASSEL,
                        Chromosome = sig_SNP_file$CHR.TASSEL,
                        Position = as.integer(sig_SNP_file$Position.TASSEL))
markers$Position <- as.integer(markers$Position)


sig_SNP_file$CHR.polyRAD <- as.character(substring(sig_SNP_file$polyRAD , 7, 9))
sig_SNP_file$Position.polyRAD <- as.numeric(substring(sig_SNP_file$polyRAD, 14, 20))
#sig_SNP_file$Position.polyRAD <- gsub("f_", "", sig_SNP_file$Position.TASSEL)
#sig_SNP_file$Position.polyRAD <- gsub("_", "", sig_SNP_file$Position.TASSEL)
#ensure your chromosome names match how they are documented in the gff file
#markers$Chromosome <- paste("Chr", formatC(markers$Chromosome, flag="0", width=2), sep = "")

markers <- data.frame(Name = sig_SNP_file$polyRAD,
                      Chromosome = sig_SNP_file$CHR.polyRAD,
                      Position = as.integer(sig_SNP_file$Position.polyRAD))
markers$Position <- as.integer(markers$Position)
markers <- markers[-c(7:11),]
markers$Name <- gsub("-[0-9].*$", "", markers$Name)
markers$Name <- as.character(markers$Name)
#make a granges object
markers <- as.data.frame(markers)
#lowercase
markers$Name <- tolower(markers$Name)
markers_granges <- GRanges(Rle(markers$Chromosome, rep(1, length(markers$Chromosome))),
                             IRanges(start = markers$Position - 1000, end = markers$Position + 1000,
                                     # range 1000bp up and downstream of SNP
                                     names = markers$Name))
  
  # range 500bp up and downstream of SNP
  #markers_granges$SNP <- names(markers_granges)
  
  #find genes in the granges
markers_genes <- transcriptsByOverlaps(Msi_TxDb, markers_granges, columns = c("tx_id", "tx_name"))
markers_genes <- as.data.frame(markers_genes)
  
  gene_distance <- distanceToNearest(markers_granges, Msi_TxDb_gene_only)
  distance_col <- as.data.frame(t(rbind(markers_granges,mcols(gene_distance)$distance)))
  markers_granges$SNP <- names(markers_granges)
  distance_col<- as.data.frame(t(rbind(markers_granges$SNP,gene_distance@elementMetadata@listData[["distance"]])))
  #polyRAD_2 = full_genes
  full_genes <- as.data.frame(cbind(markers$Name, distance_col))
  
  #located_in_gene <- which(distance_UNEAK$V2 == 0 )
  #67,789 genes annotated in the msinesis reference genome
  #49,565 of 86,580 (57.25 %) located in gene in genes
  #586,799 out of 1,024,980 (57.25 %)located in genes 
  #UNEAK 18,269 out of the 29,248 sequences that aligned were located in genes
  #29,248 tags were aligned to chromosomes
}