#download all 3 files
UNEAK_1 <- read.csv(file = "GAPIT..ZJU.Total.culms.year.3.csv", stringsAsFactors = FALSE, header = TRUE)
TASSEL_1 <- read.csv(file = "GAPIT.MLMM.Yld.ZJU.year.3.GWAS.Results.csv", stringsAsFactors = FALSE, header = TRUE)
polyRAD_1 <- read.csv(file = "polyRAD.GAPIT.MLMM.Yld.ZJU.year.3.GWAS.Results.csv", stringsAsFactors = FALSE, header = TRUE)
TASSEL_1 <- TASSEL_1[,c(1)]
remove_tag <- startsWith(as.character(TASSEL_1$SNP), "S", trim=FALSE, ignore.case=FALSE)

#read in UNEAK tags
UNEAK_sig <- read.csv(file = "marker_match_name.csv", header = TRUE)

library(GenomicFeatures)
Msi_TxDb <-makeTxDbFromGFF("Msinensis_497_v7.1.gene.gff3",
                           format = "gff3",
                           organism = "Miscanthus sinensis",
                           dataSource = "Phytozome 12")
Msi_TxDb_gene_only <- genes(Msi_TxDb)

#sig_genes is the name of the SNPS of interest
#the below helps to create a data frame with chromosome and position of each SNP
#make function to take in significant dataframe and take only the names
#UNEAK tag sequences that were identified and aligned to the reference genome

UNEAK_tag_sequences2$SNP <- UNEAK_tag_sequences2$match.name
UNEAK_tag_sequences3 <- UNEAK_tag_sequences2[-c(15),]
UNEAK_tag_sequences3$new.position <- substring(UNEAK_tag_sequences3$SNP, 7, 15 )
UNEAK_tag_sequences3$Position <- as.numeric(UNEAK_tag_sequences3$new.position)

#match UNEAK tag names to names UNEAK tag file
match_tags <- match(UNEAK_tag_file$SNP , match_all$Query)
NA_remove <- which(is.na(match_tags)) 
match_tags <- match_tags[-c(NA_remove)]

tag_names <- as.data.frame(match_all[c(match_tags),])
names(tag_names) <- "SNP"
#remove the ones that are not aligned to actual position
scaffold_remove <- which(grepl("^s", tag_names$Marker.name) == TRUE)
tag_names <- tag_names[-c(scaffold_remove),]
#tag_names <- as.data.frame(tag_names)
#names(tag_names) <- "SNP"

#add genes to polyRAD subset
add_polyRAD <- read.csv("polyRAD.GAPIT.MLMM.CmL.HU.NEF.UI.CSU.KNU.ZJU.year.3.GWAS.Results.csv", header = TRUE, stringsAsFactors = FALSE)

add_polyRAD$new.tag.names <- sub("^S", "Chr", add_polyRAD$SNP)
add_polyRAD$new.tag.names <- sub("\\_[ATCG].*$", "", add_polyRAD$new.tag.names)
add_polyRAD$new.position <- sprintf("%09d", add_polyRAD$Position)
add_polyRAD$new.tag.names <- sub("\\_", "\\-", add_polyRAD$new.tag.names)
add_polyRAD$new.tag.names <- sub("\\-.*$","", add_polyRAD$new.tag.names)
add_polyRAD$new.tag.names2 <- paste(add_polyRAD$new.tag.names, add_polyRAD$new.position, sep = "-")

add_polyRAD <- add_polyRAD[c(1,2,3),]
polyRAD_tag_sequences2 <- rbind(add_polyRAD, polyRAD_tag_sequences)
add_polyRAD$new.name <- add_polyRAD$new.tag.names2
add_polyRAD$upper.position <- add_polyRAD$Position + 1000
add_polyRAD$lower.position <- add_polyRAD$Position - 1000

add_polyRAD2 <- add_polyRAD[,-c(11,12,13)]
polyRAD_tag_sequences2 <- rbind(add_polyRAD, polyRAD_tag_sequences)

TASSEL_GWAS_tags <- read.csv(file = "TASSEL.GAPIT.MLMM.CmL.HU.NEF.UI.CSU.KNU.ZJU.year.3.GWAS.Results.csv", header = TRUE, stringsAsFactors = FALSE)
UNEAK_tags <- read.csv(file = "marker_match_name.csv", header = TRUE, stringsAsFactors = FALSE)
TASSEL_1 <- as.data.frame(TASSEL_1)
#names(UNEAK_tags)[1] <- "SNP" 
#UNEAK_tags <- as.data.frame(UNEAK_tags$SNP)
#UNEAK_tags <- UNEAK_tags[-c(4,18,54,55),]
UNEAK_tags <- as.data.frame(UNEAK_tags)
polyRAD_sig <- read.csv(file = "polyRAD_sigs_03112020.csv", header = TRUE, stringsAsFactors = FALSE)
polyRAD_sig <- as.data.frame(polyRAD_1$SNP)

candidate_genes <- function(sig_genes) {
markers <- data.frame(Name = UNEAK_tags$SNP,
                      Chromosome = as.integer(substring(UNEAK_tags$SNP , 4, 5)),
                      Position = as.integer(sapply(strsplit(as.character(UNEAK_tags$SNP ), "-"), function(x) x[2])))
remove <- which(is.na(markers$Chromosome) == TRUE)
markers <- markers[-c(remove),]
#remove <- which(is.na(markers$Chromosome)==TRUE)
#markers <- markers[-c(remove),]
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
#take marker genes in a 1 kb radius and save
#need to add distance to GWAS file and save

UNEAK_tags <- read.csv(file = "output_UNEAK_compare_1.csv", stringsAsFactors = FALSE, header = TRUE)
remove<-which(UNEAK_tags$Marker.name == "")

UNEAK_tags <- UNEAK_tags[-c(remove),]
UNEAK_tags <- UNEAK_tags$Marker.name


#distance_polyRAD <- distanceToNearest()
write.csv(polyRAD_2,file =  "polyRAD_distances_full.csv")
#names(full_genes) <- c("V1", "V2")
full_genes$V1 <- as.numeric(full_genes$V1)
in_gene <- length(which(full_genes$V1 == 0))
#polyRAD_2$V2 <- polyRAD_2$V2 -1
orange_5 <- length(which(full_genes$V1 < 5000))- in_gene
under_30 <- length(which(full_genes$V1 < 30000)) - sum(orange_5, in_gene)
over_30 <- length(which(full_genes$V2 > 30000))


#change directory to polyRAD directory
get_wrk_dir <- getwd()
GWAS_files <- list.files(path = get_wrk_dir  , pattern = "\\.GWAS") 
names(distance_UNEAK)[2] <- "Distance"



for (i in 1:length(GWAS_files)){
  GWAS_input <- GWAS_files[i]
  select_GWAS_file <- read.csv(file=GWAS_input, stringsAsFactors = FALSE, header = TRUE)
   match_UNEAK <- match(select_GWAS_file$SNP ,distance_UNEAK$Name)
  
  final_data_frame <-as.data.frame(cbind(select_GWAS_file, distance_TASSEL$Distance[c(match_UNEAK)]))
  names(final_data_frame)[10] <- "Distance"
  #remove .csv from file name
  GWAS_input <- sub("\\.csv", "", GWAS_input)
  write.csv(final_data_frame, file = paste(GWAS_input,i, ".csv" ,sep = ""))
}  

#correlation plots from distance analysis

library(ggplot2)
library(grid)
library(gridExtra)

get_wrk_dir <- getwd()
GWAS_files <- list.files(path = get_wrk_dir  , pattern = "\\.GWAS\\.Results[0-9].*") 
#grob3 = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(Distance, P.value), 4) ), x = 0.63, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))


#names(distance_UNEAK)[2] <- "Distance"

for (i in 1:length(GWAS_files)){
 
  GWAS_input <- GWAS_files[i]
  title_name <- gsub("\\.", " ", sub("\\.\\.", "", sub("year.*$","",sub("GAPIT", "", sub("\\." ," ",GWAS_input)))))
  df_correlation <- read.csv(file=GWAS_input, stringsAsFactors = FALSE, header = TRUE)
  attach(UNEAK_tag_sequences3)

  title_name <- paste("UNEAK", "(","R squared: ",as.character(round(cor(Distance, P.value), 4)),")" ,sep = "")                    
  #grob3 = grobTree(textGrob(paste("R squared :", corr_value)), x = 0.75, y = 0.85, hjust = 0, gp = gpar(col = "red", fontsize = 12, fontface = "bold"))
  grob3 = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(Distance, P.value), 4) ), x = 0.02, y = 0.92, hjust = 0, gp = gpar(col = "red", fontsize = 6, fontface = "bold")))
  title_name <- "polyRAD"
  polyRAD_plot <- ggplot(polyRAD_tag_sequences2, aes(y=-log10(P.value), x=Distance)) + geom_point() + 
    ggtitle(title_name) +
    geom_smooth(method=lm, se=FALSE) +
  scale_x_continuous(name = "Distance") +
  scale_y_continuous(name = expression(-log[10](P[value]))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  #annotation_custom(grob3) +
  scale_x_reverse()
  
  length(which(polyRAD_tag_sequences2$Distance == 0 ))
  length(UNEAK_tag_sequences3$Distance == 0)
  length(TASSEL_tag_sequences$Distance == 0)
  
  
  TASSEL_plot <- ggplot(TASSEL_tag_sequences, aes( x=Distance,y=-log10(P.value))) + geom_point() + 
    ggtitle(title_name) + 
    geom_smooth(method=lm, se=FALSE) +
   scale_x_continuous(name = "Distance") + 
    scale_y_continuous(name = expression(-log[10](P[value]))) + 
    #annotation_custom(grob3) + 
    theme(plot.title = element_text(hjust = 0.5), axis.line = element_line(color="black"), axis.line.x = element_line(color="black")) +
    scale_x_reverse()
  
  
  
  
   corr_plot  <-ggplot(UNEAK_tag_sequences3, aes(x=Distance, y=-log10(P.value))) + geom_point() + 
   ggtitle(title_name) + geom_smooth(method=lm, se=FALSE) +
   scale_x_continuous(name = "Distance (bp)") +
   scale_y_continuous(name = expression(-log[10](P[value]))) +
   theme(plot.title = element_text(hjust = 0.5)) +
   annotation_custom(grob3)  +
    scale_x_reverse() 
   #limits =c(2e+05,0)
   +
   annotation_custom(grob3) 
   #scale_x_sqrt()
   
  
   png(sub("\\.csv","",paste(GWAS_input,".png", sep ="")))
   print(corr_plot)
   dev.off()
  
  
  grid.arrange(polyRAD_plot, TASSEL_plot, UNEAK_plot,  nrow = 1, ncol = 3)
  
} 

polyRAD_GWAS_file <- read.csv("polyRAD.GAPIT.MLMM.TCmN.HU.NEF.UI.CSU.KNU.ZJU.year.3.GWAS.Results31.csv")
#Pearson correlation polyRAD -0.0399

polyRAD_tag_sequences2$Distance <- NA
match_names1 <- match(polyRAD_tag_sequences2$SNP, polyRAD_GWAS_file$SNP)
polyRAD_tag_sequences2$Distance <- polyRAD_GWAS_file$Distance[c(match_names1)]
#match_names2 <- match(UNEAK_tag_sequences3$P.value,UNEAK_tag_sequences$P.value)  
#UNEAK_tag_sequences3$Distance <- distance_col$V2
  
#print out spearman correlations 
corr_sp <-cor.test(TASSEL_tag_sequences$Distance, TASSEL_tag_sequences$P.value,  method = "spearman")
corr_sp <- cor.test(UNEAK_tag_sequences3$Distance, UNEAK_tag_sequences3$P.value,  method = "spearman")
corr_sp <- cor.test(polyRAD_tag_sequences2$Distance, polyRAD_tag_sequences2$P.value,  method = "spearman")


