#run total_RAD

load("Nsi_RAD.Rdata")
load("seperate_RAD_data_objects.Rdata")
load("allSam.Rdata")
#RAD data objects names NsiI_RAD and PstI_RAD
#subset polyRAD
keep_taxa <- which(rownames(PstI_RAD$alleleDepth) %in% rownames(NsiI_RAD$alleleDepth) == TRUE)
PstI_RAD <- SubsetByTaxon(PstI_RAD, taxa = allSam)
NsiI_RAD <- SubsetByTaxon(NsiI_RAD, taxa = allSam)
# remove any loci duplicated across the two sets

# remove any loci duplicated across the two sets
nLoci(PstI_RAD)    # 119008
nLoci(NsiI_RAD)    # 186936
nAlleles(PstI_RAD) # 490854
nAlleles(NsiI_RAD) # 955678
NsiI_keeploci <- which(!GetLoci(NsiI_RAD) %in% GetLoci(PstI_RAD))
cat(nLoci(NsiI_RAD) - length(NsiI_keeploci), 
    file = "180522Num_duplicate_loci.txt") #992 duplicate
NsiI_RAD <- SubsetByLocus(NsiI_RAD, NsiI_keeploci)



# combine allele depth into one matrix
PstI_depth <- PstI_RAD$alleleDepth
NsiI_depth <- NsiI_RAD$alleleDepth
total_depth <- matrix(0L, nrow = length(GD.all.taxa$taxa), 
                      ncol = ncol(PstI_depth) + ncol(NsiI_depth),
                      dimnames = list(GD.all.taxa$taxa, 
                                      c(colnames(PstI_depth), 
                                        colnames(NsiI_depth))))
total_depth[,colnames(PstI_depth)] <- PstI_depth[GD.all.taxa$taxa,]
total_depth[rownames(NsiI_depth),colnames(NsiI_depth)] <- NsiI_depth

# combine other slots
total_alleles2loc <- c(PstI_RAD$alleles2loc, 
                       NsiI_RAD$alleles2loc + nLoci(PstI_RAD))
total_locTable <- rbind(PstI_RAD$locTable, NsiI_RAD$locTable)
total_alleleNucleotides <- c(PstI_RAD$alleleNucleotides, 
                             NsiI_RAD$alleleNucleotides)

# build new RADdata object and save
total_RAD <- RADdata(total_depth, total_alleles2loc, total_locTable,
                     list(2L), 0.001, total_alleleNucleotides)
#remove the blank sample and any sample that we need to
#could equal 563

load("final_Taxa.Rdata")
keep_taxa2 <- which(new_allSam %in% rownames(total_RAD$alleleDepth) == TRUE)
new_allSam <- new_allSam[c(keep_taxa2)]
total_RAD <- SubsetByTaxon(total_RAD, taxa = new_allSam)
#last total 562 of 

save(total_RAD, file = "02192020_RADdata_NsiIPstI_560.RData")

#remove Ch
locTable <- as.data.frame(as.character(sub("Chr", "", total_RAD$locTable$Chr)))
names(locTable)[1] <- "Chr" 
total_RAD$locTable$Chr <- locTable$Chr

# Make groups representing pairs of chromosomes, and one group for all 
# non-assembled scaffolds.
splitlist <- list(c("^01$", "^02$"),
                  c("^03$", "^04$"),
                  c("^05$", "^06$"),
                  c("^07$", "^08$"),
                  c("^09$", "^10$"),
                  c("^11$", "^12$"),
                  c("^13$", "^14$", "^15$"),
                  c("^16$", "^17$"),
                  c("^18$", "^19$"))
# split by chromosome and save seperate objects
SplitByChromosome(total_RAD, chromlist = splitlist, 
                  chromlist.use.regex = TRUE, fileprefix = "180524splitRAD")

# files with RADdata objects
splitfiles <- grep("^180524splitRAD", list.files("."), value = TRUE)

save(total_RAD, file = "total_RAD_560_022020.Rdata")



# list to hold markers formatted for GAPIT/FarmCPU
GAPITlist <- list()
length(GAPITlist) <- length(splitfiles)

# loop through RADdata objects
for(i in 1:length(splitfiles)){
  load(splitfiles[i])
  splitRADdata <- IteratePopStructLD(splitRADdata)
  GAPITlist[[i]] <- ExportGAPIT(splitRADdata)
}
save(GAPITlist, file = "180524GAPITlist.RData")

# put together into one dataset for FarmCPU
GM.all <- rbind(GAPITlist[[1]]$GM, GAPITlist[[2]]$GM, GAPITlist[[3]]$GM,
                GAPITlist[[4]]$GM, GAPITlist[[5]]$GM, GAPITlist[[6]]$GM, 
                GAPITlist[[7]]$GM, GAPITlist[[8]]$GM,
                GAPITlist[[9]]$GM, GAPITlist[[10]]$GM)
GD.all <- cbind(GAPITlist[[1]]$GD, GAPITlist[[2]]$GD[,-1],
                GAPITlist[[3]]$GD[,-1], GAPITlist[[4]]$GD[,-1],
                GAPITlist[[5]]$GD[,-1], GAPITlist[[6]]$GD[,-1],
                GAPITlist[[7]]$GD[,-1], GAPITlist[[8]]$GD[,-1], 
                GAPITlist[[9]]$GD[,-1], GAPITlist[[10]]$GD[,-1])
save(GD.all, GM.all, file = "02142020GM_GD_all_polyRAD.RData") # 1076888 markers
