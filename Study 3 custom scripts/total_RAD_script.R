library(polyRAD)
library(VariantAnnotation)



NsiI_file <- "170705Msi_NsiI_genotypes.vcf"
PstI_file <- "170608Msi_PstI_genotypes.vcf"

#allSam <- covariate_matrix1$Genotype


# The vector allSam was defined outside of this script, and contains the 
# names of all samples that I wanted to import.  Below I find sample names
# within the VCF files that match those samples.

NsiI_RAD <- VCF2RADdata(NsiI_file, yieldSize = 5e4,
                         expectedAlleles = 1e6, expectedLoci = 2e5)
save(NsiI_RAD, file = "Nsi_RAD.Rdata")
NsiI_RAD <SubsetByTaxon
PstI_RAD <- VCF2RADdata(PstI_file, yieldSize = 5e4,
                        expectedAlleles = 1e6, expectedLoci = 2e5)
save(PstI_RAD, file = "seperate_RAD_data_objects.Rdata")
#filter both by samples
names_RAD <- rownames(NsiI_RAD$alleleDepth)
allSam <- names_RAD[names_RAD %in% GD.all7$taxa]
save(allSam, file = "allSam.Rdata")

# remove any loci duplicated across the two sets
nLoci(PstI_RAD)    # 116757
nLoci(NsiI_RAD)    # 187434
nAlleles(PstI_RAD) # 478210
nAlleles(NsiI_RAD) # 952511
NsiI_keeploci <- which(!GetLoci(NsiI_RAD) %in% GetLoci(PstI_RAD))
cat(nLoci(NsiI_RAD) - length(NsiI_keeploci), 
    file = "180522Num_duplicate_loci.txt") #992 duplicate
NsiI_RAD <- SubsetByLocus(NsiI_RAD, NsiI_keeploci)

# combine allele depth into one matrix
PstI_depth <- PstI_RAD$alleleDepth
NsiI_depth <- NsiI_RAD$alleleDepth
total_depth <- matrix(0L, nrow = length(allSam), 
                      ncol = ncol(PstI_depth) + ncol(NsiI_depth),
                      dimnames = list(allSam, 
                                      c(colnames(PstI_depth), 
                                        colnames(NsiI_depth))))
total_depth[,colnames(PstI_depth)] <- PstI_depth[allSam,]
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
save(total_RAD, file = "02142020_RADdata_NsiIPstI.RData")

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
                  c("^18$", "^194"), "^SCAFFOLD")
# split by chromosome and save seperate objects
SplitByChromosome(total_RAD, chromlist = splitlist, 
                  chromlist.use.regex = TRUE, fileprefix = "180524splitRAD")

# files with RADdata objects
splitfiles <- grep("^180524splitRAD", list.files("."), value = TRUE)

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
