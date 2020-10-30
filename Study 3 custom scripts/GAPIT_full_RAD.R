library(polyRAD)
load("total_RAD_560_022020.Rdata")
#loading in total_RAD
# Make groups representing pairs of chromosomes, and one group for all 
# non-assembled scaffolds.
#splitlist <- list(c("^01$", "^02$"),
#                  c("^03$", "^04$"),
#                  c("^05$", "^06$"),
#                  c("^07$", "^08$"),
#                  c("^09$", "^10$"),
#                  c("^11$", "^12$"),
#                  c("^13$", "^14$", "^15$"),
#                  c("^16$", "^17$"),
#                  c("^18$", "^19$", "^SCAFFOLD"))
# split by chromosome and save seperate objects
#SplitByChromosome(total_RAD, chromlist = splitlist, 
#                  chromlist.use.regex = TRUE, fileprefix = "180524splitRAD")

# files with RADdata objects
#splitfiles <- grep("^180524splitRAD", list.files("."), value = TRUE)

# list to hold markers formatted for GAPIT/FarmCPU
#GAPITlist <- list()
#length(GAPITlist) <- length(splitfiles)

# loop through RADdata objects
#for(i in 1:length(splitfiles)){
#  load(splitfiles[i])
  total_RAD <- IteratePopStructLD(total_RAD)
  GAPIT_export <- ExportGAPIT(total_RAD)
#}
save(total_RAD, file = "02212020_totalRAD_Iterate.Rdata")
save(GAPIT_export, file = "GAPIT_export_02212020.Rdata")

# put together into one dataset for FarmCPU
#GM.all <- rbind(GAPITlist[[1]]$GM, GAPITlist[[2]]$GM, GAPITlist[[3]]$GM,
#                GAPITlist[[4]]$GM, GAPITlist[[5]]$GM, GAPITlist[[6]]$GM, 
#                GAPITlist[[7]]$GM, GAPITlist[[8]]$GM,
#                GAPITlist[[9]]$GM, GAPITlist[[10]]$GM)
#GD.all <- cbind(GAPITlist[[1]]$GD, GAPITlist[[2]]$GD[,-1],
#                GAPITlist[[3]]$GD[,-1], GAPITlist[[4]]$GD[,-1],
#                GAPITlist[[5]]$GD[,-1], GAPITlist[[6]]$GD[,-1],
#                GAPITlist[[7]]$GD[,-1], GAPITlist[[8]]$GD[,-1], 
#                GAPITlist[[9]]$GD[,-1], GAPITlist[[10]]$GD[,-1])
#save(GD.all, GM.all, file = "02142020GM_GD_all_Discrete.RData") # 1076888 markers

discrete.genotypes <- GetProbableGenotypes(total_RAD)

discrete.genotypes.final <- as.data.frame(discrete.genotypes$genotypes)

col_names <- colnames(discrete.genotypes.final)
correct_names <- which(col_names %in% GM.final$Name == TRUE)

GD.final <- discrete.genotypes.final[,correct_names ]

GD.final.1 <- cbind(rownames(GD.final), GD.final)
save(GD.final.1,GM.final, file = "correct_GD_GM_02222020.Rdata")

#must load GD and GM
#load("/media/HDD/wittney2/Miscanthus_sorting_pipeline/polyRAD_variant_calling_02192020/02202020GM_GD_all_Discrete.RData")
#load in covariate matrix
#load("/media/HDD/wittney2/Miscanthus_sorting_pipeline/polyRAD_variant_calling_02192020/02202020GM_GM_all_Discrete.RData")
load("phenotype_info_02202020.Rdata")
load("correct_GD_GM_02222020.Rdata")

#load("phenotype_discrete_genotypes_02172020.Rdata")
load("covariate_matrix_discrete_02172020.Rdata")
#load("phenotype_discrete_genotypes_2.Rdata")

#covariate_matrix2 <- covariate_matrix[covariate_matrix$Genotype %in% GD.all$taxa, ]
#save(covariate_matrix2, file = "covariate_matrix_discrete_02172020.Rdata")
#phenotype_mis1 <- read.csv(file = "gcbb12620-sup-0001-datas1-1.csv", header = TRUE, stringsAsFactors = FALSE)
##phenotype_mis2 <- phenotype_mis1[phenotype_mis1$Genotype %in% covariate_matrix2$Genotype, ]
#phenotype_mis3 <- phenotype_mis2[,c(1,76,79,81,145,148,150,168,171,173,205,208,210,243,246,248,266,269,271,289,292,294,307,310,312,330, 333,335, 99, 102, 104, 29, 32, 34, 53, 56, 58, 122, 125, 127)]
#GD.all.discrete2 <- GD.all.discrete[rownames(GD.all.discrete) %in% covariate_matrix2$Genotype,]


#GD.all.discrete = GD.all.discrete2
#save(GD.all.discrete, GM.all ,file = "02142020GM_GD_all_Discrete.RData")

#need to remove 2 observations from GD.all.discrete to match covariate matrix
library(GAPIT3)

myGAPIT_MLMM <- GAPIT(
  Y=phenotype_mis4,
  GD=GD.final.1 ,
  GM= GM.final,
  model="MLMM",
  CV=covariate_matrix2,
  file.output=T
)
#GD.all.discrete <- as.data.frame(GD.all.discrete)
#load("GD_discrete_dataframe_02202020.Rdata")
#GM.all.2 <- GM.all[1:1025348,]
#GM.all = GM.all.2

#tag_names <- colnames(GD.all.discrete)
#GM.all.1 <- GM.all[GM.all$Name %in% tag_names,]
#GD.all.discrete2 <- GD.all.discrete[,colnames(GD.all.discrete) %in% GM.all.1$Name]

#unique_tags <- unique(colnames(GD.all.discrete2))
#unique_names <- unique(GM.all$Name)
#GD.all.discrete2 <- as.data.frame(GD.all.discrete)
#GD.all.taxa <- as.data.frame(rownames(GD.all.discrete2))
#colnames(GD.all.taxa)[1] <- "taxa"
#GD.all.discrete3 <- cbind(GD.all.taxa,GD.all.discrete2)
#GD.all.discrete =GD.all.discrete3

#GM.all.2 <- GM.all[GM.all$Name %in% colnames(GD.all.discrete),]
#GD.all.discrete = GD.all.discrete2
#GM.all = GM.all.2
#x <- 26:40
#phenotype_mis4 <- phenotype_mis3[,c(1,x)]
#save(phenotype_mis4, file = "phenotype_discrete_genotypes_2.Rdata")
#write.csv(GD.all.taxa, file = "sample_txt_560_02192020.csv")

#total for TASSEL SNPs 1024980


