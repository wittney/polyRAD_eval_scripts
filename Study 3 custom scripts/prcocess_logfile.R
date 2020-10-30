#process logfile from Hind/He pipeline
extractLog <- function(file, finalTemp = FALSE, tabu = FALSE){
  if(finalTemp && tabu){
    stop("Cannot have both finalTemp and tabu")
  }
  logout <- readLines(file)
  initHindHeLines <- grep("^Initial Hind/He:", logout)
  finalHindHeLines <- initHindHeLines + 2
  # filter out where isoloci were fixed
  keep <- grepl("^Final Hind/He:", logout[finalHindHeLines])
  initHindHeLines <- initHindHeLines[keep]
  finalHindHeLines <- finalHindHeLines[keep]
  
  markerLines <- initHindHeLines - 1
  initNMLines <- initHindHeLines + 1
  finalNMLines <- initHindHeLines + 3
  if(finalTemp){
    finalTempLines <- initHindHeLines + 4
  }
  if(tabu){
    repLines <- initHindHeLines + 4
  }
  
  # extract numeric values
  init_NM <- as.numeric(sub("Initial average NM: ", "", logout[initNMLines]))
  final_NM <- as.numeric(sub("Final average NM: ", "", logout[finalNMLines]))
  
  initHindHeSplit <- strsplit(logout[initHindHeLines], " ")
  finalHindHeSplit <- strsplit(logout[finalHindHeLines], " ")
  init1_HindHe <- as.numeric(sapply(initHindHeSplit,
                                    function(x){
                                      out <- x[3]
                                      out[out == "None"] <- NA
                                      return(out)}))
  init2_HindHe <- as.numeric(sapply(initHindHeSplit,
                                    function(x){
                                      out <- x[4]
                                      out[out == "None"] <- NA
                                      return(out)}))
  final1_HindHe <- as.numeric(sapply(finalHindHeSplit,
                                     function(x){
                                       out <- x[3]
                                       out[out == "None"] <- NA
                                       return(out)}))
  final2_HindHe <- as.numeric(sapply(finalHindHeSplit,
                                     function(x){
                                       out <- x[4]
                                       out[out == "None"] <- NA
                                       return(out)}))
  
  # get maximum initial and final Hind/He
  initMax_HindHe <- apply(cbind(init1_HindHe, init2_HindHe), 1, max, na.rm = TRUE)
  finalMax_HindHe <- apply(cbind(final1_HindHe, final2_HindHe), 1, max, na.rm = TRUE)
  
  # construct data frame
  outdf <- data.frame(row.names = logout[markerLines],
                      Init_NM = init_NM,
                      Final_NM = final_NM,
                      Init1_HindHe = init1_HindHe,
                      Init2_HindHe = init2_HindHe,
                      InitMax_HindHe = initMax_HindHe,
                      Final1_HindHe = final1_HindHe,
                      Final2_HindHe = final2_HindHe,
                      FinalMax_HindHe = finalMax_HindHe)
  if(finalTemp){
    outdf$Final_Temp <- as.numeric(sub("Final temperature: ", "", logout[finalTempLines]))
  }
  if(tabu){
    outdf$Reps <- as.integer(sub("Rep where best solution found: ", "", logout[repLines]))
  }
  
  return(outdf)
}

# extract values and plot ####
df_PstI <- extractLog("logfile_Msi_02192020")

df_NsiI <- extractLog("logfile_NsiI_02192020")

#combine df_PstI and df_NsiI
df_all <- rbind(df_PstI, df_NsiI)

zero_df <- which(df_all$FinalMax_HindHe == 0)
#855 Hind/HE was not able to be made


df_all <- df_all[-c(zero_df),]
not_sorted <- which(df_all$FinalMax_HindHe == "-Inf") 
NA_sorted <- which(df_all$FinalMax_HindHe != "-Inf")


df_all$not.sorted <- NA
df_all$not.sorted <- df_all$InitMax_HindHe
df_all$not.sorted[NA_sorted] <- 1

df_all$sorted <- NA
df_all$sorted <- df_all$InitMax_HindHe 
df_all$sorted[not_sorted] <- 1
#1467 not sorted out of 22,895
#total number of markers 22,895
#21,428 sorted
#855 removed from analysis



#vector the size of df_all
#HindHe_same <- vector(length = 22895L)
#HindHe_same <- rep(NA,22895 )
#for(i in 1:length(df_all$Init_NM)){
#input <- df_all$InitMax_HindHe [i] == df_all$FinalMax_HindHe[i]
#if(input == FALSE){
#HindHe_same[i] <- i
#  }
#}
#remove NA
#NA_remove <- which(is.na(HindHe_same))
#sort_tags <- which(is.na(HindHe_same))




ggplot(df_all, aes(x = InitMax_HindHe, y = FinalMax_HindHe, col = Final_NM - Init_NM)) +
  geom_point() +
  geom_hline(yintercept = 1/2) +
  geom_vline(xintercept = 1/2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_viridis() +
  coord_cartesian(xlim = c(0, 1.25), ylim = c(0, 1.25))
#fill = "#ffb79d"
#fill ="#9dffd0"

p <-ggplot(df_all) 
p <- p + geom_density(aes(x = InitMax_HindHe), colour="#005ab3", fill ="#ff7f50",alpha = 0.6, show.legend = TRUE) + 
  geom_density(aes(x = FinalMax_HindHe),  colour ="#005ab3", fill = "#00cdcd", alpha = 0.6, show.legend = TRUE ) +
  geom_vline(xintercept = 0.30, size = 1, colour = "black",linetype = "dashed") + labs(x="Hind/He", y = "Density")  
p<- p + coord_cartesian(xlim = c(0, 1.25), ylim = c(0, 7))

plot(df_all$Init_NM)
legend(x=4,y=15,c("No sorting","Sorting with Hind/He"),cex=.8,col=c("#00cdcd","#ff7f50"),pch=c(22,22), pt.cex =2, bg="lightgrey")

no_hind <- which(df_all$Reps == 0 )
df_all$No.Hind.He <- NA
df_all$No.Hind.He <- df_all$InitMax_HindHe[no_hind]
hind_he_used <- which(df_all$Reps != 0 )

df_all$No.Sorting <- NA
df_all$No.Sorting <- df_all$InitMax_HindHe
df_all$No.Sorting[NA_remove] <- 0

df_all$sorted <- NA
df_all$sorted <- df_all$InitMax_HindHe
df_all$sorted[sort_tags] <- 0

q <-ggplot(df_all) 
q <- q + geom_density(aes(x = sorted), colour="#005ab3", fill ="#ff7f50",alpha = 0.6, show.legend = TRUE) + 
  geom_density(aes(x = not.sorted),  colour ="#005ab3", fill = "#00cdcd", alpha = 0.6, show.legend = TRUE ) +
   labs(x="Hind/He", y = "Density")  
q <- q + coord_cartesian(xlim = c(0.30, 1.00), ylim = c(0, 5))

ggplot(df_all[df_all$InitMax_HindHe > 0.5,]) +
  geom_density(aes(x = FinalMax_HindHe), col = "blue") +
  geom_density(aes(x = InitMax_HindHe), col = "red") +
  geom_vline(xintercept = 0.5)

ggplot(df_all, aes(x = InitMax_HindHe, y = FinalMax_HindHe, col = Final_NM - Init_NM)) +
  geom_point() +
  geom_hline(yintercept = 1/2) +
  geom_vline(xintercept = 1/2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_viridis() +
  coord_cartesian(xlim = c(0, 1.25), ylim = c(0, 1.25))
