#loop for CMplot 

library(CMplot)
library(RAINBOW)
man.plots <- function() {
  getwd_for_function <- getwd()
  dir_files <- list.files(path = getwd_for_function , pattern = "\\.GWAS") 
  #extract trait name from file
 #extract for UNEAK header
  
  for(i in 1:length(dir_files)){
    dir_files <- dir_files[order(dir_files)]
    
    pval <- read.csv(dir_files[i], stringsAsFactors = FALSE, header = TRUE)
    bool <- grepl("^polyRAD.*$", dir_files[i])
    bool1 <-grepl("^TASSEL.*$", dir_files[i])
    {
    if (bool == TRUE) {
      sig <- 1.0e-05
      header <- "polyRAD"
    }
    else if (bool1 == TRUE) {
      #significant for TASSEL  
      sig <- 1.0-06
      header <- "TASSEL"
    }
    else {
      #significant for UNEAK
      sig <- 1.0-06
      header <- "UNEAK"
    }
    }
    
  nCHR <- unique(pval$Chromosome)
  #chrcol <- c("#b35a00", "#00b3b3", "#b3005a" )
  chrcol <- c("#008080","#00cdcd")
  pval1 <- pval[,c(1,2,3,4)]
  
  #need to add header
  
  CMplot( pval1 , plot.type="m", col=chrcol, LOG10=TRUE, ylim=NULL, threshold = NULL,threshold.lty=c(1,2),
         threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
         chr.den.col=c("#34ffff","#00e6e6", "#009a9a") ,signal.col=c("#9a0000"),signal.cex=c(1,1),
         signal.pch=c(19,19),file="tiff",memo=header,dpi=300,file.output=TRUE,verbose=TRUE,
         width=14,height=4)
  
  
    
  }
 