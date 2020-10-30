#retrieve names of markers which contain chromosome number and position from the alignment file
#outputted in file
#need to remove lines that doesnt have anything in the second column 
data_align <- read.csv("final_11032019_1_align.csv", stringsAsFactors = FALSE, header = TRUE)
#14931 tag sequences total from wehat population simulated with RADinito
remove_lines <- which(data_align$Alignment.2 == "")
new_data_align <- data_align[-c(remove_lines),]
#3,437 aligned only once
#23% percent of the tags aligned only once
#77% of tags aligned to more than one position in the wheat data
#read in the sorted tags file
sorted_tags <- read.csv("final_11032019_1_sorted.csv", stringsAsFactors = FALSE, header = TRUE)
#remove tag names that doesnt match
#this will remove tag sequences that were removed from the final output
data_sorted <- sorted_tags[!is.na(match(sorted_tags$Tag.sequence, new_data_align$Tag.sequence)),]
#25% of the tags that were aligned more than once were able to be adequately sorted by the pipeline
#file to match with data
#text file provided by radinito which contains the correct chromsome and positions of each tag sequence included
#in the study
radinitio_file <- read.delim2(file = "modified_file_7.txt", header = FALSE)
#new column to have the chrom number in radinitio file
#adjustments made for direct comparison can be made with polyRAD outuput
radinitio_file$V5 <- as.character(sub(">t.*[np]", "", radinitio_file$V1))
radinitio_file$V6 <- sprintf("%09d", radinitio_file$V3)
radinitio_file$v7 <- paste("ta_iwgsc", radinitio_file$V5, sep = "_")
radinitio_file$v8 <- paste(radinitio_file$v7, "v1", sep = "_")
radinitio_file$v9 <- paste(radinitio_file$v7, "v2", sep = "_")
radinitio_file$v10 <- paste(radinitio_file$v7, "v3", sep = "_")

radinitio_file$v8 <- paste(radinitio_file$v8, radinitio_file$V2, sep = "_")
radinitio_file$v9 <- paste(radinitio_file$v9, radinitio_file$V2, sep = "_")
radinitio_file$v10 <- paste(radinitio_file$v10, radinitio_file$V2, sep = "_")

radinitio_file$v8 <- paste(radinitio_file$v8, radinitio_file$V6, sep = "-")
radinitio_file$v9 <- paste(radinitio_file$v9, radinitio_file$V6, sep = "-")
radinitio_file$v10 <- paste(radinitio_file$v10, radinitio_file$V6, sep = "-")
#for sorted file remove -top and -bot from marker
data_sorted$new_Marker <- sub("\\-[topbo].*$", "", data_sorted$Marker)
data_sorted$Marker = data_sorted$new_Marker

match_1 <- match(data_sorted$Marker,radinitio_file$v8)
match_2 <- match(data_sorted$Marker,radinitio_file$v9)
match_3 <- match(data_sorted$Marker,radinitio_file$v10)
match_all <- as.integer(match_all)
average_cols <- colMeans(match_all)
not_aligned <- which(match_all == NA && match_2 == NA && match_3 == NA) 
match_all <- t(rbind(match_1, match_2, match_3))
data_sorted$Chr <- as.character(sub("ta_iwgsc_", "",data_sorted$Marker))
data_sorted$Chr <- as.character(sub("_.*$", "", data_sorted$Chr))
data_sorted$P1 <- as.character(sub("\\-.*$", "",data_sorted$Marker))
data_sorted$P1 <- as.integer(sub("ta_iwgsc_.*_.", "",data_sorted$P1))

data_sorted$P2 <- as.character(sub("-[topb].*$", "",data_sorted$Marker))
data_sorted$P2 <- as.integer(sub("^.*-", "",data_sorted$P2 ))

#radinito file 2
radinitio_file2 <-read.delim2(file = "reference_rad_loci.txt", header = TRUE, stringsAsFactors = FALSE )
keep_lines <- which(radinitio_file2$status == "kept")
radinitio_file2 <- radinitio_file2[c(keep_lines),]
remove_lines <- which(radinitio_file2$cut_pos == NA)
radinitio_file2$new_cut_pos <- sprintf("%09d", radinitio_file2$cut_pos)
#new colum with chrom id and new cutsite
radinitio_file2$marker <- paste(radinitio_file2$X.chrom_id, radinitio_file2$new_cut_pos, sep='-')
#if tag-dir is postive put -top, if negative -bot on marker
radinitio_file2$new_marker <- radinitio_file2$new_cut_pos
i=1
for(i in 1:nrow(radinitio_file2)){
  
  if (radinitio_file2$tag_dir[i] == "pos"){
    radinitio_file2$new_marker[i] <- paste(radinitio_file2$marker[i], "top", sep = "-")
    
  }
  
  else {
    radinitio_file2$new_marker[i] <- paste(radinitio_file2$marker[i], "bot", sep = "-")
    
  }
}
no_match <- which(data_sorted$Marker == radinitio_file2$new_marker)
matched <- data_sorted[!is.na(match(data_sorted$Marker, radinitio_file2$new_marker)),]
#file used gives wrong position.. use correct position in the AM
average_rows <- rowMeans(match_all, na.rm = TRUE)
not_aligned <-which(average_rows == "NaN")
#133 out of 2883 did not align properly