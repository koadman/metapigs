library(dplyr)
library(splitstackshape)

mash_output_all <- read_table2("~/Desktop/metapigs_dry/mash_output_all.csv", col_names = FALSE)
NROW(mash_output_all)


mash_out <- mash_output_all %>%
  dplyr::filter(!X4==1)
NROW(mash_out)
View(mash_out)


mash_out <- cSplit(mash_out,"X1","/") 
mash_out <- cSplit(mash_out,"X2","/")
mash_out <- cSplit(mash_out,"X5","/")

# remove .fa extension to match bins in checkm df 
mash_out$X2_6 <- gsub(".fa","", mash_out$X2_6)
mash_out$X1_7 <- gsub(".contigs.fasta","", mash_out$X1_7)

mash_out <- cSplit(mash_out,"X2_6","_")

mash_out <- mash_out %>%
  dplyr::mutate(match_percentage=(X5_1/X5_2)*100 ) %>%
  dplyr::select(X1_7,X2_6_1,X2_6_2,X3,X4,match_percentage)

colnames(mash_out) <- c("reference","isolation_source", "bin", "mash_dist","p_value","match_percentage")


mash_out <- mash_out %>%
  dplyr::filter(p_value==0)
head(mash_out)

View(mash_out)

