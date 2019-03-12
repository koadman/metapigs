#input file: Tiz_vs_aminoglycoside_resultDB.m8 (mmseqs search output)


levels(factor(Tiz_vs_aminoglycoside_resultDB$X1))
b <- as.data.frame(unique(Tiz_vs_aminoglycoside_resultDB$X1))
names(b)[1] <- "X1"
View(b)
#add empty column
b[,"newName"] <- NA
#gives 235 entries

#select only rows with values > 0.95
library(tidyverse)
c <- Tiz_vs_aminoglycoside_resultDB %>% filter(X3 >= 0.80)
View(c)

#then merge unique(Tiz_vs_aminoglycoside_resultDB$X1) with it
#and see which values are missing in X1 
#c_b <- merge(c,b,by="X1")
b_c <- left_join(b, c, by = c("X1"))
View(b_c)

DF_b_c <- b_c[is.na(b_c$X2),]
View(DF_b_c)

#filtering to only get values >= 0.93 and 
#merging to unique genes dataframe, I get 24 "NA", 
#which means 24 were below 0.93 coverage 
#<=0.98 gives 145 genes
#<=0.95 gives 72 genes 
#<=0.93 gives 24 genes
#<=0.90 gives 21 genes 
#<=0.80 fives 19 genes


#get only one column
levels(factor(DF_b_c$X1))
part_match <- as.data.frame(unique(DF_b_c$X1))
names(part_match)[1] <- "X1"
View(part_match)

#obtain all data for those 24 genes by joining the 24-genes dataframe with the original dataframe
part_match_Tiz <- left_join(part_match, Tiz_vs_aminoglycoside_resultDB, by = c("X1"))
View(part_match_Tiz)

#conclusion: 
#24 genes of Tiziana'd database (235 genes) have been covered below 0.93
#18 genes of Tiziana'd database (235 genes) have been covered below 0.80

write.table(part_match_Tiz, "part_match_Tiz", sep="\t")



