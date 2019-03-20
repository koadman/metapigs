#
# Using a GTDBtk summary .tsv file and the bin3C cluster report, produce a treemap figure
# of taxonomic assignment, where boxes are sized by coverage of each MAG.
#
# Note, this script does not attempt to parse your input arguments or give help.
#
# Usage: [gtdbtk summary file] [bin3C cluster_report.csv] [output pdf]
#

library(treemap)
library(d3treeR)
library(dplyr)
library(tidyr)

args = commandArgs(trailingOnly=TRUE)

# read gtdb bacterial tax
d = read.csv(args[1], sep='\t', header=TRUE)

# split the classification into separate columns
d = separate(d, classification, sep=";", into=c("domain", "phylum", "class", "order", "family", "genus", "species"))

# remove the leading part of the strings
f = function(x){gsub("^.__", "", x)}
d = mutate(d, domain=f(domain), phylum=f(phylum), class=f(class), order=f(order), family=f(family), genus=f(genus), species=f(species))

# join the classification with the cluster report
d = inner_join(d, read.csv(args[2]), by=c('user_genome' = 'name'))

# make a treemap, with boxes sized by cluster coverage
pdf(args[3])
treemap(d, index=colnames(d)[3:7], vSize="cov_expect", type="index", palette='Pastel1', title="", title.legend="", aspRatio=1)
dev.off()

