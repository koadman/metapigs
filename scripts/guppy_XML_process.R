######################################################################################################

# rationale: 

##### previous steps: 

# 1 # groups for guppy are made in guppy_group.R 

# 2 # guppy is run 

# 3 # .xml conversion to .txt:

# run forester.jar from command line to convert the .xml file to phyloXML - R readable format (.txt) : 
# this way: 
# java -cp /Users/12705859/Downloads/forester_1050.jar 
# org.forester.application.phyloxml_converter -f=dummy file.xml file.txt
# in a loop: 
# for fpath in /Users/12705859/Desktop/metapigs_base/phylosift/guppy/guppy_output/*.xml; 
# do java -cp /Users/12705859/Downloads/forester_1050.jar 
# org.forester.application.phyloxml_converter -f=dummy "$fpath" 
# "/Users/12705859/Desktop/metapigs_base/phylosift/guppy/guppy_output/$(basename "$fpath").txt"; 
# done

##### HERE : 

# 1 # .xml files are read in and parsed

# 2 # output can be used in guppy_plots.R


######################################################################################################

# 0 # set working directory & load libs

setwd("/Users/12705859/Desktop/metapigs_base/phylosift/guppy")

library(readr)
library(dplyr)
library(tidyr)
library(splitstackshape)
library(grid)
library(data.table)

###########################################################################################

my.basedir <- "~/Desktop/metapigs_base/phylosift/guppy/guppy_output"

# construct an empty dataframe to build on 
complete.df <- data.frame(
  var_explained = character(),
  branch_length = character(),
  branch_width = character(),
  blue = character(),
  green = character(),
  red = character(),
  taxa = character(),
  component = character(),
  PC_position = character(),
  file = character()
  )

my.files = list.files(my.basedir,pattern=".txt.xml.txt")
my.files <- my.files[-17] # removing problematic file (won't parse - it's Ja31 anyway)
my.files <- my.files[-26] # removing problematic file (won't parse - it's Ja31 anyway)

for (textfile in my.files) {
  
  # read in file 
  my.df <- read_csv(file.path(my.basedir,textfile), col_names = FALSE)
  
  # extract file name 
  myfilename <- basename(textfile)
  
  startOFnewPC <- grep("<phylogeny rooted", my.df$X1)
  
  z <- split(my.df, cumsum(1:nrow(my.df) %in% startOFnewPC))
  
  PC1 <- as.data.frame(z$`1`)
  PC1$component = "PC1"
  PC1$var_explained <- as.character(z$`1`[2,])
  PC2 <- as.data.frame(z$`2`)
  PC2$component = "PC2"
  PC2$var_explained <- as.character(z$`2`[2,])
  PC3 <- as.data.frame(z$`3`)
  PC3$component = "PC3"
  PC3$var_explained <- as.character(z$`3`[2,])
  PC4 <- as.data.frame(z$`4`)
  PC4$component = "PC4"
  PC4$var_explained <- as.character(z$`4`[2,])
  PC5 <- as.data.frame(z$`5`)
  PC5$component = "PC5"
  PC5$var_explained <- as.character(z$`5`[2,])
  
  my.df <- rbind(PC1,PC2,PC3,PC4,PC5)
  
  # start the grepping of useful tree info 
  mylist <- grep("blue", my.df$X1)
  myend <- lapply(mylist, print)
  
  df <- data.frame(matrix(unlist(myend), nrow=length(myend), byrow=T))
  colnames(df) <- "G"
  mysel <- df %>%
    mutate(FF=G-1) %>%
    mutate(E=FF-1) %>%
    mutate(D=E-1) %>%
    mutate(C=D-1) %>%
    mutate(B=C-1) %>%
    mutate(A=B-1) %>%
    mutate(placeholder="placeholder")%>%
    pivot_longer(cols=A:G)
  
  mysel <- as.data.frame(mysel)
  row.names(mysel) <- mysel$value
  
  almostthere <- my.df[match(rownames(mysel), rownames(my.df), nomatch=0),]
  almostthere <- cSplit(almostthere, "X1","</")
  almostthere$X1_2 <- NULL
  
  myprecious <- almostthere %>% 
    filter(X1_1 != '') %>%    # drop empty rows
    mutate(key = rep(c('taxa', 'branch_length', 'branch_width', 'color','red','green','blue'), n() / 7), 
           id = cumsum(key == 'taxa')) %>% 
    spread(key, X1_1) %>%
    dplyr::select(var_explained,branch_length,branch_width,blue,green,red,taxa,component)
  
  # remove all the symbols derived from the xml format
  myprecious <- as.data.frame(lapply(myprecious, function(y) gsub("<[^>]+>", "", y)))
  
  # round the variation explained by the PC, down to two digits 
  myprecious$var_explained <-round(as.numeric(as.character(myprecious[,1])) * 100,2)
  
  myprecious$PC_position <- paste0(myprecious$blue,"_",myprecious$green,"_",myprecious$red)
  # only two unique combos that make up for green or red (unique(myprecious$PC_position)
  
  # green standing for higher up in PC, red the opposite
  myprecious$PC_position <- gsub("165_194_102", "up",myprecious$PC_position)
  myprecious$PC_position <- gsub("98_141_252", "down",myprecious$PC_position)
  
  myprecious$file <- myfilename
  
  complete.df <- rbind(
    complete.df, 
    myprecious
  )
  
}


complete <- complete.df

unique(jplace_df$file)
unique(complete$file)

# clean file name & parse file name 
complete$file <- gsub('_sel.txt.xml.txt', '', complete$file)
complete$file <- gsub('piggies_group_A', 'groupA', complete$file)
complete$file <- gsub('piggies_group_B', 'groupB', complete$file)
complete$file <- gsub('piggies_CTRLNEO', 'groupC', complete$file)
complete$file <- gsub('piggies_NEONEOD', 'groupD', complete$file)
complete$file <- gsub('piggies_NEONEOC', 'groupE', complete$file)
complete$file <- gsub('piggies_CTRLDs', 'groupF', complete$file)
complete$file <- gsub('piggies_CTRLC', 'groupG', complete$file)
complete$file <- gsub('^piggies$', 'piggies_all', complete$file)
complete <- cSplit(complete, "file","_")

colnames(complete)[colnames(complete) == 'file_1'] <- 'sample_type'
colnames(complete)[colnames(complete) == 'file_2'] <- 'guppied_date'

NROW(complete)
unique(complete$sample_type)

###############

# simplify taxa


complete$taxa_simple <- complete$taxa %<>%
  gsub('\\{|\\}', '', .) %>% # removes curly brackets
  gsub('^_', '', .) %>% # removes the first _
  gsub('[0-9]+', '', .) %>% # removes digits
  gsub('_$', '', .) %>%  # removes the last _
  gsub('.*__', '', .)  # removes everything up to __ (keeping only the most specific)

sort(unique(complete$taxa_simple))

###############

View(complete)
colnames(complete)

# Store the minimum necessary info - no branch width included as we need this xml data 
# for the axes description only (as we don't need this info to this purpose, we can do "distinct")
simplified <- complete %>%
  select(sample_type, guppied_date,var_explained,
         component,PC_position,taxa_simple) %>%
  distinct()

# save both complete and simplified dataframes 
fwrite(x = complete, file = "guppy_xml_complete.df")
fwrite(x = simplified, file = "guppy_xml_simplified.df")



# functions to collect taxa associated with PCs:

#############################


# how to use function created below: PC_up(find_PC5(pos_controls))

# what grid.text is extpecting: e.g. : PC_up(find_PC5(df))    in this case (df = pos_controls) and derives from: 
# ```
# pos_controls <- simplified %>%
#   filter(sample_type=="pos") %>%
#   group_split(component) 
# ```
# and has as columns: 
# colnames(pos_controls[[1]])
# [1] "sample_type"   "guppied_date"  "var_explained" "component"     "PC_position"   "taxa_simple"  


######################################################################################################
######################################################################################################

# functions to find PC of interest


find_PC1 <- function(x) {
  PC <- x[[1]] %>%
    group_by(taxa_simple) %>%
    filter(n()==1) %>% # keep only taxa that appear once either up or down in PC, not both
    arrange(PC_position,taxa_simple)
  return(PC)
}
find_PC2 <- function(x) {
  PC <- x[[2]] %>%
    group_by(taxa_simple) %>%
    filter(n()==1) %>% # keep only taxa that appear once either up or down in PC, not both
    arrange(PC_position,taxa_simple)
  return(PC)
}
find_PC3 <- function(x) {
  PC <- x[[3]] %>%
    group_by(taxa_simple) %>%
    filter(n()==1) %>% # keep only taxa that appear once either up or down in PC, not both
    arrange(PC_position,taxa_simple)
  return(PC)
}
find_PC4 <- function(x) {
  PC <- x[[4]] %>%
    group_by(taxa_simple) %>%
    filter(n()==1) %>% # keep only taxa that appear once either up or down in PC, not both
    arrange(PC_position,taxa_simple)
  return(PC)
}
find_PC5 <- function(x) {
  PC <- x[[5]] %>%
    group_by(taxa_simple) %>%
    filter(n()==1) %>% # keep only taxa that appear once either up or down in PC, not both
    arrange(PC_position,taxa_simple)
  return(PC)
}

######################################################################################################
######################################################################################################

# functions to get taxa going up or down the PC

PC_down <- function(x) {
  down <- paste(as.list((x$taxa_simple[x$PC_position=="down"]),"\n"),collapse="\n")
  return(down)
}


PC_up <- function(x) {
  up <- paste(as.list((x$taxa_simple[x$PC_position=="up"]),"\n"),collapse="\n")
  return(up)
}

# interrogate the function by typing 
# PC_up(anydataframeyouwant)


######################################################################################################
######################################################################################################

# functions to get taxa going up or down the PC

get_var <- function(x) {
  var <- paste(as.list((x$var_explained)[1]))   # first item only 
  return(var)
}

# interrogate the function by typing 
# get_var(anydataframeyouwant)


######################################################################################################
######################################################################################################

######################################################################################################
######################################################################################################


