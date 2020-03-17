######################################################################################################

library(readxl)
library(data.table)

# rationale: 

##### previous steps: 


# 1 # guppy fat is run 

# 3 # .xml conversion to .txt:

# run forester.jar from command line to convert the .xml file to phyloXML - R readable format (.txt) : 
# this way: 
# java -cp /Users/12705859/Downloads/forester_1050.jar 
# org.forester.application.phyloxml_converter -f=dummy file.xml file.txt
# in a loop: 
# for fpath in /Users/12705859/Desktop/metapigs_base/phylosift/guppy/guppy_output/*.gz.xml; 
# do java -cp /Users/12705859/Downloads/forester_1050.jar 
# org.forester.application.phyloxml_converter -f=dummy "$fpath" 
# "/Users/12705859/Desktop/metapigs_base/phylosift/guppy/guppy_output/$(basename "$fpath").txt"; 
# done

##### HERE : 

# 1 # .xml files are read in and parsed

# 2 # output can be used in guppy_plots.R

######################################################################################################


basedir = "/Users/12705859/Desktop/metapigs_base/phylosift/input_files/"
my.basedir <- "~/Desktop/metapigs_base/phylosift/guppy/guppy_output"


my.files = list.files(my.basedir,pattern=".gz.xml.txt")

# construct an empty dataframe to build on 
complete.df <- data.frame(
  taxa = character(),
  branch_length = character(),
  branch_width = character(),
  red = character(),
  green = character(),
  blue = character(),
  file = character()
)

for (textfile in my.files) {
  
  # read in file 
  my.df <- read_csv(file.path(my.basedir,textfile), col_names = FALSE)
  
  # extract file name 
  myfilename <- basename(textfile)
  
  # even though in this case we only have one element, we grep "<phylogeny rooted"
  # as this will discard the unnecessary first rows of the file 
  startOFnewPC <- grep("<phylogeny rooted", my.df$X1)
  z <- split(my.df, cumsum(1:nrow(my.df) %in% startOFnewPC))
  my.df <- as.data.frame(z$`1`)
  
  # start the grepping of useful tree info(we only care about the hits containing "width")
  my_width_ref <- grep("width", my.df$X1)
  my_width_ref <- lapply(my_width_ref, print)
  
  # starting from "width" we retrieve all the useful info that belong to each node
  df <- data.frame(matrix(unlist(my_width_ref), nrow=length(my_width_ref), byrow=T))
  colnames(df) <- "C"
  
  mysel <- df %>%
    mutate(B=C-1) %>%
    mutate(A=B-1) %>%
    mutate(D=C+1) %>%
    mutate(E=D+1) %>%
    mutate(FF=E+1) %>%
    mutate(G=FF+1) %>%
    mutate(placeholder="placeholder") %>%
    select(A,B,C,D,E,FF,G) %>%
    pivot_longer(cols=A:G)
  
  mysel <- as.data.frame(mysel)
  row.names(mysel) <- mysel$value

  almostthere <- my.df[match(rownames(mysel), rownames(my.df), nomatch=0),]
  almostthere <- as.data.frame(almostthere)

  almostthere <- cSplit(almostthere, "almostthere","</")

  almostthere <- almostthere %>% 
  filter(almostthere_2 != '') # drop empty rows

  #almostthere$index <- rownames(almostthere)
  almostthere$almostthere_2 <- NULL

  # transpose this single columned df every 6 elements 
  tmp <- data.frame(
    X=almostthere$almostthere_1,
    ind=rep(1:6, nrow(almostthere)/6)
  )
  
  # unstack 
  myprecious <- unstack(tmp, X~ind)
  
  # in this case we won't need the colors as there's only one color possible combo possible 
  myprecious <- myprecious[,1:3]
  
  # rename cols
  colnames(myprecious) <- c('name', 'branch_length', 'branch_width')
  
  # remove all the symbols derived from the xml format
  myprecious <- as.data.frame(lapply(myprecious, function(y) gsub("<[^>]+>", "", y)))
  
  # save file name as an extra column
  myprecious$file <- myfilename
  
  # bind all dfs from all files 
  complete.df <- rbind(
    complete.df, 
    myprecious
  )
  
}
  

complete.df <- cSplit(complete.df, "file","_")
complete.df$DNA_plate <- paste0(complete.df$file_1,"_",complete.df$file_2)
complete.df$DNA_well <- complete.df$file_3

complete.df <- complete.df %>%
  select(name,branch_length,branch_width,DNA_plate,DNA_well)



View(test)

test <- complete.df
length(unique(test$name))


# keep everything between __ and }
test$name <- test$name %>%  
  gsub('.*__(.*)\\}.*', '', .)  

test$name

test <- cSplit(test, "name","__")
test <- test %>%
  select(name_01,name_02,name_03,name_04,name_05,name_06,name_07,name_08,name_09,name_10)






View(test)
test$name <- test$name %<>%
  gsub('\\{|\\}', '', .) %>% # removes curly brackets
  gsub('\\[|\\]', '', .) %>% # removes square brackets
  gsub('^_', '', .) %>% # removes the first _
  gsub('[0-9]+', '', .) %>% # removes digits
  gsub('_$', '', .) %>%  # removes the last _
  gsub('.*__', '', .)  # removes everyting up to __ (keeping only the most specific)

length(unique(test$name))




simplified.df <- complete.df

NROW(simplified.df)
simplified.df$name <- simplified.df$name %<>%
  gsub('\\{|\\}', '', .) %>% # removes curly brackets
  gsub('\\[|\\]', '', .) %>% # removes square brackets
  gsub('^_', '', .) %>% # removes the first _
  gsub('[0-9]+', '', .) %>% # removes digits
  gsub('_$', '', .) %>%  # removes the last _
  gsub('.*__', '', .)  # removes everyting up to __ (keeping only the most specific)

fwrite(x = simplified.df, file = "/Users/12705859/Desktop/metapigs_base/phylosift/guppy/guppyfat_simplified")
fwrite(x = complete.df, file = "/Users/12705859/Desktop/metapigs_base/phylosift/guppy/guppyfat_complete")


unique(simplified.df$name)
which(is.na(simplified.df$name))
