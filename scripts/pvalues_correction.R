








myfun_pvalue_df_prep1 <- function(df1_here,n) {
  # first, join the PCs to be in one columns, so you can easily apply the t.test function to one column only 
  a <- df1_here %>%
    select(PC1,PC2,PC3,PC4,PC5,Cohort,collection_date,guppied_date,sample_type) %>%
    pivot_longer(
      cols = 1:5,
      names_to = "component",
      values_to = "value",
      values_drop_na = FALSE
    )
  
  # look at distribution - filter out if less than n observations, per component, per collection date, per cohort 
  df1 <- setDT(a)[, .(Freq = .N), by = .(collection_date,Cohort,component)]
  df1 <- df1 %>% filter(!Freq<n)
  a <- merge(a,df1)
  a <- a %>% 
    select(Cohort,collection_date,guppied_date,sample_type,component,value)
  return(a)
}

myfun_pvalue_df_prep2 <- function(df2_here) {
  
  # create an empty df to build on
  df_pval <- data.frame(
    group_2 = character(),
    group_1 = character(),
    p_value = numeric(),
    which_PC = character(),
    which_colldate = character()
  )
  
  a <- df2_here
  
  # for each of the components (PC1 to PC5) ...
  listcompos <- unique(a$component)
  
  for (compo in listcompos) {
    df <- a %>% filter(component==compo)
    
    listcoldates <- unique(df$guppied_date)
    
    for (colldate in listcoldates) {
      
      z <- pairwise.t.test(df$value[df$guppied_date==colldate], df$Cohort[df$guppied_date==colldate], p.adjust = "none")$p.value
      z <- as.data.frame(z)
      z
      z$group_2 <- rownames(z)
      rownames(z) <- NULL
      z <- z %>%
        pivot_longer(
          cols = -group_2,
          names_to = "group_1",
          values_to = "p_value",
          values_drop_na = TRUE
        ) 
      
      z <- z %>% mutate(which_colldate = colldate,
                        which_PC = compo)
      df_pval <- rbind(df_pval,z) 
    }
  }
  return(df_pval)
}

myfun_pvalue_adjust <- function(df_pval,withgrouping) {
  
  # Final adjustments 
  
  df_pval <- df_pval %>%
    # join the two groups for which the p-value has been computed
    mutate(comparison=paste0(group_1,"_vs_",group_2)) %>%
    select(p_value,which_colldate,which_PC,comparison,group_1,group_2)
  
  
  # filtering to keep only meaningful comparisons 
  # to be kept: 
  meaningfulcomparisons <- c("Control_vs_ColiGuard", "ColiGuard_vs_Control",
                             "Control_vs_D-scour", "D-scour_vs_Control",
                             "Control_vs_Neomycin", "Neomycin_vs_Control",
                             "Neomycin_vs_Neomycin+D-scour", "Neomycin+D-scour_vs_Neomycin",
                             "Neomycin_vs_Neomycin+ColiGuard", "Neomycin+ColiGuard_vs_Neomycin")
  
  # eliminate useless comparisons
  df_pval <- df_pval[df_pval$comparison %in% meaningfulcomparisons,]
  
  if (withgrouping == "yes") {
    df_pval <- df_pval %>%
      group_by(which_colldate) %>% 
      mutate(pval.adj = p.adjust (p_value, method='fdr'))
    return(df_pval)
  } else {
    df_pval <- df_pval %>%
      mutate(pval.adj = p.adjust (p_value, method='fdr'))
    return(df_pval)
  }
  
}




# DF_piggies_time,groupA,groupB,groupC,groupD,groupE,groupF,groupG

all <- rbind(DF_piggies_time,
             groupA,
             groupB,
             groupC,
             groupD,
             groupE,
             groupF,
             groupG)

all$groupsplit <- paste0(all$sample_type,"_",all$guppied_date)

# splitting into multiple dataframes (by file name)
multi_DFs <- split( all , f = all$groupsplit )

# prep empty df to build on 
significant <- data.frame(
  group_1 = character(),
  group_2 = character(),
  p_value = numeric(),
  which_colldate = character(),
  which_PC = character(),
  comparison = character(),
  pval.adj = numeric(),
  groupsplit = character(),
  stringsAsFactors = FALSE
)

for (singl_DF in multi_DFs) {
  a1 <- myfun_pvalue_df_prep1(singl_DF,2)
  a2 <- myfun_pvalue_df_prep2(a1)
  a3 <- myfun_pvalue_adjust(a2,"no")
  a4 <- a3 %>% 
    filter(pval.adj<0.05) %>%
    mutate(groupsplit = singl_DF$groupsplit[1])
  a4 <- as.data.frame(a4)
  significant <- rbind(
    significant,
    a4)
}


significant <- cSplit(significant, "groupsplit","_")
colnames(significant)[colnames(significant)=="groupsplit_1"] <- "dataframe"
colnames(significant)[colnames(significant)=="groupsplit_2"] <- "guppied_date"

# keep only useful cols
significant <- significant %>%
  select(guppied_date,which_PC,group_1,group_2,dataframe,pval.adj)


# compare the number of significant observations obtained 
# with and without grouping by collection_date
yesss <- significant
NROW(yesss)
NROW(nooo)

# compare these values with the values from the loop 

#################################################
#################################################
#################################################
#################################################

minidf <- unique(data.frame(DF_piggies_time$guppied_date,
                     DF_piggies_time$collection_date))
colnames(minidf) <- c("guppied_date","collection_date")

df <- inner_join(significant,minidf)
df

df$dataframe <- gsub("piggies","DF_piggies",df$dataframe)


b <- eval(as.name(paste(df$dataframe[1])))  %>%
  filter(collection_date==as.character(df$collection_date[1])) %>%
  filter(Cohort==as.character(df$group_1[1])|Cohort==as.character(df$group_2[1])) %>%
  ggplot( aes(x=eval(as.name(paste(df$which_PC[1]))), fill=Cohort)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity')+
  labs(fill="") +
  ggtitle(as.character(df$guppied_date[1]))
b


l <- list()
for (A in rownames(df)) {
  A <- as.numeric(A)
  p <- eval(as.name(paste(df$dataframe[A])))  %>%
    filter(collection_date==as.character(df$collection_date[A])) %>%
    filter(Cohort==as.character(df$group_1[A])|Cohort==as.character(df$group_2[A])) %>%
    ggplot( aes(x=eval(as.name(paste(df$which_PC[A]))), fill=Cohort)) +
    geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity')+
    labs(fill="") +
    ggtitle(as.character(df$guppied_date[A]))+
    xlab(as.character(df$which_PC[A]))
  name <- paste("p",A,sep="_")
  tmp <- list(p)
  l[[name]] <- tmp
}

print(l$p_8)

getwd()
pdf("test_plots.pdf")
for (i in 1:18) {
  print(l[[i]])
}
dev.off()




#ks.test(b$PC2[b$Cohort=="Neomycin"],b$PC2[b$Cohort=="Neomycin+D-scour"])


