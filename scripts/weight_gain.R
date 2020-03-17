
# new weight measurements

setwd("/Users/12705859/Desktop/metapigs_base")
basedir = "/Users/12705859/Desktop/metapigs_base/phylosift/input_files/"

weights_final <- read_csv("~/Desktop/metapigs_base/weights_final.csv", 
                          col_types = cols(Pig = col_character(), 
                                           Room = col_character()))

weights_final <- pivot_wider(weights_final, id_cols=Pig,values_from=euth_wt,names_from =euth_day)
# get rid of rows after row 60 as they miss the value (no weights available for March 10th)
weights_final <- weights_final[1:60,]

weights <- read_csv("~/Desktop/metapigs_base/weights.csv", 
                    col_types = cols(Pig = col_character(), 
                                     Room = col_character()))


all_weights <- full_join(weights,weights_final)

colnames(all_weights)[colnames(all_weights) == 'Pig'] <- 'isolation_source'

all_weights <- all_weights[,3:12]

all_weights <- all_weights %>%
  pivot_longer(-isolation_source,
               names_to="date", 
               values_to="value",
               values_drop_na = TRUE)

# load metadata 
mdat <- read_excel(paste0(basedir,"Metagenome.environmental_20190308_2.xlsx"),
                   col_types = c("text", "numeric", "numeric", "text", "text",
                                 "text", "date", "text","text", "text", "numeric",
                                 "numeric", "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "text", "text","text", "text", "text", "text",
                                 "text","text", "text", "text", "text", "text","text", "text"),
                   skip = 12)

mdat <- mdat %>%
  dplyr::select(isolation_source, Cohort) %>%
  distinct()

# load details (breed, line, bday, mothers)
details <- read_excel(paste0(basedir, "pigTrial_GrowthWtsGE.hlsx.xlsx"),
                      "Piglet details")

# format details
colnames(details)[colnames(details) == 'STIG'] <- 'pig'
colnames(details)[colnames(details) == 'Nursing Dam'] <- 'nurse'
colnames(details)[colnames(details) == 'STIGDAM'] <- 'stig'
colnames(details)[colnames(details) == '...8'] <- 'breed'
details$pig <- gsub("G","",details$pig)
details$pig <- gsub("T","",details$pig)
colnames(details)[colnames(details) == 'pig'] <- 'isolation_source'

details$BIRTH_DAY <- as.character(details$BIRTH_DAY)
details$LINE <- as.character(details$LINE)
details <- details %>%
  dplyr::select(isolation_source,BIRTH_DAY,LINE,breed,stig,nurse)

df <- inner_join(all_weights,mdat)
df <- inner_join(df,details)
df


a <- aov(value ~ breed + date, df)
tidy(a)

a <- aov(value ~ BIRTH_DAY + date, df)
tidy(a)

a <- aov(value ~ BIRTH_DAY*breed + date, df)
tidy(a)

a <- aov(value ~ Cohort + date, df)
tidy(a)

a <- aov(value ~ Cohort*date, df)
tidy(a)


all_timepoints_cohorts_plot <- ggboxplot(df, x="date", y="value", fill = "Cohort", 
                                         ylab="weight (kg)")
all_timepoints_cohorts_plot

all_timepoints_breed_plot <- ggboxplot(df, x="date", y="value", fill = "breed", 
                                       ylab="weight (kg)")
all_timepoints_breed_plot

all_timepoints_bday_plot <- ggboxplot(df, x="date", y="value", fill = "BIRTH_DAY", 
                                      ylab="weight (kg)")
all_timepoints_bday_plot

all_timepoints_line_plot <- ggboxplot(df, x="date", y="value", fill = "LINE", 
                                      ylab="weight (kg)")
all_timepoints_line_plot

######################################################################################################

# DELTAS 


# filtering out piglets that had dysentery
df1 <- df %>%
  filter(!isolation_source=="29665"|isolation_source=="29865"|isolation_source=="29702")

pigs_1 <- df1 %>%
  filter(date == "31-Jan") %>%
  dplyr::select(isolation_source,Cohort,value,breed,BIRTH_DAY)
NROW(pigs_1)

###########################

pigs_2 <- df1 %>%
  filter(date == "7-Feb") %>%
  dplyr::select(isolation_source,Cohort,value,breed,BIRTH_DAY)
NROW(pigs_2)

###########################


pigs_3 <- df1 %>%
  filter(date == "14-Feb") %>%
  dplyr::select(isolation_source,Cohort,value,breed,BIRTH_DAY)
NROW(pigs_3)

###########################


pigs_4 <- df1 %>%
  filter(date == "21-Feb") %>%
  dplyr::select(isolation_source,Cohort,value,breed,BIRTH_DAY)
NROW(pigs_4)

###########################

pigs_5 <- df1 %>%
  filter(date == "28-Feb") %>%
  dplyr::select(isolation_source,Cohort,value,breed,BIRTH_DAY)
NROW(pigs_5)


##############################################################################

# settings for plots: 

#font size for pvalues 
your_font_size <- 2 

My_Theme = theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(), 
  axis.title.y = element_text(size = 8),
  axis.text.y = element_text(size = 8)) 

##############################################################################


df <- merge(pigs_1,pigs_2, by=c("isolation_source"))
df$diff = ((df$value.y-df$value.x)/df$value.y)*100

# reorder
df$Cohort.x <- factor(df$Cohort.x, 
                      levels=c("Control", 
                               "D-scour", 
                               "ColiGuard",
                               "Neomycin",
                               "Neomycin+D-scour",
                               "Neomycin+ColiGuard"))

res1 <- aov(diff ~ Cohort.x, data=df)
res <- TukeyHSD(res1)
aov.out <- as.data.frame(res$Cohort.x)
aov.out <- tibble::rownames_to_column(aov.out, "comparison")
aov.out$time_delta="Ja31_Fe7"
aov.out$test <- "anova"
aov.out$padj_method <- "TukeyHSD"
cohort_stats <- aov.out

cw_summary <- df %>% 
  group_by(Cohort.x) %>% 
  tally()

plot_A_B <- ggboxplot(df, x="Cohort.x", y="diff", fill = "Cohort.x", 
                      ylab="weight gain (%)", legend = "none")+
  My_Theme+
  geom_text(data = cw_summary,
            aes(Cohort.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 

plot_A_B


############

df <- merge(pigs_2,pigs_3, by=c("isolation_source"))
df$diff = ((df$value.y-df$value.x)/df$value.y)*100

# reorder
df$Cohort.x <- factor(df$Cohort.x, 
                      levels=c("Control", 
                               "D-scour", 
                               "ColiGuard",
                               "Neomycin",
                               "Neomycin+D-scour",
                               "Neomycin+ColiGuard"))

res1 <- aov(diff ~ Cohort.x, data=df)
res <- TukeyHSD(res1)
aov.out <- as.data.frame(res$Cohort.x)
aov.out <- tibble::rownames_to_column(aov.out, "comparison")
aov.out$time_delta="Fe7_Fe14"
aov.out$test <- "anova"
aov.out$padj_method <- "TukeyHSD"
cohort_stats <- rbind(cohort_stats,aov.out)

cw_summary <- df %>% 
  group_by(Cohort.x) %>% 
  tally()

plot_B_C <- ggboxplot(df, x="Cohort.x", y="diff", fill = "Cohort.x", 
                      ylab="weight gain (%)", legend = "none")+
  My_Theme+
  geom_text(data = cw_summary,
            aes(Cohort.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 

plot_B_C


############

df <- merge(pigs_3,pigs_4, by=c("isolation_source"))
df$diff = ((df$value.y-df$value.x)/df$value.y)*100

# reorder
df$Cohort.x <- factor(df$Cohort.x, 
                      levels=c("Control", 
                               "D-scour", 
                               "ColiGuard",
                               "Neomycin",
                               "Neomycin+D-scour",
                               "Neomycin+ColiGuard"))

res1 <- aov(diff ~ Cohort.x, data=df)
res <- TukeyHSD(res1)
aov.out <- as.data.frame(res$Cohort.x)
aov.out <- tibble::rownames_to_column(aov.out, "comparison")
aov.out$time_delta="Fe14_Fe21"
aov.out$test <- "anova"
aov.out$padj_method <- "TukeyHSD"
cohort_stats <- rbind(cohort_stats,aov.out)


cw_summary <- df %>% 
  group_by(Cohort.x) %>% 
  tally()

plot_C_D <- ggboxplot(df, x="Cohort.x", y="diff", fill = "Cohort.x", 
                      ylab="weight gain (%)", legend = "none")+
  My_Theme+
  geom_text(data = cw_summary,
            aes(Cohort.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 

plot_C_D


############

df <- merge(pigs_4,pigs_5, by=c("isolation_source"))
df$diff = ((df$value.y-df$value.x)/df$value.y)*100

# reorder
df$Cohort.x <- factor(df$Cohort.x, 
                      levels=c("Control", 
                               "D-scour", 
                               "ColiGuard",
                               "Neomycin",
                               "Neomycin+D-scour",
                               "Neomycin+ColiGuard"))

res1 <- aov(diff ~ Cohort.x, data=df)
res <- TukeyHSD(res1)
aov.out <- as.data.frame(res$Cohort.x)
aov.out <- tibble::rownames_to_column(aov.out, "comparison")
aov.out$time_delta="Fe21_Fe28"
aov.out$test <- "anova"
aov.out$padj_method <- "TukeyHSD"
cohort_stats <- rbind(cohort_stats,aov.out)

cw_summary <- df %>% 
  group_by(Cohort.x) %>% 
  tally()

plot_D_E <- ggboxplot(df, x="Cohort.x", y="diff", fill = "Cohort.x", 
                      ylab="weight gain (%)", legend = "none")+
  My_Theme+
  geom_text(data = cw_summary,
            aes(Cohort.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 

plot_D_E

############

df <- merge(pigs_2,pigs_4, by=c("isolation_source"))
df$diff = ((df$value.y-df$value.x)/df$value.y)*100

# reorder
df$Cohort.x <- factor(df$Cohort.x, 
                      levels=c("Control", 
                               "D-scour", 
                               "ColiGuard",
                               "Neomycin",
                               "Neomycin+D-scour",
                               "Neomycin+ColiGuard"))

res1 <- aov(diff ~ Cohort.x, data=df)
res <- TukeyHSD(res1)
aov.out <- as.data.frame(res$Cohort.x)
aov.out <- tibble::rownames_to_column(aov.out, "comparison")
aov.out$time_delta="Fe7_Fe21"
aov.out$test <- "anova"
aov.out$padj_method <- "TukeyHSD"
cohort_stats <- rbind(cohort_stats,aov.out)

cw_summary <- df %>% 
  group_by(Cohort.x) %>% 
  tally()

plot_B_D <- ggboxplot(df, x="Cohort.x", y="diff", fill = "Cohort.x", 
                      ylab="weight gain (%)", legend = "none")+
  My_Theme+
  geom_text(data = cw_summary,
            aes(Cohort.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 

plot_B_D

############

df <- merge(pigs_3,pigs_5, by=c("isolation_source"))
df$diff = ((df$value.y-df$value.x)/df$value.y)*100

# reorder
df$Cohort.x <- factor(df$Cohort.x, 
                      levels=c("Control", 
                               "D-scour", 
                               "ColiGuard",
                               "Neomycin",
                               "Neomycin+D-scour",
                               "Neomycin+ColiGuard"))

res1 <- aov(diff ~ Cohort.x, data=df)
res <- TukeyHSD(res1)
aov.out <- as.data.frame(res$Cohort.x)
aov.out <- tibble::rownames_to_column(aov.out, "comparison")
aov.out$time_delta="Fe14_Fe28"
aov.out$test <- "anova"
aov.out$padj_method <- "TukeyHSD"
cohort_stats <- rbind(cohort_stats,aov.out)

cw_summary <- df %>% 
  group_by(Cohort.x) %>% 
  tally()

plot_C_E <- ggboxplot(df, x="Cohort.x", y="diff", fill = "Cohort.x", 
                      ylab="weight gain (%)", legend = "none")+
  My_Theme+
  geom_text(data = cw_summary,
            aes(Cohort.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 

plot_C_E

############

df <- merge(pigs_1,pigs_5, by=c("isolation_source"))
df$diff = ((df$value.y-df$value.x)/df$value.y)*100

# reorder
df$Cohort.x <- factor(df$Cohort.x, 
                      levels=c("Control", 
                               "D-scour", 
                               "ColiGuard",
                               "Neomycin",
                               "Neomycin+D-scour",
                               "Neomycin+ColiGuard"))

res1 <- aov(diff ~ Cohort.x, data=df)
res <- TukeyHSD(res1)
aov.out <- as.data.frame(res$Cohort.x)
aov.out <- tibble::rownames_to_column(aov.out, "comparison")
aov.out$time_delta="Ja31_Fe28"
aov.out$test <- "anova"
aov.out$padj_method <- "TukeyHSD"
cohort_stats <- rbind(cohort_stats,aov.out)

cw_summary <- df %>% 
  group_by(Cohort.x) %>% 
  tally()

plot_A_E <- ggboxplot(df, x="Cohort.x", y="diff", fill = "Cohort.x", 
                      ylab="weight gain (%)", legend = "none")+
  My_Theme+
  geom_text(data = cw_summary,
            aes(Cohort.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 

plot_A_E

# this is for extracting the legend 
for_legend_only <- ggboxplot(df, x = "Cohort.x", y = "diff", fill = "Cohort.x", 
                             legend = "right")+
  scale_color_manual(labels = c("Control", 
                                "D-scour",
                                "ColiGuard",
                                "Neo",
                                "Neo+D",
                                "Neo+C"), 
                     values = c("#F8766D", 
                                "#B79F00",
                                "#00BA38",
                                "#00BFC4",
                                "#619CFF",
                                "#F564E3")) +
  guides(fill=guide_legend("Cohort")) +
  My_Theme
leg <- get_legend(for_legend_only)

empty_space = plot_grid(NULL, NULL, NULL, NULL, ncol=4)
top_row = plot_grid(plot_A_B,plot_B_C,plot_C_D,plot_D_E, ncol=4, 
                    rel_widths=c(0.25,0.25,0.25,0.25),
                    labels=c("A-B","B-C","C-D","D-E"),
                    label_size = 10)
bottom_row = plot_grid(plot_B_D,plot_C_E,plot_A_E,leg, ncol=4, 
                       rel_widths=c(0.25,0.25,0.25,0.25),
                       labels=c("B-D","C-E","A-E",""),
                       label_size = 10)
all_plots <- plot_grid(empty_space,
                       top_row,
                       bottom_row,
                       nrow=3)

pdf("weight_deltas_by_cohort.pdf")
ggdraw() +
  draw_image(timeline_deltas_weight, x = 0, y = 0.16) +
  draw_plot(all_plots)
dev.off()

# DELTA p-values - Cohort

# convert rownames to first column
weight_delta_cohorts <- setDT(cohort_stats, keep.rownames = TRUE)[]
# add data to workbook 
addWorksheet(wb, "weight_delta_cohorts")
writeData(wb, sheet = "weight_delta_cohorts", weight_delta_cohorts, rowNames = FALSE)

################################################################################################

# does the breed have anything to do with it??

##############################################################################

# breed

##############################################################################


df <- merge(pigs_1,pigs_2, by=c("isolation_source"))
df$diff = ((df$value.y-df$value.x)/df$value.y)*100

# reorder
df$breed.x <- factor(df$breed.x, 
                      levels=c("Landrace x Cross bred (LW x D)", 
                               "Duroc x Landrace", 
                               "Duroc x Large white",
                               "Large white x Duroc"))

res1 <- aov(diff ~ breed.x, data=df)
res <- TukeyHSD(res1)
aov.out <- as.data.frame(res$breed.x)
aov.out <- tibble::rownames_to_column(aov.out, "comparison")
aov.out$time_delta="Ja31_Fe7"
aov.out$test <- "anova"
aov.out$padj_method <- "TukeyHSD"
breed_stats <- aov.out

cw_summary <- df %>% 
  group_by(breed.x) %>% 
  tally()

plot_A_B <- ggboxplot(df, x="breed.x", y="diff", fill = "breed.x", 
                      ylab="weight gain (%)", palette = "jco", legend = "none")+
  My_Theme+
  geom_text(data = cw_summary,
            aes(breed.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 

plot_A_B

############

df <- merge(pigs_2,pigs_3, by=c("isolation_source"))
df$diff = ((df$value.y-df$value.x)/df$value.y)*100

# reorder
df$breed.x <- factor(df$breed.x, 
                     levels=c("Landrace x Cross bred (LW x D)", 
                              "Duroc x Landrace", 
                              "Duroc x Large white",
                              "Large white x Duroc"))

res1 <- aov(diff ~ breed.x, data=df)
res <- TukeyHSD(res1)
aov.out <- as.data.frame(res$breed.x)
aov.out <- tibble::rownames_to_column(aov.out, "comparison")
aov.out$time_delta="Fe7_Fe14"
aov.out$test <- "anova"
aov.out$padj_method <- "TukeyHSD"
breed_stats <- aov.out
breed_stats <- rbind(breed_stats,aov.out)

cw_summary <- df %>% 
  group_by(breed.x) %>% 
  tally()

plot_B_C <- ggboxplot(df, x="breed.x", y="diff", fill = "breed.x", 
                      ylab="weight gain (%)", palette = "jco", legend = "none")+
  My_Theme+
  geom_text(data = cw_summary,
            aes(breed.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 

plot_B_C

############

df <- merge(pigs_3,pigs_4, by=c("isolation_source"))
df$diff = ((df$value.y-df$value.x)/df$value.y)*100

# reorder
df$breed.x <- factor(df$breed.x, 
                     levels=c("Landrace x Cross bred (LW x D)", 
                              "Duroc x Landrace", 
                              "Duroc x Large white",
                              "Large white x Duroc"))


res1 <- aov(diff ~ breed.x, data=df)
res <- TukeyHSD(res1)
aov.out <- as.data.frame(res$breed.x)
aov.out <- tibble::rownames_to_column(aov.out, "comparison")
aov.out$time_delta="Fe14_Fe21"
aov.out$test <- "anova"
aov.out$padj_method <- "TukeyHSD"
breed_stats <- aov.out
breed_stats <- rbind(breed_stats,aov.out)

cw_summary <- df %>% 
  group_by(breed.x) %>% 
  tally()

plot_C_D <- ggboxplot(df, x="breed.x", y="diff", fill = "breed.x", 
                      ylab="weight gain (%)", palette = "jco", legend = "none")+
  My_Theme+
  geom_text(data = cw_summary,
            aes(breed.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 

plot_C_D


############

df <- merge(pigs_4,pigs_5, by=c("isolation_source"))
df$diff = ((df$value.y-df$value.x)/df$value.y)*100

# reorder
df$breed.x <- factor(df$breed.x, 
                     levels=c("Landrace x Cross bred (LW x D)", 
                              "Duroc x Landrace", 
                              "Duroc x Large white",
                              "Large white x Duroc"))

res1 <- aov(diff ~ breed.x, data=df)
res <- TukeyHSD(res1)
aov.out <- as.data.frame(res$breed.x)
aov.out <- tibble::rownames_to_column(aov.out, "comparison")
aov.out$time_delta="Fe21_Fe28"
aov.out$test <- "anova"
aov.out$padj_method <- "TukeyHSD"
breed_stats <- aov.out
breed_stats <- rbind(breed_stats,aov.out)

cw_summary <- df %>% 
  group_by(breed.x) %>% 
  tally()

plot_D_E <- ggboxplot(df, x="breed.x", y="diff", fill = "breed.x", 
                      ylab="weight gain (%)", palette = "jco", legend = "none")+
  My_Theme+
  geom_text(data = cw_summary,
            aes(breed.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 

plot_D_E

############

df <- merge(pigs_2,pigs_4, by=c("isolation_source"))
df$diff = ((df$value.y-df$value.x)/df$value.y)*100

# reorder
df$breed.x <- factor(df$breed.x, 
                     levels=c("Landrace x Cross bred (LW x D)", 
                              "Duroc x Landrace", 
                              "Duroc x Large white",
                              "Large white x Duroc"))

res1 <- aov(diff ~ breed.x, data=df)
res <- TukeyHSD(res1)
aov.out <- as.data.frame(res$breed.x)
aov.out <- tibble::rownames_to_column(aov.out, "comparison")
aov.out$time_delta="Fe7_Fe21"
aov.out$test <- "anova"
aov.out$padj_method <- "TukeyHSD"
breed_stats <- aov.out
breed_stats <- rbind(breed_stats,aov.out)

cw_summary <- df %>% 
  group_by(breed.x) %>% 
  tally()

plot_B_D <- ggboxplot(df, x="breed.x", y="diff", fill = "breed.x", 
                      ylab="weight gain (%)", palette = "jco", legend = "none")+
  My_Theme+
  geom_text(data = cw_summary,
            aes(breed.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 

plot_B_D

############

df <- merge(pigs_3,pigs_5, by=c("isolation_source"))
df$diff = ((df$value.y-df$value.x)/df$value.y)*100

# reorder
df$breed.x <- factor(df$breed.x, 
                     levels=c("Landrace x Cross bred (LW x D)", 
                              "Duroc x Landrace", 
                              "Duroc x Large white",
                              "Large white x Duroc"))

res1 <- aov(diff ~ breed.x, data=df)
res <- TukeyHSD(res1)
aov.out <- as.data.frame(res$breed.x)
aov.out <- tibble::rownames_to_column(aov.out, "comparison")
aov.out$time_delta="Fe14_Fe28"
aov.out$test <- "anova"
aov.out$padj_method <- "TukeyHSD"
breed_stats <- aov.out
breed_stats <- rbind(breed_stats,aov.out)


cw_summary <- df %>% 
  group_by(breed.x) %>% 
  tally()

plot_C_E <- ggboxplot(df, x="breed.x", y="diff", fill = "breed.x", 
                      ylab="weight gain (%)", palette = "jco", legend = "none")+
  My_Theme+
  geom_text(data = cw_summary,
            aes(breed.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 

plot_C_E


# this is for extracting the legend 
for_legend_only <- ggboxplot(df, x = "breed.x", y = "diff", fill = "breed.x", 
                             legend = "right")+
  scale_fill_manual(labels = c("Landrace x Cross bred (LW x D)", 
                                "Duroc x Landrace",
                                "Large white x Duroc",
                                "Duroc x Large white"
                                ), 
                     values = c("#0073C2FF",   # first 4 colors of jco palette
                                "#EFC000FF",
                                "#868686FF",
                                "#CD534CFF"
                                )) +
  guides(fill=guide_legend("breed")) +
  My_Theme

leg <- get_legend(for_legend_only)


empty_space = plot_grid(NULL, NULL, NULL, NULL, ncol=4)
top_row = plot_grid(plot_A_B,plot_B_C,plot_C_D,plot_D_E, ncol=4, 
                    rel_widths=c(0.25,0.25,0.25,0.25),
                    labels=c("A-B","B-C","C-D","D-E"),
                    label_size = 10)
bottom_row = plot_grid(plot_B_D,plot_C_E, leg, ncol=3, 
                       rel_widths=c(0.25,0.25,0.50),
                       labels=c("B-D","C-E",""),
                       label_size = 10)
all_plots <- plot_grid(empty_space,
                       top_row,
                       bottom_row,
                       nrow=3)

# pdf("weight_deltas_by_breed.pdf")
# ggdraw() +
#   draw_image(timeline_deltas, x = 0, y = 0.16) +
#   draw_plot(all_plots)
# dev.off()


# DELTA p-values - Breed

# convert rownames to first column
weight_delta_breed <- setDT(breed_stats, keep.rownames = TRUE)[]
# add data to workbook 
addWorksheet(wb, "weight_delta_breed")
writeData(wb, sheet = "weight_delta_breed", weight_delta_breed, rowNames = FALSE)

##############################################################################

# bday 

##############################################################################


df <- merge(pigs_1,pigs_2, by=c("isolation_source"))
df$diff = ((df$value.y-df$value.x)/df$value.y)*100

# reorder
df$BIRTH_DAY.x <- factor(df$BIRTH_DAY.x, 
                      levels=c("2017-01-06",
                               "2017-01-07",
                               "2017-01-08",
                               "2017-01-09",
                               "2017-01-10",
                               "2017-01-11"))

res1 <- aov(diff ~ BIRTH_DAY.x, data=df)
res <- TukeyHSD(res1)
aov.out <- as.data.frame(res$BIRTH_DAY.x)
aov.out <- tibble::rownames_to_column(aov.out, "comparison")
aov.out$time_delta="Ja31_Fe7"
aov.out$test <- "anova"
aov.out$padj_method <- "TukeyHSD"
bday_stats <- aov.out

cw_summary <- df %>% 
  group_by(BIRTH_DAY.x) %>% 
  tally()

plot_A_B <- ggboxplot(df, x="BIRTH_DAY.x", y="diff", fill = "BIRTH_DAY.x", 
                      ylab="weight gain (%)", palette = "jco", legend = "none")+
  My_Theme+
  geom_text(data = cw_summary,
            aes(BIRTH_DAY.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 

plot_A_B

############

df <- merge(pigs_2,pigs_3, by=c("isolation_source"))
df$diff = ((df$value.y-df$value.x)/df$value.y)*100

# reorder
df$BIRTH_DAY.x <- factor(df$BIRTH_DAY.x, 
                         levels=c("2017-01-06",
                                  "2017-01-07",
                                  "2017-01-08",
                                  "2017-01-09",
                                  "2017-01-10",
                                  "2017-01-11"))

res1 <- aov(diff ~ BIRTH_DAY.x, data=df)
res <- TukeyHSD(res1)
aov.out <- as.data.frame(res$BIRTH_DAY.x)
aov.out <- tibble::rownames_to_column(aov.out, "comparison")
aov.out$time_delta="Fe7_Fe14"
aov.out$test <- "anova"
aov.out$padj_method <- "TukeyHSD"
bday_stats <- rbind(bday_stats,aov.out)

cw_summary <- df %>% 
  group_by(BIRTH_DAY.x) %>% 
  tally()

plot_B_C <- ggboxplot(df, x="BIRTH_DAY.x", y="diff", fill = "BIRTH_DAY.x", 
                      ylab="weight gain (%)", palette = "jco", legend = "none")+
  My_Theme+
  geom_text(data = cw_summary,
            aes(BIRTH_DAY.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 

plot_B_C

############

df <- merge(pigs_3,pigs_4, by=c("isolation_source"))
df$diff = ((df$value.y-df$value.x)/df$value.y)*100

# reorder
df$BIRTH_DAY.x <- factor(df$BIRTH_DAY.x, 
                         levels=c("2017-01-06",
                                  "2017-01-07",
                                  "2017-01-08",
                                  "2017-01-09",
                                  "2017-01-10",
                                  "2017-01-11"))

res1 <- aov(diff ~ BIRTH_DAY.x, data=df)
res <- TukeyHSD(res1)
aov.out <- as.data.frame(res$BIRTH_DAY.x)
aov.out <- tibble::rownames_to_column(aov.out, "comparison")
aov.out$time_delta="Fe14_Fe21"
aov.out$test <- "anova"
aov.out$padj_method <- "TukeyHSD"
bday_stats <- rbind(bday_stats,aov.out)

cw_summary <- df %>% 
  group_by(BIRTH_DAY.x) %>% 
  tally()

plot_C_D <- ggboxplot(df, x="BIRTH_DAY.x", y="diff", fill = "BIRTH_DAY.x", 
                      ylab="weight gain (%)", palette = "jco", legend = "none")+
  My_Theme+
  geom_text(data = cw_summary,
            aes(BIRTH_DAY.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 

plot_C_D

############

df <- merge(pigs_4,pigs_5, by=c("isolation_source"))
df$diff = ((df$value.y-df$value.x)/df$value.y)*100

# reorder
df$BIRTH_DAY.x <- factor(df$BIRTH_DAY.x, 
                         levels=c("2017-01-06",
                                  "2017-01-07",
                                  "2017-01-08",
                                  "2017-01-09",
                                  "2017-01-10",
                                  "2017-01-11"))

res1 <- aov(diff ~ BIRTH_DAY.x, data=df)
res <- TukeyHSD(res1)
aov.out <- as.data.frame(res$BIRTH_DAY.x)
aov.out <- tibble::rownames_to_column(aov.out, "comparison")
aov.out$time_delta="Fe21_Fe28"
aov.out$test <- "anova"
aov.out$padj_method <- "TukeyHSD"
bday_stats <- rbind(bday_stats,aov.out)

cw_summary <- df %>% 
  group_by(BIRTH_DAY.x) %>% 
  tally()

plot_D_E <- ggboxplot(df, x="BIRTH_DAY.x", y="diff", fill = "BIRTH_DAY.x", 
                      ylab="weight gain (%)", palette = "jco", legend = "none")+
  My_Theme+
  geom_text(data = cw_summary,
            aes(BIRTH_DAY.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 

plot_D_E

############

df <- merge(pigs_2,pigs_4, by=c("isolation_source"))
df$diff = ((df$value.y-df$value.x)/df$value.y)*100

# reorder
df$BIRTH_DAY.x <- factor(df$BIRTH_DAY.x, 
                         levels=c("2017-01-06",
                                  "2017-01-07",
                                  "2017-01-08",
                                  "2017-01-09",
                                  "2017-01-10",
                                  "2017-01-11"))

res1 <- aov(diff ~ BIRTH_DAY.x, data=df)
res <- TukeyHSD(res1)
aov.out <- as.data.frame(res$BIRTH_DAY.x)
aov.out <- tibble::rownames_to_column(aov.out, "comparison")
aov.out$time_delta="Fe7_Fe21"
aov.out$test <- "anova"
aov.out$padj_method <- "TukeyHSD"
bday_stats <- rbind(bday_stats,aov.out)

cw_summary <- df %>% 
  group_by(BIRTH_DAY.x) %>% 
  tally()

plot_B_D <- ggboxplot(df, x="BIRTH_DAY.x", y="diff", fill = "BIRTH_DAY.x", 
                      ylab="weight gain (%)", palette = "jco", legend = "none")+
  My_Theme+
  geom_text(data = cw_summary,
            aes(BIRTH_DAY.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 

plot_B_D

############

df <- merge(pigs_3,pigs_5, by=c("isolation_source"))
df$diff = ((df$value.y-df$value.x)/df$value.y)*100

# reorder
df$BIRTH_DAY.x <- factor(df$BIRTH_DAY.x, 
                         levels=c("2017-01-06",
                                  "2017-01-07",
                                  "2017-01-08",
                                  "2017-01-09",
                                  "2017-01-10",
                                  "2017-01-11"))

res1 <- aov(diff ~ BIRTH_DAY.x, data=df)
res <- TukeyHSD(res1)
aov.out <- as.data.frame(res$BIRTH_DAY.x)
aov.out <- tibble::rownames_to_column(aov.out, "comparison")
aov.out$time_delta="Fe14_Fe28"
aov.out$test <- "anova"
aov.out$padj_method <- "TukeyHSD"
bday_stats <- rbind(bday_stats,aov.out)

cw_summary <- df %>% 
  group_by(BIRTH_DAY.x) %>% 
  tally()

plot_C_E <- ggboxplot(df, x="BIRTH_DAY.x", y="diff", fill = "BIRTH_DAY.x", 
                      ylab="weight gain (%)", palette = "jco", legend = "none")+
  My_Theme+
  geom_text(data = cw_summary,
            aes(BIRTH_DAY.x, Inf, label = n), vjust="inward", size = your_font_size)+
  stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 

plot_C_E

# this is for extracting the legend 
for_legend_only <- ggboxplot(df, x = "BIRTH_DAY.x", y = "diff", fill = "BIRTH_DAY.x", 
                             legend = "right")+
  scale_fill_manual(labels = c("2017-01-06", 
                               "2017-01-07",
                               "2017-01-08",
                               "2017-01-09",
                               "2017-01-10",
                               "2017-01-11"
                               ), 
                    values = c("#0073C2FF",   # first 6 colors of jco palette
                               "#EFC000FF",
                               "#868686FF",
                               "#CD534CFF",
                               "#7AA6DCFF",
                               "#003C67FF"
                    ))+
  guides(fill=guide_legend("birth day")) +
  My_Theme

leg <- get_legend(for_legend_only)

empty_space = plot_grid(NULL, NULL, NULL, NULL, ncol=4)
top_row = plot_grid(plot_A_B,plot_B_C,plot_C_D,plot_D_E, ncol=4, 
                    rel_widths=c(0.25,0.25,0.25,0.25),
                    labels=c("A-B","B-C","C-D","D-E"),
                    label_size = 10)
bottom_row = plot_grid(plot_B_D,plot_C_E, leg, ncol=3, 
                       rel_widths=c(0.25,0.25,0.50),
                       labels=c("B-D","C-E",""),
                       label_size = 10)
all_plots <- plot_grid(empty_space,
                       top_row,
                       bottom_row,
                       nrow=3)

# pdf("weight_deltas_by_bday.pdf")
# ggdraw() +
#   draw_image(timeline_deltas, x = 0, y = 0.16) +
#   draw_plot(all_plots)
# dev.off()


# DELTA p-values - bday

# convert rownames to first column
weight_delta_bday <- setDT(bday_stats, keep.rownames = TRUE)[]
# add data to workbook 
addWorksheet(wb, "weight_delta_bday")
writeData(wb, sheet = "weight_delta_bday", weight_delta_bday, rowNames = FALSE)


# save stats in existing workbook
saveWorkbook(wb, "/Users/12705859/Desktop/metapigs_base/phylosift/out/stats.xlsx", overwrite=TRUE)




