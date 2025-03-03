rm(list=ls())
library(tidyverse)
library(nnls)
library(glmnet)
library(ggplot2)
library(showtext)
library(gridExtra)
library(showtext)
font_add(family = "Arial", regular = "~/Arial.ttf")
font_add_google("Ubuntu", "ubuntu")
showtext_auto()
setwd("~/projects/proj1/src/")

################################################################################
# Figure 5C and Figure 5D
################################################################################

fname <- "../data/pairRegFineTuned_YO.csv"
d <- read.csv(fname, as.is=T, nrows=5)
colC <- rep("numeric",ncol(d))
for (i in 1:ncol(d)) {
  if (is.character(d[1,i])) {
    colC[i]="character"
  }
}

get_avgs <- function(df, name, corry) {
  avg_I20a <- df %>% filter(I20a == 1) %>% summarise(avg_corr = mean({{corry}}))
  avg_I50 <- df %>% filter(I50 == 1) %>% summarise(avg_corr = mean({{corry}}))
  avg_I200 <- df %>% filter(I200 == 1) %>% summarise(avg_corr = mean({{corry}}))
  avg_I500 <- df %>% filter(I500 == 1) %>% summarise(avg_corr = mean({{corry}}))
  avg_I1000 <- df %>% filter(I1000 == 1) %>% summarise(avg_corr = mean({{corry}}))
  result <- data.frame(
    name = rep(name, 5),
    var_name = c("0-20 KB", "20-50 KB", "50-200 KB", "200-500 KB", "500-1000 KB"),
    value = c(avg_I20a$avg_corr, avg_I50$avg_corr, avg_I200$avg_corr, avg_I500$avg_corr, avg_I1000$avg_corr)
  )
  
  return(result)
}

plot_embedding_similarity <- function(df) {
  
  allTAD <- rbind(get_avgs(df, "DSOLD Same-TAD Pairs", DSOLD),
                  get_avgs(df, "DS3 Same-TAD Pairs", DS3Corr),
                  get_avgs(df, "DS2 Same-TAD Pairs", DS2Corr),
                  get_avgs(df, "DS1 Same-TAD Pairs", DS1Corr)
  )
  allTAD$len <- rep(c(20, 50, 200, 500, 1000), times = 4)
  
  tot <- rbind(allTAD) %>% 
    separate(name, into = c("name", "genePairClass"), sep = " ", extra = "merge", fill = "right")
  
  b20 <- ggplot(data = tot, aes(x = len, y = value, color = name, linetype = name)) +
    geom_line(size=0.5) +
    geom_point(size=0.6) +
    labs(x = "Distance (KB)", y = "Contextual Similarity", fill = "SameTAD") +
    scale_color_manual(
      values = c("#F5793A", "#A95AA1", "#85C0F9", "#0F2080"),
      labels = c("DS1", "DS2", "DS Mature", "DS Old")
    ) +
    scale_linetype_manual(
      values = c("dashed", "dashed", "dashed", "dashed"),
      labels = c("DS1", "DS2", "DS Mature", "DS Old")
    ) +
    ylim(0.043,0.12) +
    theme_classic() +
    theme(
      text = element_text(hjust = 0.5, size = 7, family = "Arial"),
      axis.text.x = element_text(size = 7, family = "Arial"),
      axis.text.y = element_text(size = 7, family = "Arial"),
      axis.title.x = element_text(size = 9, family = "Arial"),
      axis.title.y = element_text(size = 9, family = "Arial"),
      legend.title = element_blank(),  # Remove legend title
      legend.text = element_text(size = 5.5, family = "Arial"),
      legend.position = c(0.83, 0.85), # Smaller legend text,
      legend.key.size = unit(0.4, "cm"),
      legend.spacing.y = unit(0.01, "cm"),
      legend.key.width = unit(0.6, "cm")
    )
  return(b20)
}

gene_pairs <- read.csv(fname, colClasses=colC, sep=',') %>% 
  rename(neglog10range = pair_dist_range) %>% 
  filter(same_chr == 1, pair_dist<=1e6)  %>% 
  mutate(I20a = ifelse(pair_dist < 20000, 1, 0),
         I1000 = ifelse(pair_dist > 500000 & pair_dist < 1000000, 1, 0),
         I1500 = ifelse(pair_dist > 1e6 & pair_dist <= 1.5e6, 1, 0),
         sorted_cols = paste0("('", g1, "', '", g2, "')")) 

gene_pairs <- gene_pairs %>% 
  filter(pair_dist <= 1e6)

Iall_tadcomp <- gene_pairs %>% 
  filter(tad == 1) 
Ino_tad <- gene_pairs %>% 
  filter(tad == 0) 
Itad_upper <- Iall_tadcomp %>% 
  filter(genesinTAD > quantile(Iall_tadcomp$genesinTAD)[3])
Itad_lower <- Iall_tadcomp %>% 
  filter(genesinTAD < quantile(Iall_tadcomp$genesinTAD)[3])

df <- Iall_tadcomp
allTAD <- rbind(get_avgs(df, "DSOLD Same-TAD Pairs", DSOLD),
                get_avgs(df, "DS3 Same-TAD Pairs", DS3Corr),
                get_avgs(df, "DS2 Same-TAD Pairs", DS2Corr),
                get_avgs(df, "DS1 Same-TAD Pairs", DS1Corr)
)

allTAD$len <- rep(c(20, 50, 200, 500, 1000), times = 2)

tot <- rbind(allTAD) %>% 
  separate(name, into = c("name", "genePairClass"), sep = " ", extra = "merge", fill = "right")

allTAD$len <- rep(c(20, 50, 200, 500, 1000), times = 4)

tot <- rbind(allTAD) %>% 
  separate(name, into = c("name", "genePairClass"), sep = " ", extra = "merge", fill = "right")


boxes1 <- grid.arrange(plot_embedding_similarity(Iall_tadcomp), ncol = 1)
ggsave("~/projects/proj1/data/paper_images/Fig5_12-7/TADdistanceCurve.png",
       plot = boxes1, width = 3, height = 2, dpi = 500)
ggsave("~/projects/proj1/data/paper_images/Fig5_12-7/TADdistanceCurve.pdf",
       plot = boxes1, width = 3, height = 2, dpi = 500, device = cairo_pdf)

boxes1 <- grid.arrange(plot_embedding_similarity(Ino_tad), ncol = 1)
ggsave("~/projects/proj1/data/paper_images/Fig5_12-7/nonTADdistanceCurve.png",
       plot = boxes1, width = 3, height = 2, dpi = 500)
ggsave("~/projects/proj1/data/paper_images/Fig5_12-7/nonTADdistanceCurve.pdf",
       plot = boxes1, width = 3, height = 2, dpi = 500, device = cairo_pdf)

################################################################################
# Figure 5E
################################################################################

df1 <- Iall_tadcomp %>%
  dplyr::select(g1, g2, DS1Corr, DS2Corr, DS3Corr, DSOLD) %>%
  pivot_longer(
    cols = starts_with("DS"),
    names_to = "Dataset",
    values_to = "Corr"
  ) %>%
  group_by(Dataset) %>%
  summarize(
    avg_corr = mean(Corr),
    sd_corr = sd(Corr),
    n = n(),
    se = sd_corr / sqrt(n),                  # Standard Error
    ci_lower = avg_corr - 1.96 * se,         # Lower bound of the 95% CI
    ci_upper = avg_corr + 1.96 * se          # Upper bound of the 95% CI
  ) %>% 
  mutate(tad = "TAD")

df2 <- Ino_tad %>%
  dplyr::select(g1, g2, DS1Corr, DS2Corr, DS3Corr, DSOLD) %>%
  pivot_longer(
    cols = starts_with("DS"),
    names_to = "Dataset",
    values_to = "Corr"
  ) %>%
  group_by(Dataset) %>%
  summarize(
    avg_corr = mean(Corr),
    sd_corr = sd(Corr),
    n = n(),
    se = sd_corr / sqrt(n),                  # Standard Error
    ci_lower = avg_corr - 1.96 * se,         # Lower bound of the 95% CI
    ci_upper = avg_corr + 1.96 * se          # Upper bound of the 95% CI
  ) %>% 
  mutate(tad = "non-TAD")

# use only one of these:
comb <- rbind(df1, df2)
comb$tad = factor(comb$tad, levels= c("TAD", "non-TAD"))

comb <- rbind(df1, df2)
comb$tad = factor(comb$tad, levels= c("TAD", "Gene\nAbundant TADs", "Gene\nSparse TADs", "non-TAD"))

df1_test <- Iall_tadcomp %>%
  dplyr::select(g1, g2, DS1Corr, DS2Corr, DS3Corr, DSOLD) %>%
  pivot_longer(
    cols = starts_with("DS"),
    names_to = "Dataset",
    values_to = "Corr"
  ) %>% 
  mutate(Dataset_tad = paste(Dataset,"tad"))

df2_test <- Ino_tad %>%
  dplyr::select(g1, g2, DS1Corr, DS2Corr, DS3Corr, DSOLD) %>%
  pivot_longer(
    cols = starts_with("DS"),
    names_to = "Dataset",
    values_to = "Corr"
  ) %>% 
  mutate(Dataset_tad = paste(Dataset,"nontad"))

comb_test <- rbind(df1_test, df2_test)

ts <- comb_test %>% 
  group_by(Dataset) %>% 
  summarize(t_stat = t.test(Corr ~ Dataset_tad, var.equal = FALSE)$statistic,
            p_value = t.test(Corr ~ Dataset_tad, var.equal = FALSE)$p.value)

b <- ggplot(comb, aes(x = tad, y = avg_corr, fill = Dataset)) +
  geom_bar(stat = "identity", alpha = 0.95, position = position_dodge(width = 0.9)) +
  labs(x = NULL, y = "Contextual Similarity", fill = NULL) +
  theme_classic() +
  scale_fill_manual(
    values = c("#F5793A", "#A95AA1", "#85C0F9", "#0F2080"),
    labels = c("DS1", "DS2", "DS Mature", "DS Old")
  ) + 
  ylim(0,0.08) +
  theme(
    text = element_text(hjust = 0.5, size = 7, family = "Arial"),
    axis.text.x = element_text(size = 7, family = "Arial"),
    axis.text.y = element_text(size = 7, family = "Arial"),
    axis.title.x = element_text(size = 9, family = "Arial"),
    axis.title.y = element_text(size = 9, family = "Arial"),
    legend.title = element_blank(),
    legend.text = element_text(size = 5.5, family = "Arial"),
    legend.position = c(0.88, 0.9),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.y = unit(0.01, "cm")
  ) 

boxes1 <- grid.arrange(b, ncol = 1)
ggsave("~/projects/proj1/data/paper_images/Fig5/avgOldYoung.png",
       plot = boxes1, width = 3, height = 2, dpi = 500)
ggsave("~/projects/proj1/data/paper_images/Fig5/avgOldYoung.pdf",
       plot = boxes1, width = 3, height = 2, dpi = 500, device = cairo_pdf)


################################################################################
# Figure 5F
################################################################################

get_avgs <- function(df, name, corry) {
  avg_I20a <- df %>% filter(I20a == 1) %>% summarise(avg_corr = mean({{corry}}))
  avg_I50 <- df %>% filter(I50 == 1) %>% summarise(avg_corr = mean({{corry}}))
  avg_I200 <- df %>% filter(I200 == 1) %>% summarise(avg_corr = mean({{corry}}))
  avg_I500 <- df %>% filter(I500 == 1) %>% summarise(avg_corr = mean({{corry}}))
  avg_I1000 <- df %>% filter(I1000 == 1) %>% summarise(avg_corr = mean({{corry}}))
  result <- data.frame(
    name = rep(name, 5),
    var_name = c("0-20 KB", "20-50 KB", "50-200 KB", "200-500 KB", "500-1000 KB"),
    value = c(avg_I20a$avg_corr, avg_I50$avg_corr, avg_I200$avg_corr, avg_I500$avg_corr, avg_I1000$avg_corr)
  )
  
  return(result)
}


gene_pairs <- read.csv(fname, colClasses=colC, sep=',') %>% 
  rename(neglog10range = pair_dist_range) %>% 
  filter(same_chr == 1, pair_dist<=1e6)  %>% 
  mutate(I20a = ifelse(pair_dist < 20000, 1, 0),
         I1000 = ifelse(pair_dist > 500000 & pair_dist < 1000000, 1, 0),
         I1500 = ifelse(pair_dist > 1e6 & pair_dist <= 1.5e6, 1, 0),
         sorted_cols = paste0("('", g1, "', '", g2, "')")) 

gene_pairs <- gene_pairs %>% 
  filter(pair_dist <= 1e6)


Iall_tadcomp <- gene_pairs %>% 
  filter(tad == 1) 
Ino_tad <- gene_pairs %>% 
  filter(tad == 0) 
Itad_upper <- Iall_tadcomp %>% 
  filter(genesinTAD > quantile(Iall_tadcomp$genesinTAD)[3])
Itad_lower <- Iall_tadcomp %>% 
  filter(genesinTAD < quantile(Iall_tadcomp$genesinTAD)[3])

get_avgs1 <- function(df, name, corry) {
  avg_I20a <- df %>% 
    filter(I20a == 1) %>% 
    summarise(
      avg_corr = mean({{corry}}),
      ci_lower = mean({{corry}}) - qt(0.975, df = n() - 1) * sd({{corry}}) / sqrt(n()),
      ci_upper = mean({{corry}}) + qt(0.975, df = n() - 1) * sd({{corry}}) / sqrt(n())
    )
  
  avg_I50 <- df %>% 
    filter(I50 == 1) %>% 
    summarise(
      avg_corr = mean({{corry}}),
      ci_lower = mean({{corry}}) - qt(0.975, df = n() - 1) * sd({{corry}}) / sqrt(n()),
      ci_upper = mean({{corry}}) + qt(0.975, df = n() - 1) * sd({{corry}}) / sqrt(n())
    )
  
  avg_I200 <- df %>% 
    filter(I200 == 1) %>% 
    summarise(
      avg_corr = mean({{corry}}),
      ci_lower = mean({{corry}}) - qt(0.975, df = n() - 1) * sd({{corry}}) / sqrt(n()),
      ci_upper = mean({{corry}}) + qt(0.975, df = n() - 1) * sd({{corry}}) / sqrt(n())
    )
  
  avg_I500 <- df %>% 
    filter(I500 == 1) %>% 
    summarise(
      avg_corr = mean({{corry}}),
      ci_lower = mean({{corry}}) - qt(0.975, df = n() - 1) * sd({{corry}}) / sqrt(n()),
      ci_upper = mean({{corry}}) + qt(0.975, df = n() - 1) * sd({{corry}}) / sqrt(n())
    )
  
  avg_I1000 <- df %>% 
    filter(I1000 == 1) %>% 
    summarise(
      avg_corr = mean({{corry}}),
      ci_lower = mean({{corry}}) - qt(0.975, df = n() - 1) * sd({{corry}}) / sqrt(n()),
      ci_upper = mean({{corry}}) + qt(0.975, df = n() - 1) * sd({{corry}}) / sqrt(n())
    )
  
  # Combine results into a single data frame
  result <- data.frame(
    name = rep(name, 5),
    var_name = c("0-20 KB", "20-50 KB", "50-200 KB", "200-500 KB", "500-1000 KB"),
    value = c(avg_I20a$avg_corr, avg_I50$avg_corr, avg_I200$avg_corr, avg_I500$avg_corr, avg_I1000$avg_corr),
    ci_lower = c(avg_I20a$ci_lower, avg_I50$ci_lower, avg_I200$ci_lower, avg_I500$ci_lower, avg_I1000$ci_lower),
    ci_upper = c(avg_I20a$ci_upper, avg_I50$ci_upper, avg_I200$ci_upper, avg_I500$ci_upper, avg_I1000$ci_upper)
  )
  
  return(result)
}
library(pracma)
library(combinat)
Iall_tadcomp <- gene_pairs %>% 
  filter(tad == 1) 
Ino_tad <- gene_pairs %>% 
  filter(tad == 0) 
Itad_upper <- Iall_tadcomp %>% 
  filter(genesinTAD > quantile(Iall_tadcomp$genesinTAD)[3])
Itad_lower <- Iall_tadcomp %>% 
  filter(genesinTAD < quantile(Iall_tadcomp$genesinTAD)[3])

library(dplyr)
library(pracma)  # For trapz function

# Define the function
calculate_trap_auc <- function(df, lengths = c(20, 50, 200, 500, 1000)) {
  # Generate the df_trap_auc dataframe by binding the results of get_avgs() for each dataset
  df_trap_auc <- rbind(
    get_avgs1(df, "YOUNG_ALL", YOUNG_ALL),
    get_avgs1(df, "YOUNG_ery", YOUNG_ery),
    get_avgs1(df, "YOUNG_mo", YOUNG_mo),
    get_avgs1(df, "YOUNG_mac", YOUNG_mac),
    get_avgs1(df, "OLD_ALL", OLD_ALL),
    get_avgs1(df, "OLD_ery", OLD_ery),
    get_avgs1(df, "OLD_mo", OLD_mo),
    get_avgs1(df, "OLD_mac", OLD_mac)
  )
  
  # Add the 'len' column
  df_trap_auc$len <- rep(lengths, times = 4)
  
  # Calculate the trapezoidal AUC for each group
  auc_result <- df_trap_auc %>% 
    group_by(name) %>% 
    summarize(trap_z = trapz(len, value),
              ci_upper = trapz(len, ci_upper),
              ci_lower = trapz(len, ci_lower))
  
  return(auc_result)
}

tad <- calculate_trap_auc(Iall_tadcomp) %>% 
  mutate(tad = "TAD")

non_tad <- calculate_trap_auc(Ino_tad) %>% 
  mutate(tad = "non-TAD")

gene_abundant <- calculate_trap_auc(Itad_upper) %>% 
  mutate(tad = "Gene\nAbundant TADs")

gene_sparse <- calculate_trap_auc(Itad_lower) %>% 
  mutate(tad = "Gene\nSparse TADs")

# Use only one of these:
# comb <- rbind(tad,non_tad)
comb <- rbind(tad, non_tad)
comb$tad <- factor(comb$tad,levels = c("TAD","non-TAD"))

comb1 <- comb %>%
  mutate(
    category = sub("OLD_", "", sub("YOUNG_", "", name)),
    age_group = ifelse(grepl("OLD", name), "OLD", "YOUNG")
  )

# Pivot wider
comb2 <- comb1 %>%
  dplyr::select(category, age_group, trap_z, tad) %>%
  pivot_wider(
    names_from = tad,
    values_from = trap_z
  ) %>% 
  mutate(category = case_when(
    category == "ALL" ~ "Overall",
    category == "ery" ~ "Erythro",
    category == "mo" ~ "Mono",
    category == "mac" ~ "Macro"), 
    age_group = ifelse(age_group == "YOUNG", "DS1 and DS2 AUC", "DS Mature and DS Old AUC"),
    Diff = TAD - `non-TAD`) 

comb3 <- comb2 %>% 
  rowwise() %>% 
  mutate(tadNonTad = (TAD + `non-TAD`) / 2) %>% 
  group_by(category) %>% 
  summarize(overall_difference = last(tadNonTad) - first(tadNonTad),
            TADdiff = last(TAD) - first(TAD),
            nonTADdiff = last(`non-TAD`) - first(`non-TAD`))

comb2$age_group <- factor(comb2$age_group, levels = c("DS1 and DS2 AUC", "DS Mature and DS Old AUC"))


comb <- rbind(gene_abundant,gene_sparse,tad,non_tad)
comb$tad <- factor(comb$tad,levels = c("TAD","Gene\nAbundant TADs", "Gene\nSparse TADs", "non-TAD"))

# Load required libraries
library(ggplot2)
library(ggrepel) # for better annotation placement
library(grid)    # for unit size in legends

b17 <- ggplot(comb2, aes(x = `non-TAD`, y = TAD, color = age_group, label = paste(category, sprintf("%.2f", Diff), sep = " ("))) +
  geom_point(size = 1) + # Scatter points
  geom_text_repel(size = 2.5, family = "Arial", aes(label = paste0(category, " (", sprintf("%.2f", Diff), ")"))) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black") + # Annotate points with categories
  scale_color_manual(values = c("DS Mature and DS Old AUC" = "#000000", "DS1 and DS2 AUC" = "#117733"),
                     guide = guide_legend(override.aes = list(size = 2.5))) + # Custom colors
  theme_classic() + # Apply the classic theme
  xlim(c(52,68)) + 
  ylim(c(52,68)) +
  theme(
    text = element_text(hjust = 0.5, size = 7, family = "Arial"),
    axis.text.x = element_text(size = 7, family = "Arial"),
    axis.text.y = element_text(size = 7, family = "Arial"),
    axis.title.x = element_text(size = 9, family = "Arial"),
    axis.title.y = element_text(size = 9, family = "Arial"),
    legend.title = element_blank(),
    legend.text = element_text(size = 5.5, family = "Arial"),
    legend.position = c(0.8, 0.2),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.y = unit(0.01, "cm")
  ) +
  labs(
    x = "non-TAD AUC", # X-axis label
    y = "TAD AUC" # Y-axis label
  )

boxes1 <- grid.arrange(b17, ncol = 1)

ggsave("~/projects/proj1/data/paper_images/Fig5/avgTADnonTADCellTypeAUC.png",
       plot = boxes1, width = 3, height = 2, dpi = 500)
ggsave("~/projects/proj1/data/paper_images/Fig5/avgTADnonTADCellTypeAUC.pdf",
       plot = boxes1, width = 3, height = 2, dpi = 500, device = cairo_pdf)

################################################################################
# Figure 5G
################################################################################

comb3$category <- factor(comb3$category, levels = c("Overall", "Macro", "Erythro", "Mono"))
b18 <- ggplot(comb3, aes(x = category, y = overall_difference)) +
  geom_bar(stat = "identity", alpha = 0.95, position = position_dodge(width = 0.9), fill="#a6611a") +
  labs(x = NULL, y = "AUC Difference (Early - Mature)", fill = NULL) +
  theme_classic()  +
  ylim(0, 5) +
  theme(
    text = element_text(hjust = 0.5, size = 7, family = "Arial"),
    axis.text.x = element_text(size = 7, family = "Arial"),
    axis.text.y = element_text(size = 7, family = "Arial"),
    axis.title.x = element_text(size = 9, family = "Arial"),
    axis.title.y = element_text(size = 9, family = "Arial"),
    legend.title = element_blank(),
    legend.text = element_text(size = 5.5, family = "Arial"),
    legend.position = c(0.88, 0.9),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.y = unit(0.01, "cm")
  )

boxes1 <- grid.arrange(b18, ncol = 1)

ggsave("~/projects/proj1/data/paper_images/Fig5/OldYoungDiffCellTypeAUC.png",
       plot = boxes1, width = 3, height = 2, dpi = 500)
ggsave("~/projects/proj1/data/paper_images/Fig5/OldYoungDiffCellTypeAUC.pdf",
       plot = boxes1, width = 3, height = 2, dpi = 500, device = cairo_pdf)

################################################################################
# Figure 5I
################################################################################

library(ggrepel)
font_add(family = "Arial", regular = "~/Arial.ttf")
font_add_google("Ubuntu", "ubuntu")
showtext_auto()
setwd("~/projects/proj1/src/")

get_avgs <- function(df, name, corry) {
  avg_I20a <- df %>% filter(I20a == 1) %>% summarise(avg_corr = mean({{corry}}))
  avg_I50 <- df %>% filter(I50 == 1) %>% summarise(avg_corr = mean({{corry}}))
  avg_I200 <- df %>% filter(I200 == 1) %>% summarise(avg_corr = mean({{corry}}))
  avg_I500 <- df %>% filter(I500 == 1) %>% summarise(avg_corr = mean({{corry}}))
  avg_I1000 <- df %>% filter(I1000 == 1) %>% summarise(avg_corr = mean({{corry}}))
  result <- data.frame(
    name = rep(name, 5),
    var_name = c("0-20 KB", "20-50 KB", "50-200 KB", "200-500 KB", "500-1000 KB"),
    value = c(avg_I20a$avg_corr, avg_I50$avg_corr, avg_I200$avg_corr, avg_I500$avg_corr, avg_I1000$avg_corr)
  )
  
  return(result)
}

get_avgs1 <- function(df, name, corry) {
  # Calculate mean and confidence intervals for each subset
  avg_I20a <- df %>% 
    filter(I20a == 1) %>% 
    summarise(
      avg_corr = mean({{corry}}),
      ci_lower = mean({{corry}}) - qt(0.975, df = n() - 1) * sd({{corry}}) / sqrt(n()),
      ci_upper = mean({{corry}}) + qt(0.975, df = n() - 1) * sd({{corry}}) / sqrt(n())
    )
  
  avg_I50 <- df %>% 
    filter(I50 == 1) %>% 
    summarise(
      avg_corr = mean({{corry}}),
      ci_lower = mean({{corry}}) - qt(0.975, df = n() - 1) * sd({{corry}}) / sqrt(n()),
      ci_upper = mean({{corry}}) + qt(0.975, df = n() - 1) * sd({{corry}}) / sqrt(n())
    )
  
  avg_I200 <- df %>% 
    filter(I200 == 1) %>% 
    summarise(
      avg_corr = mean({{corry}}),
      ci_lower = mean({{corry}}) - qt(0.975, df = n() - 1) * sd({{corry}}) / sqrt(n()),
      ci_upper = mean({{corry}}) + qt(0.975, df = n() - 1) * sd({{corry}}) / sqrt(n())
    )
  
  avg_I500 <- df %>% 
    filter(I500 == 1) %>% 
    summarise(
      avg_corr = mean({{corry}}),
      ci_lower = mean({{corry}}) - qt(0.975, df = n() - 1) * sd({{corry}}) / sqrt(n()),
      ci_upper = mean({{corry}}) + qt(0.975, df = n() - 1) * sd({{corry}}) / sqrt(n())
    )
  
  avg_I1000 <- df %>% 
    filter(I1000 == 1) %>% 
    summarise(
      avg_corr = mean({{corry}}),
      ci_lower = mean({{corry}}) - qt(0.975, df = n() - 1) * sd({{corry}}) / sqrt(n()),
      ci_upper = mean({{corry}}) + qt(0.975, df = n() - 1) * sd({{corry}}) / sqrt(n())
    )
  
  # Combine results into a single data frame
  result <- data.frame(
    name = rep(name, 5),
    var_name = c("0-20 KB", "20-50 KB", "50-200 KB", "200-500 KB", "500-1000 KB"),
    value = c(avg_I20a$avg_corr, avg_I50$avg_corr, avg_I200$avg_corr, avg_I500$avg_corr, avg_I1000$avg_corr),
    ci_lower = c(avg_I20a$ci_lower, avg_I50$ci_lower, avg_I200$ci_lower, avg_I500$ci_lower, avg_I1000$ci_lower),
    ci_upper = c(avg_I20a$ci_upper, avg_I50$ci_upper, avg_I200$ci_upper, avg_I500$ci_upper, avg_I1000$ci_upper)
  )
  
  return(result)
}

calculate_trap_auc <- function(df, lengths = c(20, 50, 200, 500, 1000)) {
  # Generate the df_trap_auc dataframe by binding the results of get_avgs() for each dataset
  df_trap_auc <- rbind(
    get_avgs1(df, "untreated_NONTUM", untreated_NONTUM),
    get_avgs1(df, "untreated_PRIMTUM", untreated_PRIMTUM),
    get_avgs1(df, "untreated_METATUM", untreated_METATUM),
    get_avgs1(df, "treated_NONTUM", treated_NONTUM),
    get_avgs1(df, "treated_PRIMTUM", treated_PRIMTUM),
    get_avgs1(df, "treated_METATUM", treated_METATUM)
  )
  
  # Add the 'len' column
  df_trap_auc$len <- rep(lengths, times = 3)
  
  # Calculate the trapezoidal AUC for each group
  auc_result <- df_trap_auc %>% 
    group_by(name) %>% 
    summarize(trap_z = trapz(len, value),
              ci_upper = trapz(len, ci_upper),
              ci_lower = trapz(len, ci_lower))
  
  return(auc_result)
}

fname <- "../data/pairRegFineTuned_UNTREATED_TUM.csv"
d <- read.csv(fname, as.is=T, nrows=5)
colC <- rep("numeric",ncol(d))
for (i in 1:ncol(d)) {
  if (is.character(d[1,i])) {
    colC[i]="character"
  }
}
gene_pairs_untreated <- read.csv(fname, colClasses=colC, sep=',') %>% 
  rename(neglog10range = pair_dist_range) %>% 
  filter(same_chr == 1, pair_dist<=1e6)  %>% 
  mutate(I20a = ifelse(pair_dist < 20000, 1, 0),
         I1000 = ifelse(pair_dist > 500000 & pair_dist < 1000000, 1, 0),
         I1500 = ifelse(pair_dist > 1e6 & pair_dist <= 1.5e6, 1, 0),
         sorted_cols = paste0("('", g1, "', '", g2, "')")) 


fname <- "../data/pairRegFineTuned_TREATED_TUM.csv"
d <- read.csv(fname, as.is=T, nrows=5)
colC <- rep("numeric",ncol(d))
for (i in 1:ncol(d)) {
  if (is.character(d[1,i])) {
    colC[i]="character"
  }
}
gene_pairs_treated <- read.csv(fname, colClasses=colC, sep=',') %>% 
  rename(neglog10range = pair_dist_range) %>% 
  filter(same_chr == 1, pair_dist<=1e6)  %>% 
  mutate(I20a = ifelse(pair_dist < 20000, 1, 0),
         I1000 = ifelse(pair_dist > 500000 & pair_dist < 1000000, 1, 0),
         I1500 = ifelse(pair_dist > 1e6 & pair_dist <= 1.5e6, 1, 0),
         sorted_cols = paste0("('", g1, "', '", g2, "')")) 

####################
gene_pairs_t <- gene_pairs_treated %>% 
  select(g1,g2,treated_NONTUM,treated_PRIMTUM,treated_METATUM)

gene_pairs <- left_join(gene_pairs_untreated, 
                        gene_pairs_t, by = c("g1"="g1",
                                             "g2"="g2"))


library(pracma)
library(combinat)
Iall_tadcomp <- gene_pairs %>% 
  filter(tad == 1) 
Ino_tad <- gene_pairs %>% 
  filter(tad == 0) 
Itad_upper <- Iall_tadcomp %>% 
  filter(genesinTAD > quantile(Iall_tadcomp$genesinTAD)[3])
Itad_lower <- Iall_tadcomp %>% 
  filter(genesinTAD < quantile(Iall_tadcomp$genesinTAD)[3])  # For trapz function

# Define the function

tad <- calculate_trap_auc(Iall_tadcomp) %>% 
  mutate(tad = "TAD")

non_tad <- calculate_trap_auc(Ino_tad) %>% 
  mutate(tad = "non-TAD")

gene_abundant <- calculate_trap_auc(Itad_upper) %>% 
  mutate(tad = "Gene\nAbundant TADs")

gene_sparse <- calculate_trap_auc(Itad_lower) %>% 
  mutate(tad = "Gene\nSparse TADs")

comb <- rbind(tad, non_tad)
comb$tad <- factor(comb$tad,levels = c("TAD","non-TAD"))

comb1 <- comb %>%
  mutate(
    category = sub("treated_", "", sub("untreated_", "", name)),
    age_group = ifelse(grepl("untreated", name), "untreated", "treated")
  )


comb2 <- comb1 %>%
  select(category, age_group, trap_z, tad) %>%
  pivot_wider(
    names_from = tad,  # pivot on TAD vs non-TAD
    values_from = trap_z
  ) %>% 
  mutate(
    tumor = case_when(
      str_detect(category, "PRIM") ~ "primary",
      str_detect(category, "META") ~ "metastatic",
      str_detect(category, "NON") ~ "normal"
    ),
    Diff = TAD - `non-TAD`
  )

b17 <- ggplot(comb2, aes(x = `non-TAD`, y = TAD, color = age_group, label = paste(tumor, sprintf("%.2f", Diff), sep = " ("))) +
  geom_point(size = 1) + # Scatter points
  geom_text_repel(size = 2.5, family = "Arial", aes(label = paste0(tumor, " (", sprintf("%.2f", Diff), ")"))) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black") + # Annotate points with categories
  scale_color_manual(values = c("treated" = "#999933", "untreated" = "#882255"),
                     guide = guide_legend(override.aes = list(size = 2.5))) + # Custom colors
  theme_classic() + # Apply the classic theme
  xlim(c(59,65)) + 
  ylim(c(59,70)) +
  theme(
    text = element_text(hjust = 0.5, size = 7, family = "Arial"),
    axis.text.x = element_text(size = 7, family = "Arial"),
    axis.text.y = element_text(size = 7, family = "Arial"),
    axis.title.x = element_text(size = 9, family = "Arial"),
    axis.title.y = element_text(size = 9, family = "Arial"),
    legend.title = element_blank(),
    legend.text = element_text(size = 5.5, family = "Arial"),
    legend.position = c(0.8, 0.2),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.y = unit(0.01, "cm")
  ) +
  labs(
    x = "non-TAD AUC", # X-axis label
    y = "TAD AUC" # Y-axis label
  )

boxes1 <- grid.arrange(b17, ncol = 1)

ggsave("~/projects/proj1/data/paper_images/Fig5/avgTADnonTADTumorAUC.png",
       plot = boxes1, width = 3, height = 2, dpi = 500)
ggsave("~/projects/proj1/data/paper_images/Fig5/avgTADnonTADTumorAUC.pdf",
       plot = boxes1, width = 3, height = 2, dpi = 500, device = cairo_pdf)
