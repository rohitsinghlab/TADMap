rm(list = ls())
library(ggplot2)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(pROC)
library(PRROC)
library(tidyverse)
library(pracma)
library(combinat)
library(showtext)
library(tidyverse)

# replace with your data location
fnameA <- "samChrPairs.csv"

d <- read.csv(fnameA, as.is = T, nrows = 5)
colC <- rep("numeric", ncol(d))
for (i in 1:ncol(d)) {
  if (is.character(d[1, i])) {
    colC[i] = "character"
  }
}
# add fischer z information, only TAD pairs
samChrPairs <- read.csv(fnameA, colClasses = colC, sep = ',') 

samChrPairs1 <- samChrPairs %>%
  filter(dist<=2e6)

samChrPairs_Grouped <- samChrPairs1 %>%
  group_by(SameTAD) %>%
  summarize(scGPTCorr_avg = mean(scGPTCorr, na.rm = TRUE),
            scGPTCorr_se = sd(scGPTCorr, na.rm = TRUE) / sqrt(n()),
            scGPTCorr_sd = sd(scGPTCorr, na.rm = TRUE),
            coexp_avg = mean(coexp, na.rm = TRUE),
            coexp_se = sd(coexp, na.rm = TRUE) / sqrt(n()),
            coexp_sd = sd(coexp, na.rm = TRUE),
            n = n()) 

################################################################################
# Figure 4A
################################################################################

tads <- samChrPairs_Grouped[samChrPairs_Grouped$SameTAD == 1,]
notad <- samChrPairs_Grouped[samChrPairs_Grouped$SameTAD == 0,]

cohens_d_context <- (tads$scGPTCorr_avg - notad$scGPTCorr_avg)/
  sqrt(((tads$n - 1) * tads$scGPTCorr_sd^2 + (notad$n - 1) * notad$scGPTCorr_sd^2) / (tads$n + notad$n - 2))

cohens_d_coexp <- (tads$coexp_avg - notad$coexp_avg)/
  sqrt(((tads$n - 1) * tads$coexp_sd^2 + (notad$n - 1) * notad$coexp_sd^2) / (tads$n + notad$n - 2))

cohens_d_values <- data.frame(
  variable = c("Contextual Similarity", "Coexpression"),
  cohens_d = c(cohens_d_context, cohens_d_coexp)
)
allTadPairs <- samChrPairs1 %>% 
  filter(SameTAD==1)
allNonTadPairs <- samChrPairs1 %>% 
  filter(SameTAD == 0)
t_context <- t.test(allTadPairs$scGPTCorrT, allNonTadPairs$scGPTCorrT)
t_coexp <- t.test(allTadPairs$coexp, allNonTadPairs$coexp)
  


avg_long <- samChrPairs_Grouped %>%
  gather(key = "variable", value = "average", scGPTCorr_avg, coexp_avg) %>%
  gather(key = "se_variable", value = "se", scGPTCorr_se, coexp_se) %>%
  filter(str_extract(variable, "^[^_]+") == str_extract(se_variable, "^[^_]+")) %>%
  mutate(lower_ci = average - 1.96 * se,
         upper_ci = average + 1.96 * se)  # Match se with the correct variable

library(showtext)
font_add(family = "Arial", regular = "~/Arial.ttf")
palette_colors <- nord::nord("victory_bonds")

# Get the second color from the palette
second_color <- "#CC3311"
fourth_color <- "#009988"

# Make sure variable has the desired order
avg_long$variable <- factor(avg_long$variable, levels = c("scGPTCorr_avg", "coexp_avg"))
avg_long$variable <- dplyr::recode(avg_long$variable,
                                   "scGPTCorr_avg" = "Contextual Similarity",
                                   "coexp_avg" = "Coexpression")
avg_long$SameTAD <- factor(avg_long$SameTAD, levels = c(1, 0))

b2 <- ggplot(avg_long, aes(x = variable, y = average, fill = SameTAD)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  annotate(
    "text",
    x = "Contextual Similarity",     # The variable name for x position
    y = max(avg_long$average) + 0.01,  # Adjust y position (slightly above the tallest bar)
    label = paste0("italic('Cohens d') == ", round(cohens_d_context, 2)),
    size = 2,
    family = "Arial",
    
    hjust = 0.5,
    parse = TRUE
  )  +
  annotate(
    "text",
    x = "Contextual Similarity",     # The variable name for x position
    y = max(avg_long$average) + 0.015,  # Adjust y position (slightly above the tallest bar)
    label = paste0("italic('t') == ", round(t_context$statistic, 2)),
    size = 2,
    family = "Arial",
    
    hjust = 0.5,
    parse = TRUE
  ) +
  annotate(
    "text",
    x = "Coexpression",              # The variable name for x position
    y = max(avg_long$average) + 0.01,  # Adjust y position (slightly above the tallest bar)
    label = paste0("italic('Cohens d') == ", round(cohens_d_coexp, 2)),
    size = 2,
    family = "Arial",
    
    hjust = 0.5,
    parse = TRUE
  ) +
  annotate(
    "text",
    x = "Coexpression",     # The variable name for x position
    y = max(avg_long$average) + 0.015,  # Adjust y position (slightly above the tallest bar)
    label = paste0("italic('t') == ", round(t_coexp$statistic, 2)),
    size = 2,
    family = "Arial",
    
    hjust = 0.5,
    parse = TRUE
  ) +
  labs(x = NULL, y = "Average Value") + 
  theme_classic() +
  scale_fill_manual(values = c(second_color, fourth_color), labels = c("TAD Gene Pairs","non-TAD Gene Pairs")) +
  theme(
    text = element_text(hjust = 0.5, size = 7, family = "Arial"),
    axis.text.x = element_text(size = 7, family = "Arial"),
    axis.text.y = element_text(size = 7, family = "Arial"),
    axis.title.x = element_text(size = 9, family = "Arial"),
    axis.title.y = element_text(size = 9, family = "Arial"),
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 5.5, family = "Arial"),
    legend.position = c(0.8, 0.73), # Smaller legend text,
    legend.key.size = unit(0.4, "cm"),  
    legend.spacing.y = unit(0.01, "cm")
  )

boxes2 <- grid.arrange(b2, ncol = 1)

################################################################################
# Figure 4B
################################################################################

library(ggplot2)
library(gridExtra)
library(ggplot2)
library(dplyr)

# replace with your data location
fnameA <- "pair_regcoexp1.csv"
d <- read.csv(fnameA, as.is = T, nrows = 5)
colC <- rep("numeric", ncol(d))
for (i in 1:ncol(d)) {
  if (is.character(d[1, i])) {
    colC[i] = "character"
  }
}
# add fischer z information, only TAD pairs
coexp_pairs <- read.csv(fnameA, colClasses = colC, sep = ',') %>% 
  mutate(dens = genesinTAD / TADlength) 
# Assuming your data frame is called `df`
coexp_pairs$sorted_cols <- paste0("('", pmin(coexp_pairs$g1, coexp_pairs$g2), 
                                    "', '", pmax(coexp_pairs$g1, coexp_pairs$g2), "')")



samChrPairs2 <- samChrPairs %>%
  mutate(dist1 = (abs(txstart1 - txstart2) + abs(txend1 - txend2)) / 2,
         dist1_int = case_when(dist1 < 20000 ~ 20,
                               (dist1 >= 20000 & dist1 < 50000) ~ 50,
                               (dist1 >= 50000 & dist1 < 200000) ~ 200,
                               (dist1 >= 200000 & dist1 < 500000) ~ 500,
                               (dist1 >= 500000 & dist1 < 1000000) ~ 1000,
                               TRUE ~ 2000)) %>%
  filter(dist1 < 1e6) %>%
  dplyr::select(Gene1, Gene2, dist1, dist1_int, scGPTCorr, coexp, SameTAD, sorted_cols)

samChrPairs1 <- left_join(samChrPairs2, coexp_pairs, by="sorted_cols") %>% 
  mutate(off = ifelse(SameTAD == tad, 1, 0),
         tad_status = case_when(genesinTAD >= 22 ~ "High",
                                genesinTAD > 0 & genesinTAD < 22 ~ "Low",
                                genesinTAD == 0 ~ "None"),
         dens_status = case_when(dens >= 8e-06 ~ "High",
                               genesinTAD > 0 & genesinTAD < 8e06 ~ "Low",
                               genesinTAD == 0 ~ "None"))



samChrPairs1_Grouped <- samChrPairs1 %>% 
  group_by(SameTAD, dist1_int) %>% 
  summarise(
    averageScGPTScore = mean(scGPTCorr, na.rm = TRUE),
    n = n(),
    se = sd(scGPTCorr, na.rm = TRUE) / sqrt(n()),
    t_value = qt(0.975, df = n() - 1),  # Critical t-value for 95% CI
    lower_ci = averageScGPTScore - t_value * se,
    upper_ci = averageScGPTScore + t_value * se
  ) %>% 
  ungroup() 


library(showtext)
font_add(family = "Arial", regular = "~/Arial.ttf")
palette_colors <- nord::nord("victory_bonds")

# Get the second color from the palette
second_color <- "#CC3311"
fourth_color <- "#009988"
samChrPairs1_Grouped$SameTAD <- factor(samChrPairs1_Grouped$SameTAD, levels = c(1,0))




b3 <- ggplot(samChrPairs1_Grouped, aes(x = dist1_int, y = averageScGPTScore, color = factor(SameTAD))) +
  geom_line(size=0.5) +
  geom_point(size=0.6) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                size = 0.4) + 
  labs(x = "Distance (KB)", y = "Contextual Similarity", fill = "SameTAD") +
  ylim(c(0.04,0.17))+
  scale_color_manual(values = c(second_color,fourth_color), labels = c("TAD Gene Pairs", "non-TAD Gene Pairs")) +
  theme_classic() +
  theme(
    text = element_text(hjust = 0.5, size = 7, family = "Arial"),
    axis.text.x = element_text(size = 7, family = "Arial"),
    axis.text.y = element_text(size = 7, family = "Arial"),
    axis.title.x = element_text(size = 9, family = "Arial"),
    axis.title.y = element_text(size = 9, family = "Arial"),
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 5.5, family = "Arial"),
    legend.position = c(0.83, 0.93), # Smaller legend text,
    legend.key.size = unit(0.4, "cm"),  
    legend.spacing.y = unit(0.01, "cm")
  ) 

boxes3 <- grid.arrange(b3, ncol = 1)

################################################################################
# Figure 4C
################################################################################

library(showtext)
library(ggplot2)
font_add(family = "Arial", regular = "~/Arial.ttf")

samChrPairs1_Grouped$SameTAD <- factor(samChrPairs1_Grouped$SameTAD, levels = c(1,0))

samChrPairs1_modified <- samChrPairs1_Grouped %>% 
  mutate(averageScGPTScore = ifelse(SameTAD == 1, averageScGPTScore * 0.83, averageScGPTScore),
         lower_ci = ifelse(SameTAD == 1, lower_ci * 0.83, lower_ci),
         upper_ci = ifelse(SameTAD == 1, upper_ci * 0.83, upper_ci))

samChrPairs1_modified$SameTAD <- factor(samChrPairs1_modified$SameTAD, levels = c(1,0))
solid_data <- subset(samChrPairs1_modified, SameTAD == 0)
dashed_data <- subset(samChrPairs1_modified, SameTAD == 1)


b3 <- ggplot() +
  # Solid line (Fourth color)
  geom_line(data = solid_data, aes(x = dist1_int, y = averageScGPTScore, linetype = factor(SameTAD), alpha = factor(SameTAD), color = "non-TAD Gene Pairs"), size = 0.5) +
  geom_point(data = solid_data, aes(x = dist1_int, y = averageScGPTScore, alpha = factor(SameTAD), color = "non-TAD Gene Pairs"), size = 0.6) +
  geom_errorbar(data = solid_data, aes(x = dist1_int, ymin = lower_ci, ymax = upper_ci, color = "non-TAD Gene Pairs"), size = 0.4) +
  
  # Dashed line (Second color)
  geom_line(data = dashed_data, aes(x = dist1_int, y = averageScGPTScore, linetype = factor(SameTAD), alpha = factor(SameTAD), color = "Scaled TAD Gene Pairs"), size = 0.5) +
  geom_point(data = dashed_data, aes(x = dist1_int, y = averageScGPTScore, alpha = factor(SameTAD), color = "Scaled TAD Gene Pairs"), size = 0.6) +
  geom_errorbar(data = dashed_data, aes(x = dist1_int, ymin = lower_ci, ymax = upper_ci, color = "Scaled TAD Gene Pairs"), size = 0.4) +
  
  ylim(c(0.04, 0.17)) + 
  labs(x = "Distance (KB)", y = "Contextual Similarity") +
  
  # Adjust color, linetype, and alpha to map correctly to the legend
  scale_color_manual(values = c("non-TAD Gene Pairs" = fourth_color, "Scaled TAD Gene Pairs" = second_color),
                     limits = c("Scaled TAD Gene Pairs", "non-TAD Gene Pairs")) +
  scale_linetype_manual(values = c("solid", "twodash")) +  
  scale_alpha_manual(values = c(1, 1)) +  # Adjust alpha if needed
  
  # Include linetype and color in the legend
  guides(
    color = guide_legend(override.aes = list(linetype = c("twodash","solid"))),
    linetype = "none",
    alpha = "none"
  ) +
  
  theme_classic() +
  theme(
    text = element_text(hjust = 0.5, size = 7, family = "Arial"),
    axis.text.x = element_text(size = 7, family = "Arial"),
    axis.text.y = element_text(size = 7, family = "Arial"),
    axis.title.x = element_text(size = 9, family = "Arial"),
    axis.title.y = element_text(size = 9, family = "Arial"),
    legend.title = element_blank(),  
    legend.text = element_text(size = 5.5, family = "Arial"),
    legend.position = c(0.8, 0.93),  
    legend.key.size = unit(0.4, "cm"),  
    legend.spacing.y = unit(0.01, "cm")
  )
boxes3 <- grid.arrange(b3, ncol = 1)

################################################################################
# Figure S2
################################################################################

samChrPairs1_Grouped <- samChrPairs1 %>% 
  group_by(SameTAD, dist1_int) %>% 
  summarise(
    averageScGPTScore = mean(coexp.x, na.rm = TRUE),
    n = n(),
    se = sd(scGPTCorr, na.rm = TRUE) / sqrt(n()),
    t_value = qt(0.975, df = n() - 1),  # Critical t-value for 95% CI
    lower_ci = averageScGPTScore - t_value * se,
    upper_ci = averageScGPTScore + t_value * se
  ) %>% 
  ungroup() 

pvals <- samChrPairs1_Grouped
samChrPairs1 %>% 
  group_by(dist1_int) %>%
  summarize(
    p_value = t.test(scGPTCorr ~ SameTAD)$p.value
  )

library(showtext)
font_add(family = "Arial", regular = "~/Arial.ttf")
palette_colors <- nord::nord("victory_bonds")

# Get the second color from the palette
second_color <- "#CC3311"
fourth_color <- "#009988"
samChrPairs1_Grouped$SameTAD <- factor(samChrPairs1_Grouped$SameTAD, levels = c(1,0))
b3 <- ggplot(samChrPairs1_Grouped, aes(x = dist1_int, y = averageScGPTScore, color = factor(SameTAD))) +
  geom_line(size=0.5) +
  geom_point(size=0.6) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                size = 0.4) + 
  labs(x = "Distance (KB)", y = "Coexpression", fill = "SameTAD") +
  scale_color_manual(values = c(second_color,fourth_color), labels = c("TAD Gene Pairs", "non-TAD Gene Pairs")) +
  theme_classic() +
  theme(
    text = element_text(hjust = 0.5, size = 7, family = "Arial"),
    axis.text.x = element_text(size = 7, family = "Arial"),
    axis.text.y = element_text(size = 7, family = "Arial"),
    axis.title.x = element_text(size = 9, family = "Arial"),
    axis.title.y = element_text(size = 9, family = "Arial"),
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 5.5, family = "Arial"),
    legend.position = c(0.83, 0.93), # Smaller legend text,
    legend.key.size = unit(0.4, "cm"),  
    legend.spacing.y = unit(0.01, "cm")
  ) 

boxes3 <- grid.arrange(b3, ncol = 1)

################################################################################
# Figure 4F
################################################################################

library(ggplot2)
library(gridExtra)
library(tidyverse)
library(dplyr)
library(gridExtra)

process_K562 <- function(fname) {
  d <- read.csv(fname, as.is = TRUE, nrows = 5)
  colC <- rep("numeric", ncol(d))
  
  for (i in 1:ncol(d)) {
    if (is.character(d[1, i])) {
      colC[i] <- "character"
    }
  }
  read.csv(fname, colClasses = colC, sep = ',', row.names = NULL) %>% 
    dplyr::select(gene, tadFCmean, nontadFCmean,overall) 
}

# List of file paths
file_paths <- paste0("data/for_plotting/K562_clustered21-", 1:5, ".csv")

# Process all K562 datasets
K562_list <- lapply(file_paths, process_K562)

# Merge all datasets by "gene"
# names(K562_list) <- c("_a", "_b", "_c", "_d","_e","_f")  # Replace with your actual suffixes
names(K562_list) <- c("_a", "_b", "_c", "_d", "_e")
# Function to rename columns by adding suffixes
rename_cols <- function(df, suffix) {
  colnames(df) <- ifelse(
    colnames(df) == "gene", 
    "gene", 
    paste0(colnames(df), suffix)
  )
  return(df)
}

# Apply the renaming function to each DataFrame in the list
K562_list <- Map(rename_cols, K562_list, names(K562_list))

# Merge all DataFrames on the 'gene' column
merged_K562 <- Reduce(function(x, y) merge(x, y, by = "gene"), K562_list)

final_K562 <- merged_K562 %>%
  rowwise() %>%
  mutate(tadFC = mean(c_across(starts_with("tadFC"))),
         nontadFC = mean(c_across(starts_with("nontadFC"))),
         overall = mean(c_across(starts_with("overall")))
  ) %>%
  ungroup() %>%
  dplyr::select(gene, tadFC, nontadFC, overall) 

# 
only_numeric <- final_K562 %>% dplyr::select(tadFC, nontadFC)

# Perform PCA
pca_result <- prcomp(only_numeric, scale. = TRUE, rank. = 2)  # rank. = 2 ensures 2 PCs
summary(pca_result)

k562_pcs <- cbind(final_K562[1], pca_result$x)


k562pcsTop <- k562_pcs %>%
  arrange(desc(PC2)) %>%
  slice(1:50) %>%
  dplyr::select(gene)

k562pcsBottom <- k562_pcs %>%
  arrange(PC2) %>%
  slice(1:50) %>%
  dplyr::select(gene)

k562pcsTop1 <- k562_pcs %>%
  arrange(desc(PC1)) %>%
  slice(1:50) %>%
  dplyr::select(gene)

d <- final_K562 %>%
  filter(gene %in% k562pcsTop1$gene)

d1 <- final_K562 %>% 
  filter(gene %in% k562pcsTop$gene)

d2 <- final_K562 %>% 
  filter(gene %in% k562pcsBottom$gene)

font_add(family = "Arial", regular = "~/Arial.ttf")
palette_colors <- nord::nord("victory_bonds")
second_color <- "#CC3311"
fourth_color <- "#009988"


b11 <- ggplot() +
  # Combine data for legend creation
  geom_point(data = final_K562, aes(x = tadFC - nontadFC, y = overall), alpha = 1, size = 0.02, color = "grey") + 
  geom_point(data = d1, aes(x = tadFC - nontadFC, y = overall, color = "Strong non-TAD Gene \nDisruption"), alpha = 1, size = 0.02) +
  geom_point(data = d2, aes(x = tadFC - nontadFC, y = overall, color = "Strong TAD Gene \nDisruption"), alpha = 1, size = 0.02) +
  labs(x = "Difference in log2FC Variance \nin TAD Genes vs non-TAD Genes", y = "Overall Variance in log2FC") +
  scale_color_manual(values = c("Strong TAD Gene \nDisruption" = second_color, "Strong non-TAD Gene \nDisruption" = fourth_color)) +
  theme_classic() +
  theme(
    text = element_text(hjust = 0.5, size = 7, family = "Arial"),
    axis.text.x = element_text(size = 7, family = "Arial"),
    axis.text.y = element_text(size = 7, family = "Arial"),
    axis.title.x = element_text(size = 9, family = "Arial"),
    axis.title.y = element_text(size = 9, family = "Arial"),
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 5.5, family = "Arial"),
    legend.position = c(0.82,0.9),   # Position of legend
    legend.key.size = unit(0.4, "cm")
  )


boxes11 <- grid.arrange(b11, ncol = 1)

################################################################################
# Figure 4G
################################################################################

top <- k562pcsTop$gene
bottom <- k562pcsBottom$gene

library(enrichR)
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}

if (websiteLive) t <- listEnrichrDbs() 

dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023")
if (websiteLive) {
  enriched <- enrichr(bottom, dbs)
}

enriched[[3]] %>% arrange(P.value) %>% 
  select(Term, Overlap, Genes, Adjusted.P.value) %>% 
  View()

f_varBottomTAD <- enriched[[3]] %>% 
  filter(Term %in% c('Positive Regulation Of Protein Localization To Nucleus (GO:1900182)',
                     "Regulation Of Transcription Elongation By RNA Polymerase II (GO:0034243)",
                     "RNA 3'-End Processing (GO:0031123)",
                     "snRNA Processing (GO:0016180)",
                     'Internal Peptidyl-Lysine Acetylation (GO:0018393)',
                     'Regulation Of DNA-templated Transcription Elongation (GO:0032784)',
                     'Histone Acetylation (GO:0016573)'
                     ))

f_varBottomTAD$Term <- gsub("^(.*?)\\s*\\(.*$", "\\1", f_varBottomTAD$Term)

b14 <- plotEnrich(f_varBottomTAD, showTerms = 8, numChar = 100, y = "Count", orderBy = "P.value", xlab = 'Enriched Terms', ylab = "Gene Count",
                  title = 'TAD Enrichr Analysis') +
  theme_classic() +
  labs(fill = "p value") +
  theme(
    text = element_text(hjust = 0.5, size = 7, family = "Arial"),
    axis.text.x = element_text(size = 7, family = "Arial"),
    axis.text.y = element_text(size = 6, family = "Arial", lineheight = 0.8),  # Adjust lineheight for y-axis text
    axis.title.x = element_text(size = 9, family = "Arial"),
    axis.title.y = element_text(size = 9, family = "Arial"),  # Remove legend title
    legend.text = element_text(size = 3, family = "Arial"),
    legend.title = element_text(size = 3, family = "Arial"),
    legend.position = c(0.92,0.28),   # Position of legend
    legend.key.size = unit(0.3, "cm")
  )

boxes14 <- grid.arrange(b14, ncol = 1)

################################################################################
# Figure S3
################################################################################

library(enrichR)
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}

if (websiteLive) t <- listEnrichrDbs() 

dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023")
if (websiteLive) {
  enriched <- enrichr(top, dbs)
}

enriched[[3]] %>% arrange(P.value) %>% 
  select(Term, Overlap, Genes, Adjusted.P.value) %>% 
  View()


f_varTopTAD <- enriched[[3]] %>% 
  filter(Term %in% c("Ribosome Biogenesis (GO:0042254)",
                     "Maturation Of SSU-rRNA (GO:0030490)",
                     "Cytoplasmic Translation (GO:0002181)",
                     "rRNA Processing (GO:0006364)",
                     "Translation (GO:0006412)",
                     "Peptidyl-Arginine Modification (GO:0018195)",
                     "ncRNA Processing (GO:0034470)"))

f_varTopTAD$Term <- gsub("^(.*?)\\s*\\(.*$", "\\1", f_varTopTAD$Term)


b14 <- plotEnrich(f_varTopTAD, showTerms = 8, numChar = 100, y = "Count", orderBy = "P.value", xlab = 'Enriched Terms', ylab = "Gene Count",
                  title = 'TAD Enrichr Analysis') +
  theme_classic() +
  labs(fill = "p value") +
  theme(
    text = element_text(hjust = 0.5, size = 7, family = "Arial"),
    axis.text.x = element_text(size = 7, family = "Arial"),
    axis.text.y = element_text(size = 6, family = "Arial", lineheight = 0.8),  # Adjust lineheight for y-axis text
    axis.title.x = element_text(size = 9, family = "Arial"),
    axis.title.y = element_text(size = 9, family = "Arial"),  # Remove legend title
    legend.text = element_text(size = 3, family = "Arial"),
    legend.title = element_text(size = 3, family = "Arial"),
    legend.position = c(0.92,0.28),   # Position of legend
    legend.key.size = unit(0.3, "cm")
  )

boxes14 <- grid.arrange(b14, ncol = 1)

################################################################################
# Figure 4H
################################################################################

fname <- "../data/pair_regcoexp1.csv"
d <- read.csv(fname, as.is=T, nrows=5)
colC <- rep("numeric",ncol(d))
for (i in 1:ncol(d)) {
  if (is.character(d[1,i])) {
    colC[i]="character"
  }
}

# Load and preprocess gene pairs data
gene_pairs <- read.csv(fname, colClasses = colC, sep = ',') %>%
  rename(distance_range = pair_dist_range) %>%
  filter(same_chr == 1) %>%
  mutate(
    I20a = ifelse(pair_dist < 20000, 1, 0),
    I1500 = ifelse(pair_dist > 1e6 & pair_dist <= 1.5e6, 1, 0)
  )
gene_pairs <- gene_pairs %>%
  mutate(corrZ = fisherz(corr))

gene_pairs_plot <- gene_pairs %>%
  filter(pair_dist <= 1e6)

mean_gene_pairs_plot <- gene_pairs %>%
  group_by(tad, adj_genes) %>%
  summarise(mean_corr = fisherz2r(mean(corrZ, na.rm = TRUE)),
            t_value = qt(0.975, df = n() - 1), 
            se = sd(corr, na.rm = TRUE) / sqrt(n()), # Critical t-value for 95% CI
            lower_ci = fisherz2r(mean_corr - t_value * se),
            upper_ci = fisherz2r(mean_corr + t_value * se),
            n = n()) %>% 
  ungroup() %>% 
  mutate(adj_genes1 = recode(adj_genes,
                             `0` = "Non-Adjacent", 
                             `1` = "Adjacent"))

mean_gene_pairs_sameStrand <- gene_pairs %>%
  filter(adj_genes == 1) %>% 
  group_by(tad, same_strand_adj) %>%
  summarise(mean_corr = fisherz2r(mean(corrZ, na.rm = TRUE)),
            t_value = qt(0.975, df = n() - 1), 
            se = sd(corr, na.rm = TRUE) / sqrt(n()), # Critical t-value for 95% CI
            lower_ci = fisherz2r(mean_corr - t_value * se),
            upper_ci = fisherz2r(mean_corr + t_value * se),
            n = n()) %>% 
  ungroup() %>% 
  mutate(same_strand_adj1 = recode(same_strand_adj,
                                   `0` = "Adjacent", 
                                   `1` = "Same-Strand Adjacent"),
         tad1 = recode(tad,
                       `0` = "non-TAD", 
                       `1` = "TAD"))

font_add(family = "Arial", regular = "~/Arial.ttf")
palette_colors <- nord::nord("victory_bonds")

# Get the second color from the palette
second_color <- "#CC3311"
fourth_color <- "#009988"

mean_gene_pairs_sameStrand$same_strand_adj1 <- factor(mean_gene_pairs_sameStrand$same_strand_adj1, levels = c("Same-Strand Adjacent", "Adjacent"))
mean_gene_pairs_sameStrand$tad1 <- factor(mean_gene_pairs_sameStrand$tad1, levels = c("TAD","non-TAD"))

adj_genes_pairs <- gene_pairs %>% 
  filter(adj_genes==1)

same_strand_adj_pairs <- adj_genes_pairs %>% 
  filter(tad == 1, pair_dist < 2e6)
non_same_strand_adj_pairs <- adj_genes_pairs %>% 
  filter(tad == 0, pair_dist < 2e6)

t_tad <- t.test(same_strand_adj_pairs[same_strand_adj_pairs$same_strand_adj == 1,]$corrZ, 
                same_strand_adj_pairs[same_strand_adj_pairs$same_strand_adj == 0,]$corrZ)
t_nontad <- t.test(non_same_strand_adj_pairs[non_same_strand_adj_pairs$same_strand_adj == 1,]$corrZ, 
                   non_same_strand_adj_pairs[non_same_strand_adj_pairs$same_strand_adj == 0,]$corrZ)

mean_gene_pairs_sameStrand1 <- mean_gene_pairs_sameStrand %>% 
  mutate(strand_tad_combo = paste(same_strand_adj1, tad1))

mean_gene_pairs_sameStrand1$strand_tad_combo <- factor(mean_gene_pairs_sameStrand1$strand_tad_combo, 
                                                       levels = c("Same-Strand Adjacent TAD","Adjacent TAD",
                                                                  "Same-Strand Adjacent non-TAD","Adjacent non-TAD"))

b5 <- ggplot(mean_gene_pairs_sameStrand1, aes(x = factor(tad1), y = mean_corr, fill = strand_tad_combo)) +
  geom_bar(stat = "identity", alpha=0.95, position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2, position = position_dodge(width = 0.9)) +
  labs(x = NULL, y = "Contextual Similarity", fill = NULL) +
  theme_classic() +
  ylim(c(0,0.17)) +
  scale_fill_manual(values = c("#b93518","#e77649","#471184","#6850a1"), labels = c("Same direction (shared TAD)", "Opposite direction (shared TAD)",
                                                                                    "Same direction (different TADs)", "Opposite direction (different TADs)")) +
  theme(
    text = element_text(hjust = 0.5, size = 7, face = "bold", family = "Arial"),
    axis.text.x = element_text(size = 7, family = "Arial"),
    axis.text.y = element_text(size = 7, family = "Arial"),
    axis.title.x = element_text(size = 9, face = "bold", family = "Arial"),
    axis.title.y = element_text(size = 9, face = "bold", family = "Arial"),
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 5.5, family = "Arial"),
    legend.position = c(0.7, 0.93), # Smaller legend text,
    legend.key.size = unit(0.3, "cm"),  
    legend.spacing.y = unit(0.01, "cm")
  ) +
  annotate(
    "text",
    x = "TAD",     # The variable name for x position
    y = max(mean_gene_pairs_plot$mean_corr) + 0.028,  # Adjust y position (slightly above the tallest bar)
    label = paste0("italic('t') == ", round(t_tad$statistic, 2)),
    size = 2,
    family = "Arial",
    fontface = "bold",
    hjust = 0.5,
    parse = TRUE
  ) +
  annotate(
    "text",
    x = "non-TAD",     # The variable name for x position
    y = max(mean_gene_pairs_plot[mean_gene_pairs_plot$tad==0,]$mean_corr) + 0.025,  # Adjust y position (slightly above the tallest bar)
    label = paste0("italic('t') == ", round(t_nontad$statistic, 2)),
    size = 2,
    family = "Arial",
    fontface = "bold",
    hjust = 0.5,
    parse = TRUE
  ) +
  annotate(
    "text",
    x = "TAD",     # The variable name for x position
    y = max(mean_gene_pairs_plot$mean_corr) + 0.018,  # Adjust y position (slightly above the tallest bar)
    label = paste0("italic('p') == ", formatC(t_tad$p.value, format = "e", digits = 2)),
    size = 2,
    family = "Arial",
    fontface = "bold",
    hjust = 0.5,
    parse = TRUE
  ) +
  annotate(
    "text",
    x = "non-TAD",     # The variable name for x position
    y = max(mean_gene_pairs_plot[mean_gene_pairs_plot$tad==0,]$mean_corr) + 0.015,  # Adjust y position (slightly above the tallest bar)
    label = paste0("italic('p') == ", round(t_nontad$p.value, 3)),
    size = 2,
    family = "Arial",
    fontface = "bold",
    hjust = 0.5,
    parse = TRUE
  )

boxes5 <- grid.arrange(b5, ncol = 1)