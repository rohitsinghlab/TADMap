# install the following R libraries with "install.packages("").
rm(list = ls())
library(tidyverse)      # Used for data manipulation and dplyr operations
library(psych)          # Used for fisherz() function
library(gridExtra)      # Used for grid.arrange()
library(reshape2)       # Used for melt()
library(nord)           # Used for nord() color palette
library(showtext)       # Used for font manipulation
library(MASS)           # Used for kde2d()
library(fields)         # Used for interp.surface()
library(cocor)          # Used for cocor() function

# replace with your data location
fnameA <- "../data/pair_regcoexp1.csv"
d <- read.csv(fnameA, as.is = T, nrows = 5)
colC <- rep("numeric", ncol(d))
for (i in 1:ncol(d)) {
  if (is.character(d[1, i])) {
    colC[i] = "character"
  }
}
  
gene_pairsexp <- read.csv(fnameA, colClasses = colC, sep = ',') %>%
  mutate(
    coexp1 = fisherz(coexp),
    corr1 = fisherz(corr),
    corcoexpRatio = corr / coexp,
  ) %>% 
  mutate(normed_coexp = (coexp - median(coexp,na.rm=TRUE)),
         normed_corr = (corr - median(corr, na.rm=TRUE)),
         normed_coexp1 =(coexp1 - median(coexp1, na.rm=TRUE)),
         normed_corr1 = (corr1 - median(corr1, na.rm=TRUE)),
         std_normed_coexp = (coexp - median(coexp,na.rm=TRUE))/IQR(coexp, na.rm=TRUE),
         std_normed_corr = (corr - median(corr, na.rm=TRUE))/IQR(corr, na.rm=TRUE),
         std_normed_coexp1 =(coexp1 - median(coexp1, na.rm=TRUE))/IQR(coexp1, na.rm=TRUE),
         std_normed_corr1 = (corr1 - median(corr1, na.rm=TRUE))/IQR(corr1, na.rm=TRUE))


gene_pairsexp_long <- melt(gene_pairsexp, measure.vars = c("normed_coexp", "normed_corr"), 
                           variable.name = "Type", value.name = "Correlation")

################################################################################
# Figure 3A and Figure 3B
################################################################################

palette_colors <- nord::nord("victory_bonds")

second_color <- "#0077BB"
fourth_color <- "#EE7733"

var.test(gene_pairsexp_long$normed_coexp1, gene_pairsexp_long$normed_corr1, alternative="less")
var(gene_pairsexp_long$coexp1)
var(gene_pairsexp_long$corr1)

plot1 <- ggplot(filter(gene_pairsexp_long, Type == "normed_coexp"), aes(x = Correlation, fill = Type, color = Type)) +
  geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.7, position = "identity") +
  xlab("Coexpression") +
  ylab("Frequency (Normalized)") +  
  scale_color_manual(values=c(second_color, second_color)) +
  scale_fill_manual(values = c(second_color)) +
  theme_classic() +
  theme(text = element_text(family = "Arial"),
        legend.position = "none",  # Remove legends
        ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 7, face = "bold", family = "Arial"),
    axis.text.x = element_text(size = 7, family = "Arial"),
    axis.text.y = element_text(size = 7, family = "Arial"),
    axis.title.x = element_text(size = 8, face = "bold", family = "Arial"),
    axis.title.y = element_text(size = 8, face = "bold", family = "Arial")
  )# Ensure y-axis title is blank

plot1arr <- grid.arrange(plot1, ncol=1)
ggsave("~/projects/proj1/data/paper_images/Fig3/coexp_spread_hist.png",
       plot = plot1arr, width = 3, height = 2, dpi = 500)
ggsave("~/projects/proj1/data/paper_images/Fig3/coexp_spread_hist.pdf",
       plot = plot1arr, width = 3, height = 2, dpi = 500, device = cairo_pdf)


# Create the second plot
plot2 <- ggplot(filter(gene_pairsexp_long, Type == "normed_corr"), aes(x = Correlation, fill = Type, color = Type)) +
  geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.7, position = "identity") +
  xlab("Contextual Similarity") +
  ylab("Frequency (Normalized)") +  # Remove y-axis label
  scale_color_manual(values=c(fourth_color, fourth_color)) +
  scale_fill_manual(values = c(fourth_color)) +
  theme_classic() +
  theme(text = element_text(family = "Arial"),
        legend.position = "none") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 7, face = "bold", family = "Arial"),
    axis.text.x = element_text(size = 7, family = "Arial"),
    axis.text.y = element_text(size = 7, family = "Arial"),
    axis.title.x = element_text(size = 8, face = "bold", family = "Arial"),
    axis.title.y = element_text(size = 8, face = "bold", family = "Arial")
  )# Ensure y-axis title is blank

plot2arr <- grid.arrange(plot2, ncol=1)
ggsave("~/projects/proj1/data/paper_images/Fig3/contexSim_spread_hist.png",
       plot = plot2arr, width = 3, height = 2, dpi = 500)
ggsave("~/projects/proj1/data/paper_images/Fig3/contexSim_spread_hist.pdf",
       plot = plot2arr, width = 3, height = 2, dpi = 500, device = cairo_pdf)

################################################################################
# Figure 3C
################################################################################

# replace with your data location
fnameA <- "../data/samChrPairs.csv"
library(cocor)
d <- read.csv(fnameA, as.is = T, nrows = 5)
colC <- rep("numeric", ncol(d))
for (i in 1:ncol(d)) {
  if (is.character(d[1, i])) {
    colC[i] = "character"
  }
}
# add fischer z information, only TAD pairs
samChrPairs <- read.csv(fnameA, colClasses = colC, sep = ',') %>% 
  mutate(coexpT = fisherz(coexp))

samChrPairs1 <- samChrPairs %>% 
  filter(dist < 1e6)

cor_with_ci <- function(x, y) {
  cor_value <- cor(x, y, method = 'pearson')
  n <- length(x)
  stderr <- sqrt((1 - cor_value^2) / (n - 2))
  ci_upper <- cor_value + 1.96 * stderr
  ci_lower <- cor_value - 1.96 * stderr
  return(c(cor = cor_value, lower = ci_lower, upper = ci_upper))
}
df <- samChrPairs1

cocor(~scGPTCorrT + goCorr1kT | scGPTCorrT + gene2vecT, data = df, alternative = "less")
cocor(~coexpT + goCorr1kT | coexpT + gene2vecT, data = df, alternative = "less")

cocor(~gene2vecT + scGPTCorrT | gene2vecT + coexpT, data = df, alternative = "greater")
cocor(~goCorr1kT + scGPTCorrT | goCorr1kT + coexpT, data = df, alternative = "greater")


# Create a data frame with correlations and CIs
results <- data.frame(
  group = c("Contextual Similarity vs HiG2Vec", "Contextual Similarity vs Gene2Vec", "Coexpression vs HiG2Vec", "Coexpression vs Gene2Vec"),
  variable_x = c("Contextual Similarity", "Contextual Similarity", "Coexpression", "Coexpression"),
  variable_color = c("HiG2Vec", "Gene2Vec", "HiG2Vec", "Gene2Vec"),
  cor = c(cor_with_ci(df$scGPTCorrT, df$goCorr1kT)["cor"],
          cor_with_ci(df$scGPTCorrT, df$gene2vecT)["cor"],
          cor_with_ci(df$coexp, df$goCorr1kT)["cor"],
          cor_with_ci(df$coexp, df$gene2vecT)["cor"]),
  lower_ci = c(cor_with_ci(df$scGPTCorrT, df$goCorr1kT)["lower"],
               cor_with_ci(df$scGPTCorrT, df$gene2vecT)["lower"],
               cor_with_ci(df$coexp, df$goCorr1kT)["lower"],
               cor_with_ci(df$coexp, df$gene2vecT)["lower"]),
  upper_ci = c(cor_with_ci(df$scGPTCorrT, df$goCorr1kT)["upper"],
               cor_with_ci(df$scGPTCorrT, df$gene2vecT)["upper"],
               cor_with_ci(df$coexp, df$goCorr1kT)["upper"],
               cor_with_ci(df$coexp, df$gene2vecT)["upper"])
)
library(showtext)
font_add(family = "Arial", regular = "~/Arial.ttf")
palette_colors <- nord::nord("victory_bonds")

# Get the second color from the palette
second_color <- "#0077BB"
fourth_color <- "#EE7733"

# Your existing code above...

# Adjust the order of variable_x to switch "Coexpression" and "scGPT"
results$variable_x <- factor(results$variable_x, levels = c("Contextual Similarity", "Coexpression"))

# Proceed with your ggplot code
b <- ggplot(results, aes(x = variable_color, y = cor, fill = variable_x)) +
  geom_bar(stat = "identity", alpha=0.95, position = position_dodge(width = 0.9)) +
  labs(x = NULL, y = "Correlation", fill = NULL) +
  theme_classic() +
  ylim(c(0,0.3))+
  scale_fill_manual(values = c(fourth_color, second_color), labels = c("Contextual Similarity","Coexpression")) +
  theme(
    text = element_text(hjust = 0.5, size = 7, face = "bold", family = "Arial"),
    axis.text.x = element_text(size = 7, family = "Arial"),
    axis.text.y = element_text(size = 7, family = "Arial"),
    axis.title.x = element_text(size = 9, face = "bold", family = "Arial"),
    axis.title.y = element_text(size = 9, face = "bold", family = "Arial"),
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 5.5, family = "Arial"),
    legend.position = c(0.8, 0.96), # Smaller legend text,
    legend.key.size = unit(0.3, "cm"),  
    legend.spacing.y = unit(0.01, "cm")
  )

boxes1 <- grid.arrange(b, ncol = 1)
ggsave("~/projects/proj1/data/paper_images/Fig3/embCompAllPairs.png",
       plot = boxes1, width = 3, height = 2, dpi = 500)
ggsave("~/projects/proj1/data/paper_images/Fig3/embCompAllPairs.pdf",
       plot = boxes1, width = 3, height = 2, dpi = 500, device = cairo_pdf)



################################################################################
# Figure 3D
################################################################################

# Load necessary libraries
library(ggplot2)
library(MASS)   # For kde2d()
library(fields) # For interp.surface()
d1 <- gene_pairsexp 
# Calculate density using kde2d()
density_data <- with(d1, MASS::kde2d(normed_coexp, normed_corr, n = 125))

# Interpolate density for each point
d1$density <- with(d1, interp.surface(density_data, cbind(normed_coexp, normed_corr)))

# Assuming you have calculated the Spearman correlation and saved it as 'spearman_corr'
spearman_corr <- round(cor(d1$normed_coexp, d1$normed_corr, method = "spearman"), 3)
pearson_corr <- round(cor(d1$normed_coexp, d1$normed_corr, method = "pearson"), 3)# Example value, replace with your calculation
spearman_corr
pearson_corr
p3 <- ggplot(d1, aes(x = normed_coexp, y = normed_corr)) +
  geom_point(aes(color = density), size = 0.05) +
  scale_color_viridis_c(trans = "log",
                        breaks = c(min(d1$density), max(d1$density)),  # Show only min and max
                        labels = c(0, 20)) +  # Adjust color scale
  xlab("Coexpression") +
  ylab("Contextual Similarity") +
  xlim(-0.6, 0.6) +
  labs(color = "Density") +  # Relabel the color legend
  guides(color = guide_colorbar(barwidth = 0.3, barheight = 7,
         frame.colour = "black", frame.linewidth = 0.3,)) +
  annotate("text", x = 0.07, y = 0.89, fontface = "bold", label = paste("Spearman Correlation:", spearman_corr), 
           size = 1.6, hjust = 0, family = "Arial") + 
  annotate("text", x = 0.07, y = 0.94, fontface = "bold", label = paste("Pearson Correlation:", pearson_corr), 
           size = 1.6, hjust = 0, family = "Arial") + # Add correlation text
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 7, face = "bold", family = "Arial"),
    axis.text.x = element_text(size = 7, family = "Arial"),
    axis.text.y = element_text(size = 7, family = "Arial"),
    axis.title.x = element_text(size = 8, face = "bold", family = "Arial"),
    axis.title.y = element_text(size = 8, face = "bold", family = "Arial"),
    legend.title = element_text(size = 5, family = "Arial"),  # Adjust legend title
    legend.text = element_text(size = 6, family = "Arial"),
    legend.position = "right")

# Print the plot
plot3arr <- grid.arrange(p3, ncol=1)
ggsave("~/projects/proj1/data/paper_images/Fig3/contexSim_coexp_dens.png",
       plot = plot3arr, width = 3, height = 2, dpi = 500)
ggsave("~/projects/proj1/data/paper_images/Fig3/contexSim_coexp_dens.pdf",
       plot = plot3arr, width = 3, height = 2, dpi = 500, device = cairo_pdf)

################################################################################
# Figure 3F
################################################################################

# Get the second color from the palette
second_color <- "#0077BB"
fourth_color <- "#EE7733"

# PCDH

d1 <- gene_pairsexp 
d <- gene_pairsexp %>% 
  filter(grepl('^PCDH', g1) & grepl('^PCDH', g2)) 

mean(d$normed_corr)
mean(d$normed_coexp)

median(d$normed_corr)
median(d$normed_coexp)


pcdh_spearman_corr <- round(cor(d$normed_coexp, d$normed_corr, method = "spearman"), 3)
p4 <- ggplot() +
  geom_point(data = d1, aes(x = normed_coexp, y = normed_corr), color = "grey", size = 0.05) +
  geom_point(data = d, aes(x = normed_coexp, y = normed_corr), color = "black", size = 0.05) +
  annotate("text", x = 0.17, y = 0.94, fontface = "bold", label = paste("Spearman Correlation:", pcdh_spearman_corr), 
           size = 1.6, hjust = 0, family = "Arial") +
  xlab("Coexpression") +
  ylab("Contextual Similarity") +
  ggtitle("Protocadherin Gene Pairs") +
  xlim(-0.6, 0.6) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9, face = "bold", family = "Arial"),
    axis.text.x = element_text(size = 7, family = "Arial"),
    axis.text.y = element_text(size = 7, family = "Arial"),
    axis.title.x = element_text(size = 8, face = "bold", family = "Arial"),
    axis.title.y = element_text(size = 8, face = "bold", family = "Arial")
  )
plot4arr <- grid.arrange(p4, ncol=1)
ggsave("~/projects/proj1/data/paper_images/Fig3/contexSim_coexp_pcdh.png",
       plot = plot4arr, width = 3, height = 2, dpi = 500)
ggsave("~/projects/proj1/data/paper_images/Fig3/contexSim_coexp_pcdh.pdf",
       plot = plot4arr, width = 3, height = 2, dpi = 500, device = cairo_pdf)

################################################################################
# Figure 3E
################################################################################

# OR
d <- gene_pairsexp %>%
  filter(grepl('^OR\\d', g1) & grepl('^OR\\d', g2)) 

mean(d$normed_corr)
mean(d$normed_coexp)

median(d$normed_corr)
median(d$normed_coexp)

or_spearman_corr <- round(cor(d$normed_coexp, d$normed_corr, method = "spearman"), 3)
p5 <- ggplot() +
  geom_point(data = d1, aes(x = normed_coexp, y = normed_corr), color = "grey", size = 0.05) +
  geom_point(data = d, aes(x = normed_coexp, y = normed_corr), color = "black", size = 0.05) +
  annotate("text", x = 0.17, y = 0.94, fontface = "bold", label = paste("Spearman Correlation:", or_spearman_corr),
           size = 1.6, hjust = 0, family = "Arial") +
  xlab("Coexpression") +
  ylab("Contextual Similarity") +
  ggtitle("Olfactory Receptor Gene Pairs") +
  xlim(-0.6, 0.6) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9, face = "bold", family = "Arial"),
    axis.text.x = element_text(size = 7, family = "Arial"),
    axis.text.y = element_text(size = 7, family = "Arial"),
    axis.title.x = element_text(size = 8, face = "bold", family = "Arial"),
    axis.title.y = element_text(size = 8, face = "bold", family = "Arial")
  )
plot5arr <- grid.arrange(p5, ncol=1)
ggsave("~/projects/proj1/data/paper_images/Fig3/contexSim_coexp_or.png",
       plot = plot5arr, width = 3, height = 2, dpi = 500)
ggsave("~/projects/proj1/data/paper_images/Fig3/contexSim_coexp_or.pdf",
       plot = plot5arr, width = 3, height = 2, dpi = 500, device = cairo_pdf)

################################################################################
# Figure 3G
################################################################################

# Histones
d <- gene_pairsexp %>%
  filter(grepl('^H[1-5]', g1) & grepl('^H[1-5]', g2))

mean(d$normed_corr)
mean(d$normed_coexp)

median(d$normed_corr)
median(d$normed_coexp)

histone_spearman_corr <- round(cor(d$normed_coexp, d$normed_corr, method = "spearman"), 3)

p6 <- ggplot() +
  geom_point(data = d1, aes(x = normed_coexp, y = normed_corr), color = "grey", size = 0.05) +
  geom_point(data = d, aes(x = normed_coexp, y = normed_corr), color = "black", size = 0.05) +
  xlab("Coexpression") +
  ylab("Contextual Similarity") +
  ggtitle("Histone Gene Pairs") +
  annotate("text", x = 0.17, y = 0.94, fontface = "bold", label = paste("Spearman Correlation:", histone_spearman_corr),
           size = 1.6, hjust = 0, family = "Arial") +
  xlim(-0.6, 0.6) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9, face = "bold", family = "Arial"),
    axis.text.x = element_text(size = 7, family = "Arial"),
    axis.text.y = element_text(size = 7, family = "Arial"),
    axis.title.x = element_text(size = 8, face = "bold", family = "Arial"),
    axis.title.y = element_text(size = 8, face = "bold", family = "Arial")
  )

plot6arr <- grid.arrange(p6, ncol=1)
ggsave("~/projects/proj1/data/paper_images/Fig3/contexSim_coexp_hist.png",
       plot = plot6arr, width = 3, height = 2, dpi = 500)
ggsave("~/projects/proj1/data/paper_images/Fig3/contexSim_coexp_hist.pdf",
       plot = plot6arr, width = 3, height = 2, dpi = 500, device = cairo_pdf)

################################################################################
# Figure S1
################################################################################

# KRT
d <- gene_pairsexp %>%
  filter(grepl('^KRT', g1) & grepl('^KRT', g2))
krt_spearman_corr <- round(cor(d$normed_coexp, d$normed_corr, method = "spearman"), 3)

p7 <- ggplot() +
  geom_point(data = d1, aes(x = normed_coexp, y = normed_corr), color = "grey", size = 0.05) +
  geom_point(data = d, aes(x = normed_coexp, y = normed_corr), color = "black", size = 0.05) +
  xlab("Coexpression") +
  ylab("Contextual Similarity") +
  ggtitle("Keratin Gene Pairs") +
  annotate("text", x = 0.17, y = 0.94, fontface = "bold", label = paste("Spearman Correlation:", krt_spearman_corr),
           size = 1.6, hjust = 0, family = "Arial") +
  xlim(-0.6, 0.6) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9, face = "bold", family = "Arial"),
    axis.text.x = element_text(size = 7, family = "Arial"),
    axis.text.y = element_text(size = 7, family = "Arial"),
    axis.title.x = element_text(size = 8, face = "bold", family = "Arial"),
    axis.title.y = element_text(size = 8, face = "bold", family = "Arial")
  )

plot7arr <- grid.arrange(p7, ncol=1)
ggsave("~/projects/proj1/data/paper_images/Fig3/contexSim_coexp_krt.png",
       plot = plot7arr, width = 3, height = 2, dpi = 500)
ggsave("~/projects/proj1/data/paper_images/Fig3/contexSim_coexp_krt.pdf",
       plot = plot7arr, width = 3, height = 2, dpi = 500, device = cairo_pdf)


################################################################################
# Figure 3H
################################################################################

full <- gene_pairsexp %>%
  filter((grepl('^H[1-5]', g1) & grepl('^H[1-5]', g2)) |
           (grepl('^OR\\d', g1) & grepl('^OR\\d', g2)) |
           (grepl('^PCDH', g1) & grepl('^PCDH', g2)) |
           (grepl('^KRT', g1) & grepl('^KRT', g2))
  ) %>% 
  mutate(family = case_when(
    grepl('^H[1-5]', g1) & grepl('^H[1-5]', g2) ~ "Histone",
    grepl('^OR\\d', g1) & grepl('^OR\\d', g2) ~ "OR",
    grepl('^PCDH', g1) & grepl('^PCDH', g2) ~ "PCDH",
    grepl('^KRT', g1) & grepl('^KRT', g2) ~ "KRT")) %>% dplyr::select(g1, g2, std_normed_coexp, std_normed_corr, 
                                                                      std_normed_coexp1, std_normed_corr1, family)
# Vertical Orientation
full <- full %>%
  mutate(family = factor(family, levels = c("KRT", "Histone", "PCDH","OR")))

# Reshape the data into long format and rename the metrics
full_long <- full %>%
  pivot_longer(cols = c(std_normed_corr, std_normed_coexp), 
               names_to = "metric", 
               values_to = "value") %>%
  mutate(metric = recode(metric,
                         "std_normed_corr" = "Contextual Similarity",
                         "std_normed_coexp" = "Coexpression"))

# Find the range of your data without outliers to set appropriate limits
min_val <- min(full_long$value, na.rm = TRUE)
max_val <- quantile(full_long$value, 0.9, na.rm = TRUE)  # 99th percentile to avoid extreme outliers

t_test_results <- full_long %>%
  group_by(family) %>%
  summarise(t_test = list(t.test(value ~ metric)),
            t_statistic = t_test[[1]]$statistic,
            p_value = t_test[[1]]$p.value)


box <- ggplot(full_long, aes(x = family, y = value, fill = metric)) +
  geom_boxplot(outlier.shape = NA, size = 0.3, alpha = 0.8) +  # Smaller line width for box borders
  scale_fill_manual(values = c("Contextual Similarity" = fourth_color, "Coexpression" = second_color)) +
  scale_color_manual(values = c("Contextual Similarity" = fourth_color, "Coexpression" = second_color)) +
  labs(x = "Gene Family", y = "Standardized Values", fill = NULL) +  # Adjust labels and remove legend title
  ylim(min_val, 4.1) +
  coord_flip() + # Set y-axis limits based on data
  theme_classic() +
  theme(
    text = element_text(hjust = 0.5, size = 7, face = "bold", family = "Arial"),
    axis.text.x = element_text(size = 7, family = "Arial"),
    axis.text.y = element_text(size = 7, family = "Arial"),
    axis.title.x = element_text(size = 9, face = "bold", family = "Arial"),
    axis.title.y = element_text(size = 9, face = "bold", family = "Arial"),
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 3, family = "Arial"),
    legend.position = c(-0.3,0.98), # Smaller legend text,
    legend.key.size = unit(0.4, "cm")
  )


for (i in 1:nrow(t_test_results)) {
  family_name <- t_test_results$family[i]
  t_stat <- round(t_test_results$t_statistic[i], 1)
  p_val <- t_test_results$p_value[i]
  
  # Conditional p-value formatting
  if (p_val < 2.2e-16) {
    p_label <- "<2.2e-16"
  } else {
    p_label <- formatC(p_val, format = "e", digits = 2)
  }
  
  # Enclose the p-value in quotes for the plotmath expression
  p_label_str <- paste0("'", p_label, "'")
  
  # Annotate t-statistic and p-value using atop() for vertical stacking
  box <- box + annotate(
    "text",
    x = family_name,
    y = 3.6,  # Adjust y position as needed
    label = paste0("atop(italic(t) == ", t_stat, ", italic(p) == ", p_label_str, ")"),
    size = 1.5,
    family = "Arial",
    fontface = "bold",
    hjust = 0.5,
    parse = TRUE
  )
}

boxes1 <- grid.arrange(box, ncol = 1)
ggsave("~/projects/proj1/data/paper_images/Fig3/famComp_vert.png",
       plot = boxes1, width = 2, height = 3, dpi = 500)
ggsave("~/projects/proj1/data/paper_images/Fig3/famComp_vert.pdf",
       plot = boxes1, width = 2, height = 3, dpi = 500, device = cairo_pdf)
