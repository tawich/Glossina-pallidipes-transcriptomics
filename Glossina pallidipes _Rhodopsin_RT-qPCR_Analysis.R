getwd()
setwd("C:/2025/R_Analysis/AnalysisOfGeneExpressionDATA/RT-qPCR_RelativeGeneExpression")
# Load libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(ggrepel)

# Load data
dfvg <- read.csv("visualGenes_final.csv")
# Compute mean HKGs
dfvg <- dfvg %>%
  mutate(HK_mean = sqrt(HK1_Ct * HK2_Ct))

# ΔCt
dfvg <- dfvg %>%
  mutate(Delta_Ct = Target_Ct - HK_mean)

# Mean ΔCt in control
control_means <- dfvg %>%
  filter(group == "control") %>%
  group_by(goi) %>%
  summarise(mean_Delta_Ct = mean(Delta_Ct))

View(dfvg)
# ΔΔCt and Fold Change
dfvg <- dfvg %>%
  left_join(control_means, by="goi") %>%
  mutate(
    DeltaDelta_Ct = Delta_Ct - mean_Delta_Ct,
    FoldChange = 2^(-DeltaDelta_Ct),
    Regulation = ifelse(FoldChange > 1, "Upregulated", "Downregulated")
  )
#NORMALITY
#After normality check using shapiro, 1 data not normal hence used wilcoxon/man whitney U test

# Wilcoxon test per gene and group
genes <- unique(dfvg$goi)
groups <- unique(dfvg$group)
groups <- groups[groups != "control"]

results <- data.frame()

for (g in genes) {
  for (trt in groups) {
    x <- dfvg %>% filter(goi==g, group=="control") %>% pull(Delta_Ct)
    y <- dfvg %>% filter(goi==g, group==trt) %>% pull(Delta_Ct)
    
    if(length(x) >=1 & length(y) >=1){
      wtest <- wilcox.test(x, y)
      pval <- wtest$p.value
      median_control <- median(x)
      median_trt <- median(y)
      direction <- ifelse(median_trt < median_control, "Upregulated", "Downregulated")
      
      results <- rbind(
        results,
        data.frame(
          Gene=g,
          Comparison=paste(trt, "vs Control"),
          p_value=round(pval,4),
          Median_Control=round(median_control,2),
          Median_Treatment=round(median_trt,2),
          Regulation=direction
        )
      )
    }
  }
}

# Saving results for use in statistics and putting significant stats
write.csv(dfvg, "visionGenesResults20082025.csv", row.names=FALSE)
write.csv(results, "wilcoxonVisionGenesResults20082025.csv", row.names=FALSE)
getwd()
# Merging p-values back for plotting
dfvg_plot <- dfvg %>%
  left_join(
    results %>%
      mutate(group = gsub(" vs Control","",Comparison)) %>%
      select(Gene, group, p_value),
    by = c("goi"="Gene","group")
  )

# Converting p-value to significance stars
dfvg_plot <- dfvg_plot %>%
  mutate(
    p_star = case_when(
      p_value <=0.001 ~ "***",
      p_value <=0.01 ~ "**",
      p_value <=0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# Aggregating Relative Expression per group
summary_dfvg <- dfvg_plot %>%
  group_by(goi, group) %>%
  summarise(
    mean_FC = mean(FoldChange),
    sd_FC = sd(FoldChange),
    p_star=first(p_star)
  )

# Generating boxplots per gene and adding a dotted reference line at 1

# filtering to keep only rows with significant stars
summary_dfvg_signif <- summary_dfvg %>%
  filter(p_star %in% c("*","**","***"))

# Generating boxplot of individual FoldChange values

# Get max FoldChange and p_star per gene and group
dfvg_plot_signif <- dfvg_plot %>%
  filter(p_star %in% c("*","**","***")) %>%
  group_by(goi, group, p_star) %>%
  summarise(
    max_FC = max(FoldChange)
  )

ggplot(dfvg_plot, aes(x=group, y=FoldChange, fill=group)) +
  geom_boxplot(width=0.6, outlier.shape=NA) +
  geom_jitter(width=0.2, alpha=0.5, size=1) +
  geom_hline(yintercept=1, linetype="dotted", color="black") +
  
  # Adding stars for significant comparisons
  geom_text(
    data=dfvg_plot_signif,
    aes(label=p_star, y=max_FC + 0.5),
    size=5
  ) +
  facet_wrap(~goi, scales="free_y") +
  labs(
    title="Relative Expression per Gene and Group",
    y="Relative Expression (2^-ΔΔCt)",
    x="Group",
    caption="* = p<=0.05, ** = p<=0.01, *** = p<=0.001"
  ) +
  theme_minimal() +
  theme(legend.position="none")

### Focusing ONLY on the analysis on transcripts from G. pallidipes heads with mature trypanosomes: Filtering control and GO+/- from plot and to retain only GM+/+ group

dfvg_filtered <- dfvg_plot %>% filter(group == "GM+/+")

ggplot(dfvg_filtered, aes(x = group, y = FoldChange, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  # Label each point with its 'name' in red and smaller font
  #geom_text_repel(
  # aes(label = name),
  #size = 2.5,
  #color = "black",
  #max.overlaps = Inf
  #) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
  facet_wrap(~goi, scales = "free_y") +
  labs(
    title = "Relative Expression for GM+/+ Group",
    y = "Relative Expression (2^-ΔΔCt)",
    x = "Group",
    caption = "* = p≤0.05, ** = p≤0.01, *** = p≤0.001"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

## Filter to remainonly GM+/+ group

dfvg_gm <- dfvg_plot %>% filter(group == "GM+/+")

#Classifying samples into 'high' and 'low' based on 'name'
dfvg_gm <- dfvg_gm %>%
  mutate(
    Parasitemia = case_when(
      grepl("high", name, ignore.case = TRUE) ~ "High",
      grepl("low", name, ignore.case = TRUE) ~ "Low",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Parasitemia))

# Plotting boxplots by Parasitemia

## customize plot aesthetics
ggplot(dfvg_gm, aes(x = Parasitemia, y = FoldChange, fill = Parasitemia)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "blue") +
  facet_wrap(~goi, scales = "free_y") +
  labs(
    title = "Expression Comparison based on Mouthpart parasitemia levels",
    x = "Genes per Group",
    y = "Relative Expression (2^-ΔΔCt)"
  ) +
  scale_fill_manual(values = c(
    "Low" = "grey",
    #"Medium" = "black"(only input this if you are adding this group),
    "High" = "darkred"
  )) +
  theme_minimal() +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold", color = "black", size = 10)
  )

view(dfvg_gm)

# Loading sorted data

dfvg <- read.csv("Final_BoxplotData.csv")
View(dfvg)
# Assigning significance stars
get_stars <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("ns")
}

# Re-running normality test per group
stats <- dfvg %>%
  group_by(goi, Parasitemia) %>%
  summarise(
    shapiro_p = tryCatch(shapiro.test(FoldChange)$p.value, error = function(e) NA),
    test_p = if (!is.na(shapiro_p[1]) && shapiro_p[1] > 0.05) {
      t.test(FoldChange, mu = 1)$p.value
    } else {
      wilcox.test(FoldChange, mu = 1)$p.value
    },
    stars = get_stars(test_p),
    y_pos = max(FoldChange, na.rm = TRUE) * 1.05,
    .groups = "drop"
  )

# Plotting all boxplots and facet wrapping

ggplot(dfvg, aes(x = Parasitemia, y = FoldChange, fill = Parasitemia)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), alpha = 0.7) +
  geom_text(data = stats, aes(x = Parasitemia, y = y_pos, label = stars, group = Parasitemia),
            position = position_dodge(width = 0.8), vjust = 0) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black", linewidth = 0.6) +
  facet_wrap(~ goi, scales = "free_y") +
  scale_fill_manual(values = c("high" = "lightblue", "low" = "lightgreen")) +
  labs(title = "Gene Expression per Gene",
       y = "Relative Expression (2^-ΔΔCt)", x = "Groups") +
  theme_minimal() +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 10, color = "darkred")
  )

##DONE##





