#Analysis of VISION genes
getwd()
setwd("C:/2025/R_Analysis/AnalysisOfGeneExpressionDATA/RT-qPCR_RelativeGeneExpression")
# Load libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)

# Load data
df <- read.csv("visualGenes_final.csv")
df
# Compute mean HK
df <- df %>%
  mutate(HK_mean = (HK1_Ct + HK2_Ct)/2)

# ΔCt
df <- df %>%
  mutate(Delta_Ct = Target_Ct - HK_mean)

# Mean ΔCt in control
control_means <- df %>%
  filter(group == "control") %>%
  group_by(goi) %>%
  summarise(mean_Delta_Ct = mean(Delta_Ct))

# ΔΔCt and Fold Change
df <- df %>%
  left_join(control_means, by="goi") %>%
  mutate(
    DeltaDelta_Ct = Delta_Ct - mean_Delta_Ct,
    FoldChange = 2^(-DeltaDelta_Ct),
    Regulation = ifelse(FoldChange > 1, "Upregulated", "Downregulated")
  )

# Wilcoxon test per gene and group
genes <- unique(df$goi)
groups <- unique(df$group)
groups <- groups[groups != "control"]

results <- data.frame()

for (g in genes) {
  for (trt in groups) {
    x <- df %>% filter(goi==g, group=="control") %>% pull(Delta_Ct)
    y <- df %>% filter(goi==g, group==trt) %>% pull(Delta_Ct)
    
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

# Save results
write.csv(df, "qPCR_relative_expressionFINALBOXPLOT.csv", row.names=FALSE)
write.csv(results, "qPCR_wilcoxon_resultsFINALBOXPLOT.csv", row.names=FALSE)

# Merge p-values back for plotting
df_plot <- df %>%
  left_join(
    results %>%
      mutate(group = gsub(" vs Control","",Comparison)) %>%
      select(Gene, group, p_value),
    by = c("goi"="Gene","group")
  )

# Convert p-value to significance stars
df_plot <- df_plot %>%
  mutate(
    p_star = case_when(
      p_value <=0.001 ~ "***",
      p_value <=0.01 ~ "**",
      p_value <=0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# Aggregate Relative Expression per group/gene
summary_df <- df_plot %>%
  group_by(goi, group) %>%
  summarise(
    mean_FC = mean(FoldChange),
    sd_FC = sd(FoldChange),
    p_star=first(p_star)
  )

# Barplot per gene WITHOUT error bars, WITHOUT "ns", and WITH dotted reference line at 1

# Keep only rows with significant stars
summary_df_signif <- summary_df %>%
  filter(p_star %in% c("*","**","***"))

#BOXPLOT

# Boxplot of individual FoldChange values with significance stars and dotted line at 1

# Get max FoldChange and p_star per gene and group (only significant)
df_plot_signif <- df_plot %>%
  filter(p_star %in% c("*","**","***")) %>%
  group_by(goi, group, p_star) %>%
  summarise(
    max_FC = max(FoldChange)
  )

ggplot(df_plot, aes(x=group, y=FoldChange, fill=group)) +
  geom_boxplot(width=0.6, outlier.shape=NA) +
  geom_jitter(width=0.2, alpha=0.5, size=1) +
  geom_hline(yintercept=1, linetype="dotted", color="black") +
  # Add stars for significant comparisons
  geom_text(
    data=df_plot_signif,
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

### FURTHER PLOT EDITS to remove control from plot
# Filter out 'control' group only for plotting
df_plot_no_control <- df_plot %>%
  filter(group != "control")

# Get significance annotations (excluding control from display but not from stats)
df_plot_signif <- df_plot_no_control %>%
  filter(p_star %in% c("*", "**", "***")) %>%
  group_by(goi, group, p_star) %>%
  summarise(max_FC = max(FoldChange), .groups = "drop")

# Final boxplot
ggplot(df_plot_no_control, aes(x = group, y = FoldChange, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
  geom_text(
    data = df_plot_signif,
    aes(label = p_star, y = max_FC + 0.5),
    size = 5
  ) +
  facet_wrap(~goi, scales = "free_y") +
  labs(
    title = "Relative Expression per Gene)",
    y = "Relative Expression (2^-ΔΔCt)",
    x = "Group",
    caption = "* = p≤0.05, ** = p≤0.01, *** = p≤0.001"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

#  DONE

