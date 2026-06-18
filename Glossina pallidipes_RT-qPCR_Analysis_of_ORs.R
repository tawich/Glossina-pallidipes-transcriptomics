# ANALYSIS OF RT-qPCR OLFACTORY DATA
getwd()
setwd("C:/2025/R_Analysis/AnalysisOfGeneExpressionDATA/RT-qPCR_RelativeGeneExpression")

# Load libraries
library(dplyr)
library(ggplot2)
library(ggpubr)

# Load your data
df <- read.csv("Final_BoxplotData.csv")
View(df)

# Function to assign stars
get_stars <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("ns")
}

# Run normality test + appropriate test per group
stats <- df %>%
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

# Assigning smaller size for 'ns', larger for stars
stats$size <- ifelse(stats$stars == "ns", 3, 5)

# Plot
ggplot(df, aes(x = Parasitemia, y = FoldChange, fill = Parasitemia)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), alpha = 0.7) +
  
  # Significance labels with custom size and vertical position
  geom_text(data = stats,
            aes(x = Parasitemia, y = y_pos, label = stars, group = Parasitemia, size = size),
            position = position_dodge(width = 0.8),
            vjust = -0.2, color = "red", show.legend = FALSE) +
  
  # Reference line at FoldChange = 1
  geom_hline(yintercept = 1, linetype = "dotted", color = "black", linewidth = 0.6) +
  
  facet_wrap(~ goi, scales = "free_y") +
  
  scale_fill_manual(values = c("high" = "red", "low" = "green")) +
  scale_size_identity() +
  
  labs(
    title = "",
    y = "Relative Expression (2^-ΔΔCt)",
    x = "Groups"
  ) +
  
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    # Remove all grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    # Normal axis text, bold axis titles
    axis.text.x = element_text(face = "plain", color = "black", size = 10),
    axis.text.y = element_text(face = "plain", color = "black", size = 10),
    axis.title.x = element_text(face = "bold", color = "black", size = 12),
    axis.title.y = element_text(face = "bold", color = "black", size = 12),
    
    # Solid axis lines
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    axis.line.x = element_line(color = "black", linewidth = 0.8),
    
    # Tick marks
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    
    # Facet strip text (bold, no background)
    strip.text = element_text(face = "bold", size = 9, color = "black"),
    strip.background = element_blank(),
    
    # Legend
    legend.position = "right"
  )

