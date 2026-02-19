################################################################################
# Author: Danielle Rivera
# 
# This code was used in:
# Rivera D, et al. 2026. Genomic data reveal historical introgression and 
# genetic erosion across an imperiled radiation of Caribbean Spondylurus skinks.
# Evolutionary Journal of the Linnean Society. In Review.
################################################################################

setwd("insert directory here")

# Load required libraries
library(detectRUNS)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
library(stringr)

# Set minimal theme for all plots
theme_conservation <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray92", linewidth = 0.25),
    panel.border = element_rect(fill = NA, color = "gray60", linewidth = 0.5),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "gray40", linewidth = 0.4),
    axis.text = element_text(color = "gray20"),
    legend.position = "top",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 11, hjust = 0),
    plot.subtitle = element_text(size = 9, color = "gray40", hjust = 0),
    strip.background = element_rect(fill = "gray95", color = "gray60", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 9),
    panel.spacing = unit(0.5, "lines")
  )

# UPDATED COLOR PALETTE
species_colors <- c("S. nitidus" = "#003c86",      # dark blue
                    "S. culebrae" = "#009ffa",     # light blue
                    "S. powelli" = "#028169")      # teal

island_colors <- c("Tintamarre" = "#028169",    # Dark teal
                   "Anguilla" = "#66c2a5",      # Medium teal  
                   "St. Barts" = "#a6d96a")     # Light green

# ============================================================================
# LOAD DATA FOR SPECIES
# ============================================================================

# S. nitidus
WINDOW_SIZE <- 100
THRESHOLD <- 0.05
MIN_SNP <- 100
ROH_ET <- FALSE
MAX_OPP_RUN <- 3
MAX_MISS_RUN <- 0

cat("Processing S. nitidus...\n")
map_nit <- "input_files/Spon_nitidus_filtered_for_ROH.map"
ped_nit <- "input_files/Spon_nitidus_filtered_for_ROH.ped"
runs_nit <- slidingRUNS.run(ped_nit, map_nit,
                            windowSize = WINDOW_SIZE,
                            threshold = THRESHOLD,
                            minSNP = MIN_SNP,
                            ROHet = ROH_ET,
                            maxOppRun = MAX_OPP_RUN,
                            maxMissRun = MAX_MISS_RUN)
froh_nit <- Froh_inbreeding(runs_nit, mapFile = map_nit, genome_wide = TRUE)
froh_nit$species <- "S. nitidus"
runs_nit$species <- "S. nitidus"
cat("  ROH detected:", nrow(runs_nit), "\n\n")
  
# S. culebrae
WINDOW_SIZE <- 150
MIN_SNP <- 150
MAX_OPP_RUN <- 5
MAX_MISS_RUN <- 0

cat("Processing S. culebrae...\n")
map_cul <- "input_files/Spon_culebrae_filtered_for_ROH.map"
ped_cul <- "input_files/Spon_culebrae_filtered_for_ROH.ped"
runs_cul <- slidingRUNS.run(ped_cul, map_cul,
                            windowSize = WINDOW_SIZE,
                            threshold = THRESHOLD,
                            minSNP = MIN_SNP,
                            ROHet = ROH_ET,
                            maxOppRun = MAX_OPP_RUN,
                            maxMissRun = MAX_MISS_RUN)
froh_cul <- Froh_inbreeding(runs_cul, mapFile = map_cul, genome_wide = TRUE)
froh_cul$species <- "S. culebrae"
runs_cul$species <- "S. culebrae"
cat("  ROH detected:", nrow(runs_cul), "\n\n")

# S. powelli
WINDOW_SIZE <- 50
MIN_SNP <- 75
MAX_OPP_RUN <- 2
MAX_MISS_RUN <- 0
  
cat("Processing S. powelli...\n")
map_pow <- "input_files/Spon_powelli_filtered_strict.map"
ped_pow <- "input_files/Spon_powelli_filtered_strict.ped"
runs_pow <- slidingRUNS.run(ped_pow, map_pow,
                            windowSize = WINDOW_SIZE,
                            threshold = THRESHOLD,
                            minSNP = MIN_SNP,
                            ROHet = ROH_ET,
                            maxOppRun = MAX_OPP_RUN,
                            maxMissRun = MAX_MISS_RUN)
froh_pow <- Froh_inbreeding(runs_pow, mapFile = map_pow, genome_wide = TRUE)
froh_pow$species <- "S. powelli"
runs_pow$species <- "S. powelli"
cat("  ROH detected:", nrow(runs_pow), "\n\n")

# ADD ISLAND INFORMATION FOR S. POWELLI
add_island_info <- function(data) {
  data %>%
    mutate(
      island = case_when(
        str_detect(id, "27477") ~ "Tintamarre",
        str_detect(id, "267291") ~ "Anguilla",
        str_detect(id, "MNHN|0844|0843") ~ "St. Barts",
        TRUE ~ "Unknown"
      )
    )
}

froh_pow <- add_island_info(froh_pow)
runs_pow <- add_island_info(runs_pow)

cat("Island assignments for S. powelli:\n")
print(table(froh_pow$island))
cat("\n")

# Combine data
froh_combined <- rbind(
  froh_nit[, c("id", "Froh_genome", "species")],
  froh_cul[, c("id", "Froh_genome", "species")],
  froh_pow[, c("id", "Froh_genome", "species")])

runs_combined <- rbind(
  runs_nit[, c("id", "chrom", "nSNP", "lengthBps", "species")],
  runs_cul[, c("id", "chrom", "nSNP", "lengthBps", "species")],
  runs_pow[, c("id", "chrom", "nSNP", "lengthBps", "species")])

# Calculate ROH length in Mb
runs_combined$lengthMb <- runs_combined$lengthBps / 1e6

# ADD ROH CATEGORIES 
runs_combined <- runs_combined %>%
  mutate(roh_category = case_when(
    lengthMb < 1 ~ "<1 Mb",              
    lengthMb >= 1 & lengthMb < 4 ~ "1-4 Mb",     
    lengthMb >= 4 & lengthMb < 8 ~ "4-8 Mb",     
    lengthMb >= 8 & lengthMb < 16 ~ "8-16 Mb",
    lengthMb >= 16 ~ ">16 Mb"            
  )) %>%
  mutate(roh_category = factor(roh_category, 
                               levels = c("<1 Mb", "1-4 Mb", "4-8 Mb", "8-16 Mb", ">16 Mb")))

cat("ROH category distribution:\n")
print(table(runs_combined$roh_category, runs_combined$species))
cat("\n")

# ============================================================================
# FIGURE 1: F_ROH COMPARISON
# ============================================================================

# Calculate stats for mean bar
froh_stats <- froh_combined %>%
  group_by(species) %>%
  summarise(
    mean_val = mean(Froh_genome, na.rm = TRUE),
    median_val = median(Froh_genome, na.rm = TRUE),
    n = n()
  )

cat("F_ROH statistics by species:\n")
print(froh_stats)
cat("\n")

# ============================================================================
# OUTLIER TESTING - Test if Spon5 and Spon13 are significantly different
# ============================================================================

cat("OUTLIER ANALYSIS\n")

# Test Spon5 (S. culebrae outlier)
spon5_froh <- froh_combined %>% filter(id == "Sculebrae_Spon5" | id == "Spon5")
other_cul_froh <- froh_combined %>% 
  filter(species == "S. culebrae" & id != "Sculebrae_Spon5" & id != "Spon5")

if(nrow(spon5_froh) > 0) {
  cat("S. culebrae Outlier Test (Spon5):\n")
  cat("  Spon5 F_ROH:", round(spon5_froh$Froh_genome[1], 4), "\n")
  cat("  Other S. culebrae mean F_ROH:", round(mean(other_cul_froh$Froh_genome), 4), "\n")
  
  # Grubbs test for outlier
  if(requireNamespace("outliers", quietly = TRUE)) {
    library(outliers)
    grubbs_test <- grubbs.test(froh_combined$Froh_genome[froh_combined$species == "S. culebrae"])
    cat("  Grubbs test p-value:", format.pval(grubbs_test$p.value, digits = 3), "\n")
  }
  
  # Z-score
  cul_mean <- mean(other_cul_froh$Froh_genome)
  cul_sd <- sd(other_cul_froh$Froh_genome)
  spon5_zscore <- (spon5_froh$Froh_genome[1] - cul_mean) / cul_sd
  cat("  Z-score:", round(spon5_zscore, 2), "\n")
  cat("  Interpretation:", ifelse(abs(spon5_zscore) > 3, "EXTREME outlier (>3 SD)", 
                                  ifelse(abs(spon5_zscore) > 2, "Significant outlier (>2 SD)", 
                                         "Not a strong outlier")), "\n\n")
}

# Test Spon13 (S. nitidus phylogenetic outlier)
spon13_froh <- froh_combined %>% filter(id == "Spon13")
other_nit_froh <- froh_combined %>% 
  filter(species == "S. nitidus" & id != "Spon13")

if(nrow(spon13_froh) > 0) {
  cat("S. nitidus Outlier Test (Spon13):\n")
  cat("  Spon13 F_ROH:", round(spon13_froh$Froh_genome[1], 4), "\n")
  cat("  Other S. nitidus mean F_ROH:", round(mean(other_nit_froh$Froh_genome), 4), "\n")
  
  # Grubbs test for outlier
  if(requireNamespace("outliers", quietly = TRUE)) {
    grubbs_test <- grubbs.test(froh_combined$Froh_genome[froh_combined$species == "S. nitidus"])
    cat("  Grubbs test p-value:", format.pval(grubbs_test$p.value, digits = 3), "\n")
  }
  
  # Z-score
  nit_mean <- mean(other_nit_froh$Froh_genome)
  nit_sd <- sd(other_nit_froh$Froh_genome)
  spon13_zscore <- (spon13_froh$Froh_genome[1] - nit_mean) / nit_sd
  cat("  Z-score:", round(spon13_zscore, 2), "\n")
  cat("  Interpretation:", ifelse(abs(spon13_zscore) > 3, "EXTREME outlier (>3 SD)", 
                                  ifelse(abs(spon13_zscore) > 2, "Significant outlier (>2 SD)", 
                                         "Not a strong outlier")), "\n\n")
}

# Create a display category for coloring
froh_combined <- froh_combined %>%
  mutate(display_group = ifelse(id == "Spon13" | id == "Snitidus_Spon13", 
                                "Spon13", 
                                species))

# Extend color palette to include Spon13
display_colors <- c(
  "S. nitidus" = "#003c86",
  "S. culebrae" = "#009ffa", 
  "S. powelli" = "#028169",
  "Spon13" = "#ff70fd"
)

p1 <- ggplot(froh_combined, aes(x = species, y = Froh_genome)) +
  # Individual points - colored by display_group instead of species
  geom_jitter(aes(fill = display_group), width = 0.15, 
              size = 5, shape = 21, color = "black", stroke = 0.5) +
  # Mean horizontal bar
  geom_segment(data = froh_stats,
               aes(x = as.numeric(factor(species)) - 0.35,
                   xend = as.numeric(factor(species)) + 0.35,
                   y = mean_val, yend = mean_val),
               linewidth = 1.2, color = "gray20") +
  scale_fill_manual(values = display_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)),
                     labels = number_format(accuracy = 0.01)) +
  labs(x = NULL, 
       y = expression(bold(F[ROH])),
       title = "A. Genomic Inbreeding Levels",
       subtitle = "Proportion of genome in runs of homozygosity") +
  theme_conservation +
  theme(legend.position = "none",
        plot.title = element_text(size = 25, face = "bold"),
        plot.subtitle = element_text(size = 15),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(face = "italic", size = 15, color = "gray20"),
        axis.text.y = element_text(size = 15, color = "gray20"))

p1
# ============================================================================
# FIGURE 2: RECENT INBREEDING INDICATOR
# ============================================================================

# Calculate proportion of long ROH (>5 Mb) per individual
recent_inbreeding <- runs_combined %>%
  mutate(is_recent = lengthMb > 5) %>%
  group_by(id, species) %>%
  summarise(
    prop_recent = sum(is_recent) / n(),
    n_recent = sum(is_recent),
    n_total = n(),
    mean_roh_length = mean(lengthMb),
    .groups = 'drop'
  )
# Add display category for coloring Spon13
recent_inbreeding <- recent_inbreeding %>%
  mutate(display_group = ifelse(id == "Spon13" | id == "Snitidus_Spon13", 
                                "Spon13", 
                                species))
  
  # Color palette
  display_colors <- c(
    "S. nitidus" = "#003c86",
    "S. culebrae" = "#009ffa", 
    "S. powelli" = "#028169",
    "Spon13" = "#ff70fd"
  )
  
  # Calculate species means for annotation
  recent_means <- recent_inbreeding %>%
    group_by(species) %>%
    summarise(mean_prop = mean(prop_recent))
  
  max_prop <- max(recent_inbreeding$prop_recent)
  
  if(max_prop < 0.01) {
    p2 <- ggplot(recent_inbreeding, aes(x = species, y = mean_roh_length)) +
      geom_boxplot(aes(fill = species), alpha = 0.4, outlier.shape = NA,
                   color = "gray40", linewidth = 0.5, width = 0.6) +
      geom_jitter(aes(fill = display_group), width = 0.12,
                  size = 5, shape = 21, color = "black", stroke = 0.5) +
      scale_fill_manual(values = display_colors) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
      labs(x = NULL,
           y = "Mean ROH Length (Mb)",
           title = "B. ROH Length Distribution",
           subtitle = "Longer ROH indicate more recent inbreeding events") +
      theme_conservation +
      theme(legend.position = "none",
            plot.title = element_text(size = 25, face = "bold"),
            plot.subtitle = element_text(size = 15),
            axis.title.y = element_text(size = 18, face = "bold"),
            axis.text.x = element_text(face = "italic", size = 15, color = "gray20"),
            axis.text.y = element_text(size = 15, color = "gray20"))
  } else {
    y_format <- if(max_prop < 0.1) percent_format(accuracy = 0.1) else percent_format(accuracy = 1)
    
    p2 <- ggplot(recent_inbreeding, aes(x = species, y = prop_recent)) +
      geom_boxplot(aes(fill = species), outlier.shape = NA,
                   color = "gray40", linewidth = 0.5, width = 0.6) +
      geom_jitter(aes(fill = display_group), width = 0.12,
                  size = 5, shape = 21, color = "black", stroke = 0.5) +
      scale_fill_manual(values = display_colors) +
      scale_y_continuous(labels = y_format,
                         expand = expansion(mult = c(0.05, 0.1))) +
      labs(x = NULL,
           y = "Proportion of ROH > 5 Mb",
           title = "B. Recent Inbreeding Signal",
           subtitle = "Higher values indicate more recent population bottlenecks") +
      theme_conservation +
      theme(legend.position = "none",
            plot.title = element_text(size = 25, face = "bold"),
            plot.subtitle = element_text(size = 15),
            axis.title.y = element_text(size = 18, face = "bold"),
            axis.text.x = element_text(face = "italic", size = 15, color = "gray20"),
            axis.text.y = element_text(size = 15, color = "gray20"))
  }
  
  p2
# ============================================================================
# S. POWELLI ISLAND TESTING
# ============================================================================

# Test for differences between islands
froh_pow_detailed <- froh_pow %>%
  filter(!is.na(island) & island != "Unknown")

if(nrow(froh_pow_detailed) >= 3) {
  cat("S. POWELLI ISLAND COMPARISON\n")
  
  # Summary by island
  island_stats <- froh_pow_detailed %>%
    group_by(island) %>%
    summarise(
      n = n(),
      mean_FROH = mean(Froh_genome),
      sd_FROH = sd(Froh_genome),
      median_FROH = median(Froh_genome)
    )
  
  print(island_stats)
  cat("\n")
  
  # Statistical test if we have enough samples per island
  if(n_distinct(froh_pow_detailed$island) > 1 && 
     all(table(froh_pow_detailed$island) >= 1)) {
    
    if(n_distinct(froh_pow_detailed$island) == 2) {
      # Wilcoxon test for 2 groups
      test_result <- wilcox.test(Froh_genome ~ island, data = froh_pow_detailed)
      cat("Wilcoxon rank-sum test:\n")
      cat("  W =", round(test_result$statistic, 2), "\n")
      cat("  p-value =", format.pval(test_result$p.value, digits = 3), "\n")
    } else {
      # Kruskal-Wallis for >2 groups
      test_result <- kruskal.test(Froh_genome ~ island, data = froh_pow_detailed)
      cat("Kruskal-Wallis test:\n")
      cat("  Chi-squared =", round(test_result$statistic, 3), "\n")
      cat("  p-value =", format.pval(test_result$p.value, digits = 3), "\n")
    }
    
    cat("\nInterpretation:", 
        ifelse(test_result$p.value < 0.05, 
               "Significant difference between islands",
               "No significant difference between islands"), "\n")
  } else {
    cat("Insufficient data for statistical comparison between islands\n")
  }
}

# ============================================================================
# COMBINE PLOTS AND SAVE
# ============================================================================

# Create multi-panel figure - VERTICAL ARRANGEMENT
combined_fig_main <- p1 / p2
combined_fig_main
# Save individual plots
ggsave("Figure_FROH_comparison.pdf", p1, width = 6, height = 5, dpi = 300)
ggsave("Figure_FROH_comparison.png", p1, width = 6, height = 5, dpi = 300)
ggsave("Figure_FROH_comparison.svg", p1, width = 6, height = 5, dpi = 300)

ggsave("Figure_recent_inbreeding.pdf", p2, width = 6, height = 5, dpi = 300)
ggsave("Figure_recent_inbreeding.png", p2, width = 6, height = 5, dpi = 300)
ggsave("Figure_recent_inbreeding.svg", p2, width = 6, height = 5, dpi = 300)

# Combined 2-panel figure - VERTICAL
ggsave("Figure_ROH_combined.pdf", combined_fig_main, width = 8, height = 12, dpi = 300)
ggsave("Figure_ROH_combined.png", combined_fig_main, width = 8, height = 12, dpi = 300)
ggsave("Figure_ROH_combined.svg", combined_fig_main, width = 8, height = 12, dpi = 300)

# S. powelli island detail figures
if(!is.null(p1b) && !is.null(p3b)) {
  combined_fig_islands <- p1b / p3b
  ggsave("Figure_Spowelli_islands.pdf", combined_fig_islands, width = 7, height = 8, dpi = 300)
  ggsave("Figure_Spowelli_islands.png", combined_fig_islands, width = 7, height = 8, dpi = 300)
  ggsave("Figure_Spowelli_islands.svg", combined_fig_islands, width = 7, height = 8, dpi = 300)
}

# ============================================================================
# SUMMARY STATISTICS TABLE
# ============================================================================

cat("Calculating summary statistics...\n")

# Calculate total ROH per individual
ind_metrics <- runs_combined %>%
  group_by(id, species) %>%
  summarise(
    total_roh_mb = sum(lengthMb),
    n_roh = n(),
    .groups = 'drop'
  )

# Calculate quartiles for each species for reference
quartiles <- ind_metrics %>%
  group_by(species) %>%
  summarise(
    q25 = quantile(total_roh_mb, 0.25),
    q50 = quantile(total_roh_mb, 0.50),
    q75 = quantile(total_roh_mb, 0.75)
  )

# Multi-species summary
summary_table <- froh_combined %>%
  group_by(species) %>%
  summarise(
    n_individuals = n(),
    mean_FROH = round(mean(Froh_genome, na.rm = TRUE), 4),
    sd_FROH = round(sd(Froh_genome, na.rm = TRUE), 4),
    median_FROH = round(median(Froh_genome, na.rm = TRUE), 4),
    min_FROH = round(min(Froh_genome, na.rm = TRUE), 4),
    max_FROH = round(max(Froh_genome, na.rm = TRUE), 4)
  )

# ROH characteristics by species
roh_summary <- runs_combined %>%
  group_by(species) %>%
  summarise(
    total_ROH = n(),
    mean_ROH_per_ind = round(n() / n_distinct(id), 1),
    mean_ROH_length_Mb = round(mean(lengthMb), 2),
    median_ROH_length_Mb = round(median(lengthMb), 2),
    prop_short_ROH = round(sum(lengthMb < 1) / n(), 3),    
    prop_long_ROH = round(sum(lengthMb > 8) / n(), 3)      
  )

# ROH by category
roh_by_category <- runs_combined %>%
  group_by(species, roh_category) %>%
  summarise(n = n(), .groups = 'drop')

# Combine summaries
final_summary <- left_join(summary_table, roh_summary, by = "species")

# Save summary tables
write.csv(final_summary, "ROH_summary_statistics.csv", row.names = FALSE)
write.csv(roh_by_category, "ROH_by_category.csv", row.names = FALSE)

# Save combined ROH and FROH data
write.csv(froh_combined, "Combined_FROH.csv", row.names = FALSE)
write.csv(runs_combined, "Combined_ROH.csv", row.names = FALSE)

# Save individual-level metrics (ROH count and total length per sample)
write.csv(ind_metrics, "Individual_ROH_metrics.csv", row.names = FALSE)

# Create comprehensive per-individual dataset (FROH + ROH metrics)
comprehensive_data <- froh_combined %>%
  left_join(ind_metrics, by = c("id", "species")) %>%
  mutate(
    mean_roh_length_mb = total_roh_mb / n_roh,
    froh_category = case_when(
      Froh_genome < 0.0625 ~ "Low (<0.0625)",
      Froh_genome >= 0.0625 & Froh_genome < 0.125 ~ "Moderate (0.0625-0.125)",
      Froh_genome >= 0.125 & Froh_genome < 0.25 ~ "High (0.125-0.25)",
      Froh_genome >= 0.25 ~ "Very High (>0.25)"
    )
  ) %>%
  arrange(species, -Froh_genome)

write.csv(comprehensive_data, "Comprehensive_Individual_Data.csv", row.names = FALSE)

cat("\nSaved combined data files:\n")
cat("  - Combined_FROH.csv (F_ROH per individual)\n")
cat("  - Combined_ROH.csv (all ROH segments)\n")
cat("  - Individual_ROH_metrics.csv (ROH count & total length per individual)\n")
cat("  - Comprehensive_Individual_Data.csv (F_ROH + ROH metrics merged)\n")

# ============================================================================
# PRINT RESULTS
# ============================================================================

cat("MULTI-SPECIES ROH ANALYSIS SUMMARY\n")
print(final_summary)

cat("ROH BY CATEGORY\n")
print(roh_by_category)

# ============================================================================
# STATISTICAL TESTS
# ============================================================================

cat("STATISTICAL TESTS\n")

# Test for F_ROH difference across species
test_result <- kruskal.test(Froh_genome ~ species, data = froh_combined)
cat("Kruskal-Wallis test for species F_ROH difference:\n")
cat("  Chi-squared =", round(test_result$statistic, 3), "\n")
cat("  p-value =", format.pval(test_result$p.value, digits = 3), "\n\n")

# Pairwise comparisons if significant
if(test_result$p.value < 0.05) {
  cat("Pairwise Wilcoxon tests (Bonferroni correction):\n")
  pairwise_result <- pairwise.wilcox.test(froh_combined$Froh_genome, 
                                          froh_combined$species,
                                          p.adjust.method = "bonferroni")
  print(pairwise_result)
  cat("\n")
}

# Recent inbreeding comparison
if(max_prop >= 0.01) {
  recent_test <- kruskal.test(prop_recent ~ species, data = recent_inbreeding)
  cat("Kruskal-Wallis test for recent inbreeding (prop ROH >5 Mb):\n")
  cat("  Chi-squared =", round(recent_test$statistic, 3), "\n")
  cat("  p-value =", format.pval(recent_test$p.value, digits = 3), "\n\n")
}

# Total ROH burden comparison
burden_test <- kruskal.test(total_roh_mb ~ species, data = ind_metrics)
cat("Kruskal-Wallis test for total ROH burden:\n")
cat("  Chi-squared =", round(burden_test$statistic, 3), "\n")
cat("  p-value =", format.pval(burden_test$p.value, digits = 3), "\n\n")


cat("ANALYSIS COMPLETE\n")

