#!/bin/bash

# Genetic load analysis - 85% sample missing data threshold
# Calculates realized, masked, and total load per sample
# Produces species-level and temporal summaries for Snitidus and Smonae

set -e

echo "=== GENETIC LOAD ANALYSIS ==="
echo "Started at: $(date)"
echo ""

# The site missingness threshold is passed as an argument 
# Usage: bash run_genetic_load_analysis.sh 30site
SITE_THRESHOLD="${1:-30site}"
echo "Site missingness threshold: $SITE_THRESHOLD"
echo ""

# --- INPUT FILE: coding_variants_85pct_${SITE_THRESHOLD}.vcf.gz ---

# Only variants in coding regions (CDS) are retained so that the load
# calculation reflects functional (expressed) genetic burden rather than
# background neutral variation.

VCF_85="coding_variants_85pct_${SITE_THRESHOLD}.vcf.gz"

if [ ! -f "$VCF_85" ]; then
    echo "ERROR: VCF not found: $VCF_85"
    echo "Available VCF files:"
    ls -1 coding_variants_*.vcf.gz 2>/dev/null || echo "  None found"
    exit 1
fi

echo "Input VCF: $VCF_85"
echo ""

# Write the R analysis script to disk, then run it
cat > genetic_load_analysis.R << 'RSCRIPT'

args <- commandArgs(trailingOnly = TRUE)
site_threshold <- args[1]

cat("=== GENETIC LOAD ANALYSIS ===\n")
cat("Site threshold:", site_threshold, "\n\n")

# --- FUNCTION: Read a gzipped or plain VCF into a list ---
# Returns: header lines, fixed fields (CHROM-INFO), genotype matrix, sample names
read_vcf_base <- function(vcf_file) {
    cat("Reading VCF:", vcf_file, "\n")
    
    con <- if (grepl("\\.gz$", vcf_file)) gzfile(vcf_file, "rt") else file(vcf_file, "rt")
    
    # Skip metadata lines; stop at the column header line
    header_lines <- c()
    while (TRUE) {
        line <- readLines(con, n = 1)
        if (length(line) == 0) break
        if (startsWith(line, "##")) {
            header_lines <- c(header_lines, line)
        } else if (startsWith(line, "#CHROM")) {
            column_line <- line
            break
        }
    }
    
    col_names <- strsplit(gsub("^#", "", column_line), "\t")[[1]]
    
    # Read variant data in chunks of 10,000 lines for memory efficiency
    all_data <- list()
    chunk_size <- 10000
    chunk_num <- 0
    repeat {
        chunk <- readLines(con, n = chunk_size)
        if (length(chunk) == 0) break
        chunk_num <- chunk_num + 1
        if (chunk_num %% 10 == 0) cat(sprintf("  Read %d variants...\r", chunk_num * chunk_size))
        all_data <- c(all_data, lapply(chunk, function(x) strsplit(x, "\t")[[1]]))
    }
    close(con)
    
    cat(sprintf("\nTotal variants: %d\n", length(all_data)))
    
    data_matrix <- do.call(rbind, all_data)
    colnames(data_matrix) <- col_names
    
    list(
        header  = header_lines,
        fix     = data_matrix[, 1:8, drop = FALSE],
        gt      = data_matrix[, 9:ncol(data_matrix), drop = FALSE],
        samples = col_names[10:length(col_names)]
    )
}

# --- FUNCTION: Classify a variant by REF/ALT length ---
# SNP: single base change
# indel_frameshift: length difference not divisible by 3 (disrupts reading frame)
# indel_inframe: length difference divisible by 3 (preserves reading frame)
# complex: multi-base substitution of equal length
classify_variant <- function(ref, alt) {
    ref_len <- nchar(ref)
    alt_len <- nchar(alt)
    if (ref_len == 1 && alt_len == 1) return("SNP")
    if (ref_len != alt_len) {
        return(if (abs(alt_len - ref_len) %% 3 == 0) "indel_inframe" else "indel_frameshift")
    }
    return("complex")
}

# --- FUNCTION: Assign deleteriousness weight by variant type ---
# Weights reflect expected severity of functional impact:
#   frameshift = 1.0 (most severe: disrupts reading frame)
#   inframe    = 0.8 (moderate: alters protein length but preserves frame)
#   complex    = 0.7 (moderate: multi-nucleotide substitution)
#   SNP        = 0.5 (least severe: single amino acid change)
assign_weight <- function(variant_type) {
    switch(variant_type,
        "indel_frameshift" = 1.0,
        "indel_inframe"    = 0.8,
        "complex"          = 0.7,
        "SNP"              = 0.5,
        0.1
    )
}

# --- FUNCTION: Extract the GT field from a FORMAT/SAMPLE string pair ---
extract_gt <- function(format_str, sample_str) {
    if (is.na(sample_str) || sample_str == ".") return(NA)
    format_fields <- strsplit(format_str, ":")[[1]]
    sample_fields <- strsplit(sample_str, ":")[[1]]
    gt_idx <- which(format_fields == "GT")
    if (length(gt_idx) == 0 || gt_idx > length(sample_fields)) return(NA)
    sample_fields[gt_idx]
}

# Historical (museum) samples - used to assign temporal category
all_historical_samples <- c(
    'Sanegadae_242057', 'Sanegadae_242060', 'Smagnacruzae_242174',
    'Smonae_UPRMR365',  'Smonae_UPRMR695',  'Smonae_UPRM696',
    'Snitidus_UPRMR694', 'Snitidus_UPRMR710'
)

# --- FUNCTION: Process a VCF and calculate genetic load per sample ---
# Realized load  = sum of weights for homozygous alt genotypes (2/2)
# Masked load    = sum of (weights * 0.5) for heterozygous genotypes (0/1)
#                  multiplied by 0.5 to reflect partial dominance
# Total load     = realized + masked
# Load ratio     = realized / total (proportion of load currently expressed)
# Het ratio      = heterozygous sites / total genotyped sites
process_vcf <- function(vcf_file, threshold_label) {
    cat("\n--- Processing", threshold_label, "---\n")
    
    vcf_data    <- read_vcf_base(vcf_file)
    sample_names <- vcf_data$samples
    museum_samples       <- sample_names[sample_names %in% all_historical_samples]
    contemporary_samples <- sample_names[!sample_names %in% museum_samples]
    
    cat("Historical samples:", length(museum_samples),
        "| Contemporary samples:", length(contemporary_samples), "\n")
    
    # Build variant table with type and weight
    variant_info <- data.frame(
        REF = vcf_data$fix[, "REF"],
        ALT = vcf_data$fix[, "ALT"],
        stringsAsFactors = FALSE
    )
    variant_info$type   <- mapply(classify_variant, variant_info$REF, variant_info$ALT)
    variant_info$weight <- sapply(variant_info$type, assign_weight)
    
    cat("Variant composition:\n")
    for (vtype in names(table(variant_info$type))) {
        cat(sprintf("  %-18s: %8d (weight=%.1f)\n",
            vtype, sum(variant_info$type == vtype), assign_weight(vtype)))
    }
    
    # Calculate load metrics for each sample
    all_results <- data.frame()
    cat("Processing samples...\n")
    
    for (i in seq_along(sample_names)) {
        sample <- sample_names[i]
        if (i %% 10 == 0) cat(sprintf("  %d/%d...\r", i, length(sample_names)))
        
        sample_idx <- which(colnames(vcf_data$gt) == sample)
        gts <- mapply(extract_gt, vcf_data$gt[, 1], vcf_data$gt[, sample_idx],
                      USE.NAMES = FALSE)
        
        # Count alternate alleles per site (0 = ref/ref, 1 = het, 2 = hom alt, NA = missing)
        n_alt <- sapply(gts, function(x) {
            if (is.na(x) || x == "./.") return(NA)
            sum(strsplit(x, "/|\\|")[[1]] != "0", na.rm = TRUE)
        })
        
        realized_load  <- sum((n_alt == 2) * variant_info$weight, na.rm = TRUE)
        masked_load    <- sum((n_alt == 1) * variant_info$weight * 0.5, na.rm = TRUE)
        total_load     <- realized_load + masked_load
        total_hom      <- sum(n_alt == 2, na.rm = TRUE)
        total_het      <- sum(n_alt == 1, na.rm = TRUE)
        total_genotyped <- total_hom + total_het
        
        all_results <- rbind(all_results, data.frame(
            species           = sub("_.*", "", sample),
            sample            = sample,
            temporal_category = ifelse(sample %in% museum_samples, "Historical", "Contemporary"),
            n_variants        = sum(!is.na(n_alt)),
            n_genotyped       = total_genotyped,
            total_hom         = total_hom,
            total_het         = total_het,
            realized_load     = realized_load,
            masked_load       = masked_load,
            total_load        = total_load,
            load_ratio        = ifelse(total_load > 0, realized_load / total_load, NA),
            het_ratio         = ifelse(total_genotyped > 0, total_het / total_genotyped, NA),
            stringsAsFactors  = FALSE
        ))
    }
    
    cat("\n")
    return(all_results)
}

# --- RUN ANALYSIS ---
vcf_85   <- paste0("coding_variants_85pct_", site_threshold, ".vcf.gz")
results  <- process_vcf(vcf_85, "85%")

write.csv(results, paste0("genetic_load_results_85pct_", site_threshold, ".csv"),
          row.names = FALSE)

# --- SPECIES-LEVEL SUMMARY ---
cat("\n=== SPECIES-LEVEL SUMMARY ===\n")
species_summary <- aggregate(
    cbind(het_ratio, load_ratio, total_load, n_genotyped) ~ species,
    data = results,
    FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
)

for (i in 1:nrow(species_summary)) {
    sp <- species_summary$species[i]
    n  <- sum(results$species == sp)
    cat(sprintf("  %-15s (n=%2d): Het=%.3f, Load=%.3f\n",
        sp, n,
        species_summary$het_ratio[i, "mean"],
        species_summary$load_ratio[i, "mean"]))
}

# --- TEMPORAL ANALYSIS: Snitidus and Smonae ---
cat("\n=== TEMPORAL ANALYSIS ===\n")
temporal_species <- c("Snitidus", "Smonae")
temporal_results <- list()

for (sp in temporal_species) {
    cat("\nSpecies:", sp, "\n")
    
    sp_data   <- results[results$species == sp, ]
    hist_data <- sp_data[sp_data$temporal_category == "Historical", ]
    cont_data <- sp_data[sp_data$temporal_category == "Contemporary", ]
    
    cat("  Historical samples:", nrow(hist_data), "\n")
    cat("  Contemporary samples:", nrow(cont_data), "\n")
    
    if (nrow(hist_data) > 0 && nrow(cont_data) > 0) {
        hist_het  <- mean(hist_data$het_ratio,  na.rm = TRUE)
        cont_het  <- mean(cont_data$het_ratio,  na.rm = TRUE)
        hist_load <- mean(hist_data$load_ratio, na.rm = TRUE)
        cont_load <- mean(cont_data$load_ratio, na.rm = TRUE)
        delta_het  <- (cont_het  - hist_het)  * 100
        delta_load <-  cont_load - hist_load
        
        cat(sprintf("  Historical:   Het=%.3f, Load ratio=%.3f\n", hist_het,  hist_load))
        cat(sprintf("  Contemporary: Het=%.3f, Load ratio=%.3f\n", cont_het,  cont_load))
        cat(sprintf("  Change:       dHet=%+.1f pp, dLoad=%+.3f\n", delta_het, delta_load))
        
        temporal_results[[sp]] <- list(
            species = sp, hist_het = hist_het, cont_het = cont_het,
            hist_load = hist_load, cont_load = cont_load,
            delta_het = delta_het, delta_load = delta_load
        )
    } else {
        cat("  Not enough data for temporal comparison\n")
    }
}

# --- SAVE TEMPORAL SUMMARY ---
if (length(temporal_results) > 0) {
    temporal_df <- do.call(rbind, lapply(temporal_results, function(x) {
        data.frame(
            species          = x$species,
            historical_het   = x$hist_het,
            contemporary_het = x$cont_het,
            delta_het        = x$delta_het,
            historical_load  = x$hist_load,
            contemporary_load = x$cont_load,
            delta_load       = x$delta_load,
            stringsAsFactors = FALSE
        )
    }))
    write.csv(temporal_df,
              paste0("temporal_comparison_summary_", site_threshold, ".csv"),
              row.names = FALSE)
}

# Save species-level summary
write.csv(
    aggregate(cbind(het_ratio, load_ratio, total_load, n_genotyped) ~
              species + temporal_category,
              data = results,
              FUN = function(x) c(mean = mean(x, na.rm=TRUE),
                                  sd   = sd(x,   na.rm=TRUE),
                                  n    = length(x))),
    paste0("species_summary_", site_threshold, ".csv"),
    row.names = FALSE
)

# --- VISUALIZATION: Species-level load ratio ---
# Individual points per sample; circles = contemporary, stars = historical
# Each species gets a unique color; dashed orange line at load ratio 0.6
cat("\n=== CREATING PLOT ===\n")

if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)

# Order species alphabetically on x-axis
results$species <- factor(results$species, levels = sort(unique(results$species)))

png(paste0("genetic_load_", site_threshold, ".png"), width = 1400, height = 900, res = 150)

p <- ggplot(results, aes(x = species, y = load_ratio, color = species)) +
    # Contemporary samples: filled circles
    geom_point(data = subset(results, temporal_category == "Contemporary"),
               shape = 16, size = 3, alpha = 0.85) +
    # Historical samples: stars
    geom_point(data = subset(results, temporal_category == "Historical"),
               shape = 11, size = 4, alpha = 0.85) +
    # Dashed reference line at 0.6
    geom_hline(yintercept = 0.6, linetype = "dashed", color = "orange", linewidth = 0.8) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = scales::hue_pal()(length(levels(results$species)))) +
    labs(
        title = "Load Ratio by Species",
        x     = NULL,
        y     = "Load Ratio"
    ) +
    theme_bw() +
    theme(
        plot.title      = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.text.x     = element_text(angle = 45, hjust = 1, face = "italic", size = 10),
        axis.title.y    = element_text(size = 12),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),
        legend.position = "none"
    )

print(p)
dev.off()

cat("\n=== DONE ===\n")
cat("Output files:\n")
cat("  genetic_load_results_85pct_", site_threshold, ".csv\n", sep="")
cat("  temporal_comparison_summary_", site_threshold, ".csv\n", sep="")
cat("  species_summary_", site_threshold, ".csv\n", sep="")
cat("  genetic_load_", site_threshold, ".png\n", sep="")

RSCRIPT

# Run the R script
echo "Running R analysis..."
Rscript genetic_load_analysis.R "$SITE_THRESHOLD"

if [ $? -eq 0 ]; then
    echo ""
    echo "=== ANALYSIS COMPLETE ==="
    echo "Completed at: $(date)"
    echo ""
    echo "Output files:"
    ls -1 genetic_load_results_85pct_${SITE_THRESHOLD}.csv 2>/dev/null | sed 's/^/  /'
    ls -1 temporal_comparison_summary_${SITE_THRESHOLD}.csv 2>/dev/null | sed 's/^/  /'
    ls -1 species_summary_${SITE_THRESHOLD}.csv 2>/dev/null | sed 's/^/  /'
    ls -1 genetic_load_${SITE_THRESHOLD}.png 2>/dev/null | sed 's/^/  /'
else
    echo "ERROR: R analysis failed"
    exit 1
fi
