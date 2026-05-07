# =============================================================================
# Environment Setup and Argument Parsing
# =============================================================================

# environment
args = commandArgs(trailingOnly = TRUE)

workdir = args[1] # OUTDIR
strain = args[2]
prefix = args[3]
rna_type = args[4]

length_file = paste0(workdir, "/coverage/", strain, "_", rna_type, ".length")
probe_file = paste0(workdir, "/coverage/", rna_type, ".probe.tsv")
depth_file = paste0(workdir, "/", prefix, "/1_rrna/", strain, ".", rna_type, ".depth")

#setwd(workdir)

# =============================================================================
# Package Imports
# =============================================================================

# packages
library(tidyverse)
library(purrr)

# =============================================================================
# Data Loading
# =============================================================================

# Read input data
rna_length = read_tsv(length_file, col_names = c("rna_full", "rna_length")) %>%
    mutate(rna_id = str_split(rna_full, " ", simplify = T)[,1])

probe_region = read_tsv(probe_file, col_names = c("rna_id", "probe_region")) %>%
    separate(probe_region, "-", into = c("probe_start", "probe_end"), convert = TRUE) %>%
    arrange(rna_id, probe_start)

# Merge contiguous probe regions
merged_regions = probe_region %>%
  group_by(rna_id) %>%
  mutate(
    is_continuous = probe_start == lag(probe_end, default = first(probe_start)) + 1,
    merge_group = cumsum(!is_continuous)
  ) %>%
  group_by(rna_id, merge_group) %>%
  summarize(
    start = min(probe_start),
    end = max(probe_end),
    .groups = "drop"
  )


rna_depth = read_tsv(depth_file, col_names = c("rna_id", "rna_position","seq_depth"))

# =============================================================================
# Visualization Function
# =============================================================================

#' Plot Sequencing Coverage for a Given RNA
#'
#' Generates a bar plot of sequencing depth across an RNA sequence, highlighting
#' probe-targeted regions in red and non-probe regions in black.
#'
#' @param rid Character. RNA ID to plot.
#'
#' @return None. Saves a PDF plot to the output directory.
PlotCov = function(rid){
    # Get rna length
    rna_len = rna_length %>% filter(rna_id == rid) %>% pull(rna_length)

    # Get probe regions
    region = merged_regions %>%
        filter(rna_id == rid) %>%
        mutate(region_expr = map2(start, end, ~seq(.x, .y))) %>%
        pull(region_expr) %>%
        unlist()

    depth = rna_depth %>% filter(rna_id == rid) %>%
        mutate(col_group = ifelse(rna_position %in% region, "Probe", "Non-probe"))

    p = ggplot(depth, aes(x = rna_position, y = seq_depth, fill = col_group)) +
        geom_col() +
        scale_x_continuous(limits = c(0, rna_len), name = paste(rid, "of", prefix)) +
        scale_y_continuous(name = "Sequencing depth") +
        scale_fill_manual(values = c("Non-probe" = "black", "Probe" = "red")) +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5),
                text = element_text(size = 10),
                axis.text.x = element_text(hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
                axis.text.y = element_text(size = 8,color = "black"),
                legend.text = element_text(size = 6,color = "black"),
                legend.title = element_blank(),
                legend.key.size = unit(11, "pt"),
                legend.position = "top")

    # Save plot
    ggsave(
        filename = paste0(workdir, "/", prefix, "/1_rrna/", rid, ".pdf"),
        plot = p,
        units = "mm", width = 120, height = 70
    )
}


# =============================================================================
# Execution
# =============================================================================

# Apply the function to all unique RNA IDs
selected_rna = unique(rna_depth$rna_id)
walk(selected_rna, PlotCov)
