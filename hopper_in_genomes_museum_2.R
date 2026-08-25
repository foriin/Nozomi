
suppressPackageStartupMessages({
  library(tidyverse)
  library(zoo)
})

#### inputs ####

coverage_dir <- "~/Work/projects/DrosoTE/Hopper_A6/museum/coverage"

#### read depth files ####
#### expected samtools depth format:
#### ref_name    position    coverage

depth_files <- list.files(
  coverage_dir,
  pattern = "\\.norm\\.depth\\.tsv$",
  full.names = TRUE
)

if (length(depth_files) == 0) {
  stop("No *.depth.tsv files found in: ", coverage_dir)
}

depth_df <- map_dfr(depth_files, function(f) {
  sample <- basename(f) |>
    sub("\\.depth\\.tsv$", "", x = _)
  
  read_tsv(
    f,
    col_names = c("ref", "pos", "depth", "norm_depth"),
    show_col_types = FALSE
  ) |>
    mutate(sample = sample)
})

depth_df <- depth_df |>
  mutate(
    sample = factor(sample),
    ref = factor(ref)
  )

# Remove signal from polyT
depth_df$depth[depth_df$pos > 142 & depth_df$pos < 179] <- 0
depth_df$norm_depth[depth_df$pos > 142 & depth_df$pos < 179] <- 0

#### basic summary ####

summary_df <- depth_df |>
  group_by(sample, ref) |>
  summarise(
    ref_length = max(pos),
    mean_depth = mean(norm_depth),
    median_depth = median(norm_depth),
    max_depth = max(norm_depth),
    covered_bp = sum(norm_depth > 0),
    frac_covered = covered_bp / ref_length,
    .groups = "drop"
  )

# write_tsv(summary_df, file.path(outdir, "depth_summary.tsv"))

#### optional smoothing ####

window_size <- 50

depth_smooth <- depth_df |>
  arrange(sample, ref, pos) |>
  group_by(sample, ref) |>
  mutate(
    depth_smooth = zoo::rollmean(norm_depth, k = window_size, fill = NA, align = "center")
  ) |>
  ungroup()

#### plot 1: raw depth, one PDF page with facets ####

p_raw <- ggplot(depth_df, aes(x = pos, y = depth)) +
  geom_line(linewidth = 0.25) +
  facet_wrap(~ sample, ncol = 5, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Position along reference",
    y = "Read depth",
    title = "Base-pair coverage over Hopper/reference"
  )

ggsave(
  filename = file.path(outdir, "depth_raw_facets.pdf"),
  plot = p_raw,
  width = 16,
  height = 12,
  limitsize = FALSE
)

#### plot 2: smoothed depth, using normalized counts, one PDF page with facets ####

p_smooth <- ggplot(depth_smooth, aes(x = pos, y = depth_smooth)) +
  geom_line(linewidth = 0.35, na.rm = TRUE) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = 2)+
  facet_wrap(~ sample, ncol = 5, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Position along reference",
    y = paste0("Smoothed read depth, ", window_size, " bp window"),
    title = "Smoothed coverage over Hopper/reference"
  )

ggsave(
  filename = file.path(outdir, "depth_smoothed_facets.pdf"),
  plot = p_smooth,
  width = 16,
  height = 12,
  limitsize = FALSE
)

#### plot 3: normalized depth per sample ####
#### useful if libraries have very different sequencing depth ####


p_norm <- ggplot(depth_df, aes(x = pos, y = norm_depth)) +
  geom_line(linewidth = 0.25) +
  facet_wrap(~ sample, ncol = 5, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Position along reference",
    y = "Depth per million dm6 mapped bases",
    title = "Normalized base-pair coverage over Hopper/reference"
  )

ggsave(
  filename = file.path(outdir, "depth_normalized_facets.pdf"),
  plot = p_norm,
  width = 16,
  height = 12,
  limitsize = FALSE
)

#### plot 4: all samples overlaid, normalized ####

p_overlay <- ggplot(depth_smooth, aes(x = pos, y = depth_smooth, group = sample)) +
  geom_line(linewidth = 0.25, alpha = 0.5) +
  facet_wrap(~ ref, scales = "free_x") +
  theme_classic() +
  labs(
    x = "Position along reference",
    y = "Depth (smooth window 50 bp)",
    title = "Museum samples overlaid"
  )

ggsave(
  filename = file.path(outdir, "depth_normalized_overlay.pdf"),
  plot = p_overlay,
  width = 12,
  height = 5,
  limitsize = FALSE
)

#### plot 5: heatmap along reference ####

bin_size <- 25

depth_binned <- depth_df |>
  mutate(bin = ceiling(pos / bin_size)) |>
  group_by(sample, ref, bin) |>
  summarise(
    start = min(pos),
    end = max(pos),
    mid = mean(c(start, end)),
    mean_depth = mean(depth),
    .groups = "drop"
  ) |>
  group_by(sample) |>
  mutate(
    mean_depth_norm = mean_depth / sum(mean_depth) * 1e6
  ) |>
  ungroup()

p_heat <- ggplot(depth_binned, aes(x = mid, y = sample, fill = mean_depth_norm)) +
  geom_tile() +
  facet_wrap(~ ref, scales = "free_x") +
  theme_classic() +
  labs(
    x = "Position along reference",
    y = NULL,
    fill = "Norm. depth",
    title = paste0("Coverage heatmap, ", bin_size, " bp bins")
  )

ggsave(
  filename = file.path(outdir, "depth_heatmap_binned.pdf"),
  plot = p_heat,
  width = 14,
  height = 8,
  limitsize = FALSE
)

#### save processed tables ####

write_tsv(depth_df, file.path(outdir, "all_samples.depth.long.tsv"))
write_tsv(depth_norm, file.path(outdir, "all_samples.depth.normalized.long.tsv"))
write_tsv(depth_binned, file.path(outdir, "all_samples.depth.binned.tsv"))

message("Done. Plots written to: ", outdir)