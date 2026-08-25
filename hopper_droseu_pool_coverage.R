suppressPackageStartupMessages({
  library(tidyverse)
  library(zoo)
})

#### Hopper-Nakanuki coverage from DrosEU pool-seq ####

coverage_dir <- "~/Work/projects/DrosoTE/Hopper_A6/DrosEU_pool/coverage/"

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
  sample <- basename(f) %>%
    sub("\\.depth\\.tsv$", "", x = .)
  
  read_tsv(
    f,
    col_names = c("ref", "pos", "depth", "norm_depth"),
    show_col_types = FALSE
  ) %>%
    mutate(sample = sample)
}) %>%
  group_by(sample, ref) %>%
  mutate(
    max_norm_depth = max(norm_depth, na.rm = TRUE),
    norm_depth_max = ifelse(
      max_norm_depth > 0,
      norm_depth / max_norm_depth,
      NA_real_
    )
  ) %>%
  ungroup()

depth_df <- depth_df |>
  mutate(
    sample = factor(sample),
    ref = "A6 Nozomi",
    ref = factor(ref)
  )

# Remove signal from polyT
# depth_df$depth[depth_df$pos > 142 & depth_df$pos < 179] <- 0
# depth_df$norm_depth[depth_df$pos > 142 & depth_df$pos < 179] <- 0

#### basic summary ####

summary_df <- depth_df |>
  group_by(sample, ref) |>
  summarise(
    ref_length = max(pos),
    
    mean_depth = mean(norm_depth, na.rm = TRUE),
    median_depth = median(norm_depth, na.rm = TRUE),
    max_depth = max(norm_depth, na.rm = TRUE),
    covered_bp = sum(norm_depth > 0, na.rm = TRUE),
    frac_covered = covered_bp / ref_length,
    
    central_1000_depth = sum(norm_depth[pos >= 1000 & pos <= 2000], na.rm = TRUE),
    
    flank_depth = sum(
      norm_depth[pos <= 500 | pos > (ref_length - 500)],
      na.rm = TRUE
    ),
    
    middle_loss_index = ifelse(
      flank_depth > 0,
      central_1000_depth / flank_depth,
      NA_real_
    ),
    
    .groups = "drop"
  )

q90 <- quantile(summary_df$middle_loss_index, 0.90, na.rm = TRUE)

high_middle <- summary_df |>
  filter(!is.na(middle_loss_index)) |>
  filter(middle_loss_index >= q90) |>
  pull(sample) %>% as.character()


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


p_norm <- ggplot(depth_df%>% filter(sample %in% high_middle), aes(x = pos, y = norm_depth)) +
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

p_overlay_high <- ggplot(
  depth_smooth %>% filter(sample %in% high_middle),
  aes(x = pos, y = depth_smooth, group = sample, color = sample)
) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = 2) +
  
  # regions used for the middle-loss index
  geom_vline(xintercept = c(500, 1000, 2000, 2347), linewidth = 0.3, linetype = 3) +
  
  geom_line(linewidth = 0.35, alpha = 0.75, na.rm = TRUE) +
  facet_wrap(~ ref, scales = "free_x") +
  theme_classic() +
  guides(color = "none") +
  labs(
    x = "Position along reference",
    y = "Depth (smooth window 50 bp)",
    title = "Dmel pool-seq samples with highest middle part coverage overlaid"
  )

ref_len <- max(coverage_df$pos)
# x_marks <- c(1, 438, 1819, ref_len)
x_marks <- c(1, 1040, 2419, ref_len)

p_overlay <- ggplot(
  depth_smooth,
  aes(x = rev(pos), y = depth_smooth, group = sample)
) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = 2) +
  geom_line(linewidth = 0.35, alpha = 0.25, na.rm = TRUE) +
  scale_x_continuous(
    breaks = x_marks,
    labels = x_marks,
    limits = c(1, ref_len),
    expand = expansion(mult = c(0.1, 0.1))
  )+
  facet_wrap(~ ref, scales = "free_x") +
  theme_classic() +
  labs(
    x = "Position along reference",
    y = "Depth (smooth window 50 bp)",
    title = "Dmel pool-seq all samples overlaid"
  )

p_overlay_norm <- ggplot(
  depth_df,
  aes(x = rev(pos), y = norm_depth_max, group = sample)
) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = 2) +
  geom_line(linewidth = 0.35, alpha = 0.15, na.rm = TRUE) +
  facet_wrap(~ ref, scales = "free_x") +
  scale_x_continuous(
    breaks = x_marks,
    labels = x_marks,
    limits = c(1, ref_len),
    expand = expansion(mult = c(0.1, 0.1))
  )+
  theme_classic() +
  labs(
    x = "Position along reference",
    y = "Depth normalized on max value",
    title = "DrosEU all samples from 48 populations overlaid"
  )

pdf("outputs/plots/nozomi_droseu_poolseq_cov.pdf", width = 6, height = 4)
p_overlay_norm
dev.off()
