suppressPackageStartupMessages({
  library(tidyverse)
  library(zoo)
})

#### Dsuzuki Hopper TE (protein MSA similarity = 73%) ####

#### 1. DNA-seq coverage across > 200 sequencing datasets from all continet pops ####

coverage_dir <- "~/Work/projects/DrosoTE/Hopper_A6/Dsuz_hopp_cov/coverage"

###### read depth files ######
# expected samtools depth format:
# ref_name    position    coverage_raw    coverage_norm

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
    ref = factor(ref)
  )

###### basic summary ######

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

depth_df <- depth_df %>% filter(sample %in% summary_df$sample[summary_df$frac_covered > 0.95])

q10 <- quantile(summary_df$middle_loss_index, 0.10, na.rm = TRUE)

low_middle <- summary_df |>
  filter(!is.na(middle_loss_index)) |>
  filter(middle_loss_index <= q10) |>
  pull(sample) %>% as.character()


# write_tsv(summary_df, file.path(outdir, "depth_summary.tsv"))

###### optional smoothing ######

window_size <- 50

depth_smooth <- depth_df |>
  arrange(sample, ref, pos) |>
  group_by(sample, ref) |>
  mutate(
    depth_smooth = zoo::rollmean(norm_depth, k = window_size, fill = NA, align = "center")
  ) |>
  ungroup()

###### plot 1: raw depth, one PDF page with facets ######

p_raw <- ggplot(depth_df, aes(x = pos, y = depth)) +
  geom_line(linewidth = 0.25) +
  facet_wrap(~ sample, ncol = 5, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Position along reference",
    y = "Read depth",
    title = "Base-pair coverage over Hopper/reference"
  )



###### plot 2: smoothed depth, using normalized counts, one PDF page with facets ######

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



###### plot 3: normalized depth per sample ######
# useful if libraries have very different sequencing depth


p_norm <- ggplot(depth_df%>% filter(sample %in% low_middle), aes(x = pos, y = norm_depth)) +
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

###### plot 4: all samples overlaid, normalized ######

p_overlay_low <- ggplot(
  depth_smooth %>% filter(sample %in% low_middle),
  aes(x = pos, y = depth_smooth, group = sample)
) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = 2) +
  
  # regions used for the middle-loss index
  geom_vline(xintercept = c(500, 1000, 2000, 2500), linewidth = 0.3, linetype = 3) +
  
  geom_line(linewidth = 0.35, alpha = 0.75, na.rm = TRUE) +
  facet_wrap(~ ref, scales = "free_x") +
  theme_classic() +
  guides(color = "none") +
  labs(
    x = "Position along reference",
    y = "Depth (smooth window 50 bp)",
    title = "Dsuzukii samples with lowest middle part coverage overlaid"
  )

p_overlay <- ggplot(
  depth_smooth,
  aes(x = pos, y = depth_smooth, group = sample)
) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = 2) +
  geom_line(linewidth = 0.35, alpha = 0.65, na.rm = TRUE) +
  facet_wrap(~ ref, scales = "free_x") +
  theme_classic() +
  labs(
    x = "Position along reference",
    y = "Depth (smooth window 50 bp)",
    title = "Dsuzukii all samples overlaid"
  )

depth_df <- depth_df %>% filter(!grepl("WA[1-5]", sample)) %>% 
  mutate(ref = "Hayabusa")

p_overlay_norm <- ggplot(
  depth_df,
  aes(x = pos, y = norm_depth_max, group = sample)
) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = 2) +
  geom_line(linewidth = 0.35, alpha = 0.25, na.rm = TRUE) +
  facet_wrap(~ ref, scales = "free_x") +
  theme_classic() +
  labs(
    x = "Position along reference",
    y = "Depth normalized on max value",
    title = "Dsuzukii all samples overlaid"
  )

ggsave("outputs/plots/dsuz_hayabusa_cov.pdf",
       p_overlay_norm,
       width = 6,
       height = 4)

depth_df_dsub <- depth_df %>% filter(grepl("WA[1-5]", sample))

p_overlay_dsub <- ggplot(
  depth_df_dsub,
  aes(x = pos, y = depth, group = sample)
) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = 2) +
  geom_line(linewidth = 0.35, alpha = 0.65, na.rm = TRUE) +
  facet_wrap(~ ref, scales = "free_x") +
  theme_classic() +
  labs(
    x = "Position along reference",
    y = "Depth (raw reads)",
    title = "Dsub pops from 2017 samples overlaid"
  )

#### 2. DNA sequences of the insertions ####

dsuz_hop_true.bed <- import.bed("/Users/artemilin/Work/projects/DrosoTE/Hopper_A6/dsuz_hopper_true_genome.bed")
seqlevels(dsuz_hop_true.bed) <- paste0(seqlevels(dsuz_hop_true.bed), ".")
dsuz_hop_true.bed <- resize(dsuz_hop_true.bed, width = width(dsuz_hop_true.bed) + 20, 
                            fix = 'center')

hist(width(dsuz_hop_true.bed))
dsuz_genome <- "/Users/artemilin/Work/projects/DrosoTE/fifteen_genomes/Dsuz/Dsuz_2.1_genome.fasta"

dsuzfa <- FaFile(dsuz_genome)
open(dsuzfa)

dsuz_ht_fa <- getSeq(dsuzfa, dsuz_hop_true.bed)
dsuz_ht_aln <- AlignSeqs(dsuz_ht_fa)
BrowseSeqs(dsuz_ht_aln)

# take first 30 and last 30 alignment columns
dsuz_aln_ends <- DNAStringSet(
  paste0(
    subseq(dsuz_ht_aln, start = 1, width = 30),
    subseq(dsuz_ht_aln, start = width(dsuz_ht_aln) - 29, width = 30)
  )
)

names(dsuz_aln_ends) <- names(dsuz_ht_aln)

BrowseSeqs(dsuz_aln_ends)

subseq(dsuz_ht_fa[1], 8, 44)


#### Another hopper variant in Dsuz - hopster (protein MSA similarity = 82%) ####
hopster <- readDNAStringSet("~/Work/projects/DrosoTE/Hopper_A6/dsuz_hopster.fa")
hopster <- subseq(hopster, 3, 2863)
subseq(hopster, 1, 42)
reverseComplement(subseq(hopster, width(hopster) - 45, width(hopster)))
writeXStringSet(hopster, "~/Work/projects/DrosoTE/Hopper_A6/dsuz_hopster.fa")

#### 1. DNA-seq coverage across > 200 sequencing datasets from all continet pops ####

coverage_dir <- "~/Work/projects/DrosoTE/Hopper_A6/Dsuz_hopp_cov/hopster_coverage/"

###### read depth files ######
# expected samtools depth format:
# ref_name    position    coverage_raw    coverage_norm

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
    sub("\\.hopper\\.R1\\.norm\\.depth\\.tsv$", "", x = .)
  
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
    ref = factor(ref)
  )

###### basic summary ######

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

depth_df <- depth_df %>% filter(sample %in% summary_df$sample[summary_df$frac_covered > 0.95])

q90 <- quantile(summary_df$middle_loss_index, 0.90, na.rm = TRUE)

high_middle <- summary_df |>
  filter(!is.na(middle_loss_index)) |>
  filter(middle_loss_index >= q90) |>
  pull(sample) %>% as.character()


# write_tsv(summary_df, file.path(outdir, "depth_summary.tsv"))

###### optional smoothing ######

window_size <- 50

depth_smooth <- depth_df |>
  arrange(sample, ref, pos) |>
  group_by(sample, ref) |>
  mutate(
    depth_smooth = zoo::rollmean(norm_depth, k = window_size, fill = NA, align = "center")
  ) |>
  ungroup()

###### plot 1: raw depth, one PDF page with facets ######

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

###### plot 2: smoothed depth, using normalized counts, one PDF page with facets ######

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

###### plot 3: normalized depth per sample ######
# useful if libraries have very different sequencing depth


p_norm <- ggplot(depth_df%>% filter(sample %in% low_middle), aes(x = pos, y = norm_depth)) +
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

###### plot 4: all samples overlaid, normalized ######

p_overlay_high <- ggplot(
  depth_smooth %>% filter(sample %in% high_middle),
  aes(x = pos, y = depth_smooth, group = sample)
) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = 2) +
  
  # regions used for the middle-loss index
  geom_vline(xintercept = c(500, 1000, 2000, 2500), linewidth = 0.3, linetype = 3) +
  
  geom_line(linewidth = 0.35, alpha = 0.75, na.rm = TRUE) +
  facet_wrap(~ ref, scales = "free_x") +
  theme_classic() +
  guides(color = "none") +
  labs(
    x = "Position along reference",
    y = "Depth (smooth window 50 bp)",
    title = "Dsuzukii samples with lowest middle part coverage overlaid"
  )

p_overlay <- ggplot(
  depth_smooth,
  aes(x = pos, y = depth_smooth, group = sample)
) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = 2) +
  geom_line(linewidth = 0.35, alpha = 0.65, na.rm = TRUE) +
  facet_wrap(~ ref, scales = "free_x") +
  theme_classic() +
  labs(
    x = "Position along reference",
    y = "Depth (smooth window 50 bp)",
    title = "Dsuzukii all samples overlaid"
  )

p_overlay_norm <- ggplot(
  depth_df,
  aes(x = pos, y = norm_depth_max, group = sample)
) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = 2) +
  geom_line(linewidth = 0.35, alpha = 0.25, na.rm = TRUE) +
  facet_wrap(~ ref, scales = "free_x") +
  theme_classic() +
  labs(
    x = "Position along reference",
    y = "Depth normalized on max value",
    title = "Dsuzukii all samples overlaid"
  )

#### 2. Hopster sequences and TSD ####
dsuz_hopster.bed <- import.bed("/Users/artemilin/Work/projects/DrosoTE/Hopper_A6/dsuz_hopster_genome.bed")
seqlevels(dsuz_hopster.bed) <- paste0(seqlevels(dsuz_hopster.bed), ".")
dsuz_hopster.bed <- dsuz_hopster.bed[!grepl("NC", seqnames(dsuz_hopster.bed))]
dsuz_hopster.bed <- GenomicRanges::reduce(dsuz_hopster.bed, min.gapwidth = 300)
dsuz_hopster.bed <- resize(dsuz_hopster.bed, width = width(dsuz_hopster.bed) + 20, 
                            fix = 'center')

dsuz_genome <- "/Users/artemilin/Work/projects/DrosoTE/fifteen_genomes/Dsuz/Dsuz_2.1_genome.fasta"

dsuzfa <- FaFile(dsuz_genome)
open(dsuzfa)

dsuz_hopst_fa <- getSeq(dsuzfa, dsuz_hopster.bed)
dsuz_hopst <- readDNAStringSet("/Users/artemilin/Work/projects/DrosoTE/Hopper_A6/DNA/dsuz_hopster.fa")
dsuz_hopst_fa <- c(dsuz_hopst_fa, dsuz_hopst)

dsuz_hopst_fa <- getSeq(dsuzfa, dsuz_hopster.bed[width(dsuz_hopster.bed) > 1500])
dsuz_hopst_aln <- AlignSeqs(dsuz_hopst_fa)
BrowseSeqs(dsuz_hopst_aln)
hist(width(dsuz_hopster.bed))
