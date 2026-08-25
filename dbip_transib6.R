library(DECIPHER)
library(Biostrings)
library(dplyr)
library(ggplot2)
library(ape)
library(Rsamtools)
library(rtracklayer)

dbip.fasta <- "~/Work/projects/DrosoTE/fifteen_genomes/Dbip/Dbip.15g.fasta"
dbip_transib6 <- import.bed("~/Work/projects/DrosoTE/Hopper_A6/dbip_transib6.bed")
dbip_transib6 <-   GenomicRanges::reduce(dbip_transib6, ignore.strand = T, min.gapwidth = 300)
dbip_transib6$id <- paste("t6", 1:length(dbip_transib6), sep = "_")
names(dbip_transib6) <- dbip_transib6$id

genomio <- FaFile(dbip.fasta)
open(genomio)
dbip_t6.fa <- getSeq(genomio, dbip_transib6)
names(dbip_t6.fa) <- dbip_transib6$id

dbt6_wtirend.fa <- dbip_t6.fa[endsWith(as.character(dbip_t6.fa), "CATAGTGC")]
dbt6_wtirstart.fa <- dbip_t6.fa[startsWith(as.character(dbip_t6.fa), "GCACTATG")]
dbt6_bothtirs.fa <- dbip_t6.fa[startsWith(as.character(dbip_t6.fa), "GCACTATGGGGT") & endsWith(as.character(dbip_t6.fa), "CCCCATAGTGC")]

dbt6_tirs.gr <- dbip_transib6[names(dbt6_bothtirs.fa)]
hist(width(dbt6_bothtirs.fa), breaks = 50, main = "Histogram of Transib6 insertions\nwith intact TIRs sizes\nin D bipectinata genome")
hist(width(dbip_t6.fa), breaks = 50, main = "Histogram of Transib6 insertions\nin D bipectinata genome")

dbip_t6.aln <- AlignSeqs(dbip_t6.fa[width(dbip_t6.fa) > 1300])
BrowseSeqs(dbip_t6.aln)

dbip_t6_ends.aln <- AlignSeqs(dbt6_wtirend.fa)
BrowseSeqs(dbip_t6_ends.aln)

dbip_t6_starts.aln <- AlignSeqs(dbt6_wtirstart.fa)
BrowseSeqs(dbip_t6_starts.aln)

dbipt6_tir_names <- unique(c(names(dbt6_wtirstart.fa), names(dbt6_wtirend.fa)))

dbip_t6_tir.aln <- AlignSeqs(dbip_t6.fa[dbipt6_tir_names])
BrowseSeqs(dbip_t6_tir.aln)

hist(width(dbip_transib6), breaks = 50)

subseq(dbip_t6_starts.aln, 1, 60)
reverseComplement(subseq(dbip_t6_ends.aln, 3068 - 59, 3068))

#### Histos for Dbip Transib-6 ####


wdf <- bind_rows(
  data.frame(
    width = width(dbip_t6.fa),
    set = "Dbip Transib-6 all insertions"
  ),
  data.frame(
    width = width(dbt6_bothtirs.fa),
    set = "Dbip Transib-6 insertions with both TIRs"
  )
)

wdf$set <- factor(
  wdf$set,
  levels = c(
    "Dbip Transib-6 all insertions",
    "Dbip Transib-6 insertions with both TIRs"
  )
)

## Common x-axis and binning
xmax <- ceiling(max(wdf$width, na.rm = TRUE) / 500) * 500
xlim_use <- c(0, xmax)

bins_use <- 50
binwidth_use <- diff(xlim_use) / bins_use

## Add one extra boundary so the rightmost value is included
breaks_use <- seq(
  xlim_use[1],
  xlim_use[2] + binwidth_use,
  by = binwidth_use
)

## Pre-bin counts
hist_df <- wdf %>%
  filter(
    width >= xlim_use[1],
    width <= xlim_use[2]
  ) %>%
  mutate(
    bin = cut(
      width,
      breaks = breaks_use,
      include.lowest = TRUE,
      right = FALSE
    ),
    bin_left = breaks_use[as.integer(bin)],
    bin_right = bin_left + binwidth_use
  ) %>%
  count(
    set,
    bin_left,
    bin_right,
    name = "n"
  ) %>%
  filter(n > 0)

## Median insertion length only for insertions with both TIRs
med_df <- wdf %>%
  filter(set == "Dbip Transib-6 insertions with both TIRs") %>%
  group_by(set) %>%
  summarise(
    median_width = median(width),
    .groups = "drop"
  )

transib6_cols <- c(
  "Dbip Transib-6 all insertions" = "#8E44AD",
  "Dbip Transib-6 insertions with both TIRs" = "#8E44AD"
)

p_t6 <- ggplot(hist_df) +
  geom_rect(
    aes(
      xmin = bin_left,
      xmax = bin_right,
      ymin = 0.5,
      ymax = n,
      fill = set
    ),
    color = "grey20",
    linewidth = 0.18,
    alpha = 0.85
  ) +
  geom_vline(
    data = med_df,
    aes(xintercept = median_width),
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey15"
  ) +
  facet_wrap(
    ~set,
    axes = "all_x",
    axis.labels = "all_x",
    nrow = 1,
    ncol = 2
  ) +
  scale_fill_manual(
    values = transib6_cols
  ) +
  scale_x_continuous(
    limits = xlim_use,
    breaks = seq(0, xmax, by = 500),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_log10(
    breaks = c(1, 10, 100),
    labels = c("1", "10", "100"),
    limits = c(
      0.5,
      max(hist_df$n, na.rm = TRUE) * 1.2
    ),
    expand = expansion(mult = c(0, 0.06))
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = "Insertion length (bp)",
    y = "Copy number"
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      size = 12,
      face = "bold"
    ),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.line = element_line(
      color = "black",
      linewidth = 0.35
    ),
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.35
    ),
    axis.ticks.length = unit(2, "pt"),
    panel.grid = element_blank(),
    panel.spacing.x = unit(1, "lines"),
    legend.position = "none",
    plot.margin = margin(6, 8, 6, 8)
  )

p_t6
med_df$median_width[1]

ggsave(
  "outputs/plots/dbip_transib6_width_histograms_counts.pdf",
  p_t6,
  width = 10,
  height = 4,
  useDingbats = FALSE
)

#### Decay hist ####

source("functions/te_decay.R")

dbipt6.all.decay <- get_te_decay_to_ref(
  dbip_t6.fa,
  ref = "t6_1",
  set_name = "Dbip_Transib-6",
  include_ref = TRUE
)

p_decay_dbt6 <- ggplot(
  dbipt6.all.decay,
  aes(x = te_decay, weight = n)
) +
  geom_histogram(
    binwidth = 0.025,
    boundary = 0,
    fill = "#8E44AD",
    color = "grey20",
    alpha = 0.75,
    linewidth = 0.25
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = "TE decay index",
    y = "Copy number per bin"
  ) +
  ggtitle("Transib-6 D. bipectinata")+
  theme(
    axis.text = element_text(color = "black"),
    legend.position = "none"
  )

p_decay_dbt6

pdf("outputs/plots/dbip_transib6_decay_index.pdf", width = 5, height = 8)
p_t6/p_decay_dbt6+patchwork::plot_layout(heights = c(3, 1))
dev.off()

