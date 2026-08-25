library(Biostrings)
library(ggplot2)
library(dplyr)

#### Histos for Hopper ####

wdf <- bind_rows(
  data.frame(
    width = width(a6hop_bothtirs),
    set = "A6"
  ),
  data.frame(
    width = width(dm6hop_bothtirs),
    set = "dm6"
  ),
  data.frame(
    width = width(barhop_bothtirs),
    set = "BL2969"
  )
)

wdf$set <- factor(wdf$set, levels = c("dm6", "A6", "BL2969"))

## Common binning
xlim_use <- c(0, 3000)
bins_use <- 50
binwidth_use <- diff(xlim_use) / bins_use

breaks_use <- seq(xlim_use[1], xlim_use[2], by = binwidth_use)

## Pre-bin manually so log-scale counts behave cleanly
hist_df <- wdf %>%
  filter(width >= xlim_use[1], width <= xlim_use[2]) %>%
  mutate(
    bin = cut(
      width,
      breaks = breaks_use,
      include.lowest = TRUE,
      right = FALSE
    ),
    bin_left = breaks_use[as.integer(bin)],
    bin_right = bin_left + binwidth_use,
    bin_mid = bin_left + binwidth_use / 2
  ) %>%
  count(set, bin_left, bin_right, bin_mid, name = "n") %>%
  filter(n > 0)

med_df <- wdf %>%
  group_by(set) %>%
  summarise(
    median_width = median(width),
    .groups = "drop"
  )

hist_cols <- c(
  "dm6"    = "#2E7D32", 
  "A6"     = "#4E79A7",
  "BL2969" = "#D95F5F"
)

p <- ggplot(hist_df) +
  geom_rect(
    aes(
      xmin = bin_left,
      xmax = bin_right,
      ymin = 1,
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
  facet_wrap(~set, nrow = 3) +
  scale_fill_manual(values = hist_cols) +
  scale_x_continuous(
    limits = xlim_use,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_log10(
    breaks = c(1, 10, 100),
    labels = c("1", "10", "100"),
    limits = c(1, max(hist_df$n, na.rm = TRUE) * 1.2),
    expand = expansion(mult = c(0, 0.06))
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = "Insertion length (bp)",
    y = "Copy number"
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.35),
    axis.ticks = element_line(color = "black", linewidth = 0.35),
    axis.ticks.length = unit(2, "pt"),
    
    panel.grid = element_blank(),
    
    ## more space between dm6 / A6 / BL2969 panels
    panel.spacing.x = unit(1.2, "lines"),
    
    legend.position = "none",
    
    ## extra outer whitespace, helps PDF not look cramped
    plot.margin = margin(6, 8, 6, 8)
  )

p

ggsave(
  "outputs/plots/hopper_width_histograms_log_counts_vertical.pdf",
  p,
  width = 4,
  height = 4,
  useDingbats = FALSE
)

