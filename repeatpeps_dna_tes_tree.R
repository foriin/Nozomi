library(ape)
library(dplyr)
library(ggtree)
library(ggplot2)
library(Polychrome)

outdir <- "outputs/phylo/dde_tree2"

tr <- read.tree(
  file.path(outdir, "dde_superfamily_tree.contree")
)

if (is.rooted(tr)) {
  tr <- unroot(tr)
}

meta <- read.delim(
  file.path(outdir, "dde_superfamily_representatives.tsv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  distinct(tip_label, .keep_all = TRUE) %>%
  transmute(
    label = tip_label,
    superfamily = superfamily
  )

missing_meta <- setdiff(tr$tip.label, meta$label)

if (length(missing_meta) > 0) {
  stop(
    "Missing metadata for:\n",
    paste(missing_meta, collapse = "\n")
  )
}

superfamilies <- sort(unique(meta$superfamily))

sf_cols <- setNames(
  Polychrome::glasbey.colors(length(superfamilies)),
  superfamilies
)

sf_cols["CMC-Transib"] <- "#8B0000"

p0 <- ggtree(
  tr,
  layout = "unrooted",
  color = "grey70",
  linewidth = 0.35
)

p0$data <- p0$data %>%
  left_join(meta, by = "label")

tip_other <- p0$data %>%
  filter(
    isTip,
    superfamily != "CMC-Transib"
  )

tip_transib <- p0$data %>%
  filter(
    isTip,
    superfamily == "CMC-Transib"
  )

p_tree <- p0 +
  geom_tippoint(
    data = tip_other,
    aes(color = superfamily),
    size = 2,
    alpha = 0.9
  ) +
  geom_tippoint(
    data = tip_transib,
    shape = 21,
    fill = "#8B0000",
    color = "black",
    stroke = 0.5,
    size = 4
  ) +
  scale_color_manual(
    values = sf_cols,
    breaks = superfamilies,
    name = NULL
  ) +
  guides(
    color = guide_legend(
      ncol = 4,
      byrow = TRUE,
      override.aes = list(
        size = 4,
        alpha = 1
      )
    )
  ) +
  coord_equal(clip = "off") +
  labs(
    title = "Phylogeny of DNA transposase superfamilies"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.key.width = unit(10, "pt"),
    legend.spacing.x = unit(3, "pt"),
    plot.margin = margin(15, 15, 15, 15)
  )

print(p_tree)

# ggsave(
#   file.path(outdir, "dde_superfamily_tree_unrooted.pdf"),
#   p_tree,
#   width = 11,
#   height = 10,
#   device = cairo_pdf
# )
# 
# ggsave(
#   file.path(outdir, "dde_superfamily_tree_unrooted.png"),
#   p_tree,
#   width = 11,
#   height = 10,
#   dpi = 300
# )