library(ape)
library(dplyr)
library(ggtree)
library(ggplot2)

outdir <- "outputs/phylo/cmc_rag1_tree"

tr <- read.tree(
  file.path(outdir, "cmc_rag1_tree.contree")
)

if (is.rooted(tr)) {
  tr <- unroot(tr)
}

meta <- read.delim(
  file.path(outdir, "cmc_rag1_metadata.tsv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  distinct(tip_label, .keep_all = TRUE) %>%
  transmute(
    label = tip_label,
    sequence_id = sequence_id,
    is_rag1 = group == "RAG1"
  )

p0 <- ggtree(
  tr,
  layout = "unrooted",
  branch.length = "none",
  color = "grey60",
  linewidth = 0.35
)

p0$data <- p0$data %>%
  left_join(meta, by = "label")

tip_other <- p0$data %>%
  filter(isTip, !is_rag1)

tip_rag1 <- p0$data %>%
  filter(isTip, is_rag1)

p_tree <- p0 +
  geom_tippoint(
    data = tip_other,
    shape = 21,
    fill = "grey70",
    color = "black",
    stroke = 0.2,
    size = 1.5,
    alpha = 0.65
  ) +
  geom_tippoint(
    data = tip_rag1,
    shape = 21,
    fill = "#2166AC",
    color = "black",
    stroke = 0.35,
    size = 5
  ) +
  geom_tiplab(
    data = tip_other,
    aes(label = sequence_id),
    size = 1.35,
    color = "grey20",
    offset = 0.003
  ) +
  geom_tiplab(
    data = tip_rag1,
    aes(label = sequence_id),
    size = 2.8,
    color = "#2166AC",
    fontface = "bold",
    offset = 0.006
  ) +
  coord_equal(clip = "off") +
  labs(
    title = "CMC (CACTA/Mirage/Chapaev) and Transib DDE Transposon Lineages"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    plot.margin = margin(20, 20, 20, 20)
  )

print(p_tree)

ggsave(
  "outputs/plots/cmc_rag1_tree_unrooted_labels.pdf",
  p_tree,
  width = 12,
  height = 12,
  device = cairo_pdf
)



# ggsave(
#   file.path(outdir, "cmc_rag1_tree_unrooted_labels.png"),
#   p_tree,
#   width = 12,
#   height = 12,
#   dpi = 300
# )