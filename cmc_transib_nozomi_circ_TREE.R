library(ape)
library(dplyr)
library(ggtree)
library(ggplot2)

outdir <- "outputs/phylo/cmc_extended_tree"

tr <- read.tree(
  file.path(outdir, "cmc_extended_tree.contree")
)

if (is.rooted(tr)) {
  tr <- unroot(tr)
}

pink_ids <- c(
  "Hztransib",
  "Dmel_transib1",
  "Transib-6_DBipectinata"
)

hopper_id <- "Dmel_hopper"

meta <- read.delim(
  file.path(outdir, "cmc_extended_metadata.tsv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  distinct(tip_label, .keep_all = TRUE) %>%
  transmute(
    label = tip_label,
    sequence_id = sequence_id,
    highlight = case_when(
      group == "RAG1" |
        grepl("RAG[-_. ]?1", sequence_id, ignore.case = TRUE) ~ "RAG1",
      sequence_id == hopper_id ~ "Dmel Hopper",
      sequence_id %in% pink_ids ~ "Selected Transib",
      TRUE ~ "Other"
    )
  )

missing_meta <- setdiff(tr$tip.label, meta$label)

if (length(missing_meta) > 0) {
  stop(
    "Metadata missing for:\n",
    paste(missing_meta, collapse = "\n")
  )
}

missing_highlights <- setdiff(
  c(pink_ids, hopper_id),
  meta$sequence_id
)

if (length(missing_highlights) > 0) {
  warning(
    "Highlighted IDs not found in metadata:\n",
    paste(missing_highlights, collapse = "\n")
  )
}

highlight_cols <- c(
  "Other" = "grey70",
  "RAG1" = "#2166AC",
  "Selected Transib" = "#F28CB1",
  "Dmel Hopper" = "#FF0000"
)

p0 <- ggtree(
  tr,
  layout = "unrooted",
  branch.length = "none",
  color = "grey55",
  linewidth = 0.35
)

p0$data <- p0$data %>%
  left_join(meta, by = "label")

tip_other <- p0$data %>%
  filter(isTip, highlight == "Other")

tip_rag1 <- p0$data %>%
  filter(isTip, highlight == "RAG1")

tip_pink <- p0$data %>%
  filter(isTip, highlight == "Selected Transib")

tip_hopper <- p0$data %>%
  filter(isTip, highlight == "Dmel Hopper")

p_tree <- p0 +
  geom_tiplab(
    data = tip_other,
    aes(label = sequence_id),
    size = 1.25,
    color = "grey20",
    offset = 0.003
  ) +
  geom_tippoint(
    data = tip_rag1,
    shape = 21,
    fill = highlight_cols["RAG1"],
    color = "black",
    stroke = 0.35,
    size = 4.5
  ) +
  geom_tiplab(
    data = tip_rag1,
    aes(label = sequence_id),
    size = 3,
    color = highlight_cols["RAG1"],
    fontface = "bold",
    offset = 0.006
  ) +
  geom_tippoint(
    data = tip_pink,
    shape = 21,
    fill = highlight_cols["Selected Transib"],
    color = "black",
    stroke = 0.35,
    size = 4
  ) +
  geom_tiplab(
    data = tip_pink,
    aes(label = sequence_id),
    size = 2.7,
    color = highlight_cols["Selected Transib"],
    fontface = "bold",
    offset = 0.006
  ) +
  geom_tippoint(
    data = tip_hopper,
    shape = 21,
    fill = highlight_cols["Dmel Hopper"],
    color = "black",
    stroke = 0.4,
    size = 5
  ) +
  geom_tiplab(
    data = tip_hopper,
    aes(label = sequence_id),
    size = 3.1,
    color = highlight_cols["Dmel Hopper"],
    fontface = "bold",
    offset = 0.007
  ) +
  coord_equal(clip = "off") +
  labs(
    title = paste0(
      "CMC (CACTA/Mirage/Chapaev), Transib, and related ",
      "DDE transposase lineages"
    )
  ) +
  theme_void() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    plot.margin = margin(30, 70, 30, 70)
  )

print(p_tree)

ggsave(
  "outputs/plots/cmc_rag1_hopper_tree_unrooted_labels.pdf",
  p_tree,
  width = 12,
  height = 12,
  device = cairo_pdf
)
# 
# ggsave(
#   file.path(outdir, "cmc_extended_tree_unrooted_highlighted.png"),
#   p_tree,
#   width = 16,
#   height = 16,
#   dpi = 400,
#   limitsize = FALSE
# )