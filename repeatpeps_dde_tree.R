library(ape)
library(phangorn)
library(dplyr)
library(ggtree)
library(ggplot2)
library(Polychrome)
library(ggforce)

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

tip_group <- setNames(
  meta$superfamily[match(tr$tip.label, meta$label)],
  tr$tip.label
)

ntip <- Ntip(tr)
all_nodes <- sort(unique(c(tr$edge)))

get_descendant_tips <- function(node) {
  if (node <= ntip) {
    return(node)
  }
  
  phangorn::Descendants(
    tr,
    node,
    type = "tips"
  )[[1]]
}

node_group <- vapply(
  all_nodes,
  function(node) {
    descendant_tips <- get_descendant_tips(node)
    
    groups <- unique(
      tip_group[tr$tip.label[descendant_tips]]
    )
    
    groups <- groups[!is.na(groups)]
    
    if (length(groups) == 1) {
      groups
    } else {
      NA_character_
    }
  },
  character(1)
)

parent_lookup <- setNames(
  tr$edge[, 1],
  tr$edge[, 2]
)

pure_clade_roots <- data.frame(
  node = all_nodes,
  superfamily = node_group,
  stringsAsFactors = FALSE
) %>%
  mutate(
    parent = parent_lookup[as.character(node)],
    parent_group = node_group[match(parent, all_nodes)]
  ) %>%
  filter(
    !is.na(superfamily),
    is.na(parent_group) | parent_group != superfamily
  ) %>%
  rowwise() %>%
  mutate(
    n_tips = length(get_descendant_tips(node))
  ) %>%
  ungroup() %>%
  filter(n_tips >= 3) %>%
  mutate(
    cluster_id = paste0(superfamily, "_", node)
  )

p0 <- ggtree(
  tr,
  layout = "unrooted",
  color = NA
)

p0$data <- p0$data %>%
  left_join(meta, by = "label")

hull_df <- bind_rows(
  lapply(seq_len(nrow(pure_clade_roots)), function(i) {
    root_node <- pure_clade_roots$node[i]
    
    descendant_nodes <- c(
      root_node,
      phangorn::Descendants(
        tr,
        root_node,
        type = "all"
      )[[1]]
    )
    
    p0$data %>%
      filter(node %in% descendant_nodes) %>%
      transmute(
        x = x,
        y = y,
        superfamily = pure_clade_roots$superfamily[i],
        cluster_id = pure_clade_roots$cluster_id[i]
      )
  })
)

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
  ggforce::geom_mark_hull(
    data = hull_df,
    inherit.aes = FALSE,
    aes(
      x = x,
      y = y,
      group = cluster_id,
      fill = superfamily,
      color = superfamily
    ),
    concavity = 2,
    expand = unit(2.5, "mm"),
    alpha = 0.10,
    linewidth = 0.35,
    show.legend = FALSE
  ) +
  geom_tree(
    color = "grey60",
    linewidth = 0.35
  ) +
  geom_tippoint(
    data = tip_other,
    aes(color = superfamily),
    size = 2,
    alpha = 0.95
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
  scale_fill_manual(
    values = sf_cols,
    guide = "none"
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

