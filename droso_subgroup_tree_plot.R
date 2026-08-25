library(ape)
library(phangorn)
library(dplyr)
library(tibble)
library(ggtree)
library(ggplot2)

host_tree <- read.tree("inputs/trees/full_astral_fullAnnotation.tree")

species_groups <- list(
  melanogaster = c(
    "D_MELANOGASTER",
    "D_MAURITIANA",
    "D_YAKUBA",
    "D_SIMULANS",
    "D_ERECTA",
    "D_SECHELLIA"
  ),
  suzukii = c(
    "D_SUZUKII",
    "D_MIMETICA",
    "D_LUCIPENNIS",
    "D_BIARMIPES",
    "D_PULCHRELLA",
    "D_SUBPULCHRELLA",
    "D_OSHIMAI"
  ),
  montium = c(
    "D_KIKKAWAI",
    "D_AURARIA",
    "D_TRIAURARIA",
    "D_RUFA",
    "D_TRAPEZIFRONS",
    "D_WATANABEI"
  ),
  ananassae = c(
    "D_ANANASSAE",
    "D_BIPECTINATA",
    "D_VARIANS",
    "D_ANOMALATA",
    "D_IRONENSIS"
  )
)

focal_species <- c(
  "D_MELANOGASTER",
  "D_MAURITIANA",
  "D_SUZUKII",
  "D_SUBPULCHRELLA",
  "D_OSHIMAI",
  "D_AURARIA",
  "D_TRIAURARIA",
  "D_RUFA",
  "D_BIPECTINATA"
)

wanted_species <- unique(unlist(species_groups))

missing_species <- setdiff(wanted_species, host_tree$tip.label)
missing_species

host_tree_pruned <- keep.tip(
  host_tree,
  intersect(wanted_species, host_tree$tip.label)
)

host_tree_pruned <- ladderize(host_tree_pruned)

host_meta <- bind_rows(
  lapply(names(species_groups), function(group_name) {
    data.frame(
      label = species_groups[[group_name]],
      group = group_name,
      stringsAsFactors = FALSE
    )
  })
) %>%
  filter(label %in% host_tree_pruned$tip.label) %>%
  mutate(
    focal = label %in% focal_species,
    species_label = paste0(
      "D. ",
      tolower(sub("^D_", "", label))
    ),
    group = factor(
      group,
      levels = c(
        "melanogaster",
        "suzukii",
        "montium",
        "ananassae"
      )
    )
  )

group_cols <- c(
  melanogaster = "#4C78A8",
  suzukii = "#F58518",
  montium = "#54A24B",
  ananassae = "#B279A2"
)

p0 <- ggtree(
  host_tree_pruned,
  layout = "rectangular",
  branch.length = "none"
)

p0$data <- p0$data %>%
  left_join(host_meta, by = "label")

tip_nonfocal <- p0$data %>%
  filter(isTip, !focal)

tip_focal <- p0$data %>%
  filter(isTip, focal)

p_host <- p0 +
  geom_tippoint(
    data = tip_nonfocal,
    aes(fill = group),
    shape = 21,
    size = 2.2,
    color = "black",
    stroke = 0.25
  ) +
  geom_tippoint(
    data = tip_focal,
    aes(fill = group),
    shape = 21,
    size = 4.2,
    color = "black",
    stroke = 0.8
  ) +
  geom_tiplab(
    data = tip_nonfocal,
    aes(
      label = species_label,
      color = group
    ),
    fontface = "italic",
    size = 3,
    offset = 0.15
  ) +
  geom_tiplab(
    data = tip_focal,
    aes(
      label = species_label,
      color = group
    ),
    fontface = "bold.italic",
    size = 3.4,
    offset = 0.15
  ) +
  scale_fill_manual(
    values = group_cols,
    name = "Subgroup"
  ) +
  scale_color_manual(
    values = group_cols,
    name = "Subgroup"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.02, 0.3))
  ) +
  guides(
    color = "none",
    fill = guide_legend(
      override.aes = list(
        shape = 21,
        size = 4,
        color = "black"
      )
    )
  ) +
  theme_tree() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

pdf("outputs/plots/droso_spec_subgroups.pdf", width = 6.5, height = 5)
print(p_host)
dev.off()


#### Add B dorsalis ####
host_tree_pruned <- keep.tip(
  host_tree,
  intersect(wanted_species, host_tree$tip.label)
)

drosophila_newick <- sub(
  ";$",
  "",
  write.tree(host_tree_pruned)
)

host_tree_pruned <- read.tree(
  text = paste0(
    "(",
    drosophila_newick,
    ",B_DORSALIS);"
  )
)

host_tree_pruned$edge.length <- NULL
host_tree_pruned <- ladderize(host_tree_pruned)


host_meta <- bind_rows(
  lapply(names(species_groups), function(group_name) {
    data.frame(
      label = species_groups[[group_name]],
      group = group_name,
      stringsAsFactors = FALSE
    )
  }),
  data.frame(
    label = "B_DORSALIS",
    group = "outgroup",
    stringsAsFactors = FALSE
  )
) %>%
  filter(label %in% host_tree_pruned$tip.label) %>%
  mutate(
    focal = label %in% c(
      focal_species,
      "B_DORSALIS"
    ),
    species_label = case_when(
      label == "B_DORSALIS" ~ "B. dorsalis",
      TRUE ~ paste0(
        "D. ",
        tolower(sub("^D_", "", label))
      )
    ),
    group = factor(
      group,
      levels = c(
        "melanogaster",
        "suzukii",
        "montium",
        "ananassae",
        "outgroup"
      )
    )
  )


group_cols <- c(
  melanogaster = "#4C78A8",
  suzukii = "#F58518",
  montium = "#54A24B",
  ananassae = "#B279A2",
  outgroup = "grey30"
)

p0 <- ggtree(
  host_tree_pruned,
  layout = "rectangular",
  branch.length = "none"
)

p0$data <- p0$data %>%
  left_join(host_meta, by = "label")

tip_nonfocal <- p0$data %>%
  filter(isTip, !focal)

tip_focal <- p0$data %>%
  filter(isTip, focal)

p_host <- p0 +
  geom_tippoint(
    data = tip_nonfocal,
    aes(fill = group),
    shape = 21,
    size = 2.2,
    color = "black",
    stroke = 0.25
  ) +
  geom_tippoint(
    data = tip_focal,
    aes(fill = group),
    shape = 21,
    size = 4.2,
    color = "black",
    stroke = 0.8
  ) +
  geom_tiplab(
    data = tip_nonfocal,
    aes(
      label = species_label,
      color = group
    ),
    fontface = "italic",
    size = 3,
    offset = 0.15
  ) +
  geom_tiplab(
    data = tip_focal,
    aes(
      label = species_label,
      color = group
    ),
    fontface = "bold.italic",
    size = 3.4,
    offset = 0.25
  ) +
  scale_fill_manual(
    values = group_cols,
    name = "Drosophila\nsubgroup"
  ) +
  scale_color_manual(
    values = group_cols,
    name = "Drosophila\nsubgroup"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.02, 0.3))
  ) +
  guides(
    color = "none",
    fill = guide_legend(
      override.aes = list(
        shape = 21,
        size = 4,
        color = "black"
      )
    )
  ) +
  theme_tree() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

pdf("outputs/plots/droso_spec_subgroups_bd.pdf", width = 6.5, height = 5)
print(p_host)
dev.off()
