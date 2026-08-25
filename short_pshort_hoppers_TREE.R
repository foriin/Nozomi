library(ape)
library(phangorn)
library(dplyr)
library(tibble)
library(ggtree)
library(ggplot2)

a6lh.uq <- unique(a6lh.seq)
a6lh.l <- subseq(a6lh.uq, 1, c(1034, 1016, 1016, 1011, 1034, 1034))
a6.lh.r <- subseq(a6lh.uq, c(2415, 2373, 2374, 2369, 2415, 2415), width(a6lh.uq))
a6.pshop <- DNAStringSet(paste0(a6lh.l, a6.lh.r))[1:4]
names(a6.pshop) <- c(
  "A6_longhop",
  "chr3Lfix_1",
  "chr3Lfix_2",
  "chr3Lfix_3"
)
a6sl.aln <- AlignSeqs(c(a6lh.uq, reverseComplement(unique(shops.a6))))
BrowseSeqs(a6sl.aln)

writeXStringSet(c(a6.pshop, reverseComplement(unique(shops.a6))), "outputs/fasta/a6_short_hoppers_unique_pshort_added.fa")
# dir.create("outputs/msa")

# mafft --auto --thread 8 outputs/fasta/a6_short_hoppers_unique_pshort_added.fa > outputs/msa/a6_short_hoppers_pshort_add.mafft.fa

# iqtree3 \
# -s outputs/msa/a6_short_hoppers_pshort_add.mafft.fa \
# -m GTR+G \
# -B 1000 \
# -alrt 1000 \
# -nt AUTO \
# --prefix outputs/msa/a6_short_pshort_hoppers_iqtree


tr <- read.tree("outputs/msa/a6_short_pshort_hoppers_iqtree.contree")

p_tree <- ggtree(tr, layout = "unrooted") +
  geom_tiplab(size = 2, align = FALSE) +
  theme_tree2() +
  ggtitle("Phylogeny of A6 short Hopper insertions")

print(p_tree)

p_tree_circ <- ggtree(tr, layout = "circular") +
  geom_tiplab(size = 1.8) +
  ggtitle("Phylogeny of A6 short Hopper insertions")

print(p_tree_circ)


#### Color by K80 ####
aln <- read.dna("outputs/msa/a6_short_hoppers_pshort_add.mafft.fa", format = "fasta")

k80_mat <- dist.dna(
  aln,
  model = "K80",
  pairwise.deletion = FALSE,
  as.matrix = TRUE
)

hc <- hclust(as.dist(k80_mat), method = "average")

cluster_df <- tibble(
  label = names(cutree(hc, k = 3)),
  hopper_cluster = factor(cutree(hc, k = 3))
)

p_tree_cluster <- ggtree(tr, layout = "circular") %<+% cluster_df +
  geom_tippoint(aes(color = hopper_cluster), size = 2) +
  geom_tiplab(size = 1.8) +
  ggtitle("A6 short Hopper tree colored by sequence cluster") +
  theme(legend.position = "right")

print(p_tree_cluster)

#### Add bootstrap, color by group size (i.e. in how many genomes this insertion was found) ####

setdiff(tr$tip.label, names(shops.a6))

a6mcols <- mcols(shops.a6[tr$tip.label]) %>%
  as.data.frame() %>%
  mutate(group = as.character(group)) %>% 
  rownames_to_column("label") %>%
  left_join(group_haplo_stats, by = "group")

tip_meta <- a6mcols %>%
  mutate(
    group_size = n_genomes,
    tip_label = paste0("group: ", group, " group size: ", group_size)
  ) %>%
  select(label, group, group_size, tip_label) %>%
  distinct(label, .keep_all = TRUE)

stopifnot(!anyDuplicated(tip_meta$label))
stopifnot(all(tr$tip.label %in% tip_meta$label))

p0 <- ggtree(tr, layout = "unrooted")

p0$data <- p0$data %>%
  left_join(tip_meta, by = "label") %>%
  mutate(node_support = suppressWarnings(as.numeric(label)))

p_tree_boot <- p0 +
  geom_tippoint(
    aes(fill = group_size),
    shape = 21,
    color = "black",
    stroke = 0.25,
    size = 2.6,
    na.rm = TRUE
  ) +
  geom_tiplab(
    aes(label = tip_label),
    size = 1.8,
    na.rm = TRUE
  ) +
  geom_text2(
    aes(
      subset = !isTip & !is.na(node_support) & node_support >= 70,
      label = label
    ),
    size = 2,
    hjust = -0.2
  ) +
  scale_fill_gradient(
    low = "white",
    high = "darkred",
    limits = c(0, 44),
    name = "Group size",
    na.value = "grey80"
  ) +
  ggtitle("A6 short Hopper phylogeny colored by haplotype group size") +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

print(p_tree_boot)


#### Compare clusters ####
cl <- cutree(hc, k = 3)

pair_k80_df <- as.data.frame(as.table(k80_mat)) %>%
  dplyr::rename(ins1 = Var1, ins2 = Var2, k80 = Freq) %>%
  filter(as.character(ins1) < as.character(ins2)) %>%
  mutate(
    cl1 = cl[as.character(ins1)],
    cl2 = cl[as.character(ins2)],
    comparison = ifelse(cl1 == cl2, "within_cluster", "between_cluster")
  )

pair_k80_df %>%
  group_by(comparison) %>%
  summarise(
    n = n(),
    median_k80 = median(k80, na.rm = TRUE),
    mean_k80 = mean(k80, na.rm = TRUE),
    min_k80 = min(k80, na.rm = TRUE),
    max_k80 = max(k80, na.rm = TRUE),
    .groups = "drop"
  )


ggplot(pair_k80_df, aes(x = comparison, y = k80)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 1) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Pairwise K80 divergence",
    title = "Within- vs between-cluster divergence of short Hopper insertions"
  )

#### Check groups among 2 clusters ####
tree_dist <- cophenetic.phylo(tr)

hc2 <- hclust(as.dist(tree_dist), method = "average")

tip_clusters <- data.frame(
  label = names(cutree(hc2, k = 2)),
  tree_cluster = factor(cutree(hc2, k = 2))
)

tip_meta2 <- tip_meta %>%
  left_join(tip_clusters, by = "label")

####### plot tree itself ######

p0 <- ggtree(tr, layout = "unrooted")

p0$data <- p0$data %>%
  left_join(tip_meta2, by = "label") %>%
  mutate(node_support = suppressWarnings(as.numeric(label)))

p_tree_cluster <- p0 +
  geom_tippoint(
    aes(fill = tree_cluster),
    shape = 21,
    color = "black",
    stroke = 0.25,
    size = 2.8,
    na.rm = TRUE
  ) +
  geom_text2(
    aes(
      subset = !isTip & !is.na(node_support) & node_support >= 70,
      label = label
    ),
    size = 2,
    hjust = -0.2
  ) +
  ggtitle("A6 short Hopper phylogeny split by cutree(k = 2)") +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

print(p_tree_cluster)

#### Box plot of group sizes per cluster ####
ggplot(tip_meta2, aes(x = tree_cluster, y = group_size, fill = tree_cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, alpha = 0.6, size = 2) +
  theme_bw() +
  labs(
    x = "Tree cluster",
    y = "Group size",
    title = "Haplotype group size by Hopper phylogenetic cluster"
  ) +
  theme(legend.position = "none")


#### Build unrooted tree, add gen loc groups, haplo size ####
seq_df <- data.frame(
  label = names(shops.a6),
  seq = as.character(shops.a6),
  stringsAsFactors = FALSE
)

collapse_df <- seq_df %>%
  dplyr::group_by(seq) %>%
  dplyr::summarise(
    representative = dplyr::first(label),
    n_identical = n(),
    members = paste(label, collapse = ";"),
    .groups = "drop"
  )

unique_shops <- shops.a6[match(collapse_df$representative, names(shops.a6))]
names(unique_shops) <- collapse_df$representative

mcols(unique_shops)$n_identical <- collapse_df$n_identical
mcols(unique_shops)$members <- collapse_df$members

# writeXStringSet(unique_shops, "outputs/fasta/a6_short_hoppers_unique.fa")

tip_counts <- data.frame(
  label = collapse_df$representative,
  n_identical = collapse_df$n_identical,
  members = collapse_df$members,
  stringsAsFactors = FALSE
)

stopifnot(all(tr$tip.label %in% tip_counts$label))

tip_meta <- tip_meta %>%
  left_join(tip_counts, by = "label")

p0 <- ggtree(tr, layout = "unrooted")

p0$data <- p0$data %>%
  left_join(tip_meta, by = "label") %>%
  mutate(node_support = suppressWarnings(as.numeric(label)))

p_tree_counts <- p0 +
  geom_tippoint(
    aes(size = n_identical, fill = group_size),
    shape = 21,
    color = "black",
    stroke = 0.25,
    alpha = 0.9,
    na.rm = TRUE
  ) +
  geom_text2(
    aes(
      subset = !isTip & !is.na(node_support) & node_support >= 70,
      label = label
    ),
    size = 2,
    hjust = -0.2
  ) +
  scale_fill_gradient(
    low = "white",
    high = "darkred",
    limits = c(0, 44),
    name = "Group size",
    na.value = "grey80"
  ) +
  scale_size_continuous(
    name = "Identical copies",
    range = c(1.5, 7)
  ) +
  ggtitle("A6 short Hopper phylogeny: collapsed identical sequences") +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

print(p_tree_counts)


#### Recalc group size ####
seq_meta <- cbind(
  data.frame(
    label = names(shops.a6),
    seq = as.character(shops.a6),
    stringsAsFactors = FALSE
  ),
  as.data.frame(mcols(shops.a6)) %>% mutate(group = as.character(group))
) %>%
  left_join(group_haplo_stats, by = "group") %>%
  mutate(group_size_original = n_genomes)

collapsed_meta <- seq_meta %>%
  group_by(seq) %>%
  summarise(
    label = dplyr::first(label),
    n_identical = n(),
    group_size = max(group_size_original, na.rm = TRUE),
    groups = paste(unique(group), collapse = ";"),
    members = paste(label, collapse = ";"),
    .groups = "drop"
  )

stopifnot(all(tr$tip.label %in% collapsed_meta$label))

tip_meta_collapsed <- collapsed_meta %>%
  filter(label %in% tr$tip.label) %>%
  mutate(
    tip_label = paste0(
      "max group size: ", group_size,
      " identical: ", n_identical
    )
  )

p0 <- ggtree(tr, layout = "unrooted")

# one of 'rectangular', 'dendrogram', 'slanted', 'ellipse', 'roundrect', 'fan', 'circular', 'inward_circular', 'radial', 'equal_angle', 'daylight' or 'ape'
p0 <- ggtree(tr, layout = "unrooted")

p0$data <- p0$data %>%
  left_join(tip_meta_collapsed, by = "label") %>%
  mutate(node_support = suppressWarnings(as.numeric(label)))

p_tree_collapsed <- p0 +
  geom_tippoint(
    aes(fill = group_size, size = n_identical),
    shape = 21,
    color = "black",
    stroke = 0.25,
    alpha = 0.9,
    na.rm = TRUE
  ) +
  geom_text2(
    aes(
      subset = !isTip & !is.na(node_support) & node_support >= 70,
      label = label
    ),
    size = 2,
    hjust = -0.2
  ) +
  scale_fill_gradient(
    low = "white",
    high = "darkred",
    limits = c(0, 44),
    name = "Max group size",
    na.value = "grey80"
  ) +
  scale_size_continuous(
    name = "Identical copies",
    range = c(1.5, 7)
  ) +
  ggtitle("A6 short Hopper phylogeny: collapsed identical sequences") +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

print(p_tree_collapsed)

pdf("outputs/plots/A6_short_hoppers_good_tree_locgroup_size_id_copies.pdf",
    width = 9, height = 8)
print(p_tree_collapsed)
dev.off()

collapsed_meta %>%
  arrange(desc(group_size), desc(n_identical)) %>%
  select(label, group_size, n_identical, groups, members) %>%
  head(30)


