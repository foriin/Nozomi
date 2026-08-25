library(Biostrings)
library(DECIPHER)
library(ape)
library(dplyr)

## shops.fa = DNAStringSet
## group is stored as mcols(shops.fa)$group

load('outputs/RData/short_hoppers_seq_dm6loc_group.RData', verb = T)
longhop <- readDNAStringSet("~/Work/projects/DrosoTE/popTE/earlgrey/a6_full_hopper_putat.fa")

stopifnot("group" %in% colnames(mcols(shops.fa)))

## align all sequences once
shops.aln <- AlignSeqs(shops.fa)

## keep group metadata matched by names
group_vec <- mcols(shops.fa)$group
names(group_vec) <- names(shops.fa)

group_vec <- group_vec[names(shops.aln)]

stopifnot(identical(names(group_vec), names(shops.aln)))

k80_stats_group <- function(x) {
  x_mat <- as.matrix(x)
  x_bin <- as.DNAbin(x_mat)
  
  d <- dist.dna(
    x_bin,
    model = "K80",
    pairwise.deletion = TRUE,
    as.matrix = TRUE
  )
  
  vals <- d[upper.tri(d)]
  
  data.frame(
    mean_k80 = mean(vals, na.rm = TRUE),
    median_k80 = median(vals, na.rm = TRUE),
    max_k80 = max(vals, na.rm = TRUE)
  )
}

## split sequence indices by group
idx_by_group <- split(seq_along(shops.aln), group_vec)

## keep only groups with >2 members
idx_by_group <- idx_by_group[lengths(idx_by_group) > 2]

k80_by_group <- lapply(names(idx_by_group), function(g) {
  idx <- idx_by_group[[g]]
  x <- shops.aln[idx]
  
  cbind(
    data.frame(
      group = g,
      n = length(x),
      stringsAsFactors = FALSE
    ),
    k80_stats_group(x)
  )
}) |> bind_rows()

k80_by_group <- k80_by_group |> arrange(desc(mean_k80))

k80_by_group

k80_by_group_loc <- merge(ahred.df, k80_by_group, by = "group")

#### Box Plots for distances between copies at the same loc and copies from diff locs ####

gr <- allshops.dm6
seqs <- shops.fa

common <- intersect(names(gr), names(seqs))

gr <- gr[common]
seqs <- seqs[common]
gr <- gr[names(seqs)]

stopifnot(identical(names(gr), names(seqs)))
stopifnot("group" %in% colnames(mcols(gr)))

seqs.aln <- AlignSeqs(seqs)

d <- as.matrix(stringDist(seqs.aln, method = "hamming"))
d <- d / width(seqs.aln)[1]

meta <- data.frame(
  name = names(seqs.aln),
  group = as.factor(mcols(gr)$group),
  stringsAsFactors = FALSE
)

pair.df <- as.data.frame(as.table(d)) |>
  rename(seq1 = Var1, seq2 = Var2, dist = Freq) |>
  filter(as.character(seq1) < as.character(seq2)) |>
  left_join(meta, by = c("seq1" = "name")) |>
  rename(group1 = group) |>
  left_join(meta, by = c("seq2" = "name")) |>
  rename(group2 = group) |>
  mutate(same_group = group1 == group2)

ggplot(pair.df, aes(same_group, dist)) +
  geom_boxplot(outlier.alpha = 0.15) +
  theme_classic() +
  labs(
    x = "Same genomic group",
    y = "Sequence distance"
  )

wilcox.test(dist ~ same_group, data = pair.df)

pair.df |>
  group_by(same_group) |>
  summarise(
    n = n(),
    median_dist = median(dist),
    mean_dist = mean(dist),
    .groups = "drop"
  )



library(phangorn)
library(ape)

haplo.fa <- unique(shops.fa)
names(haplo.fa) <- mcols(haplo.fa)$strict_haplotype
names(haplo.fa) <- make.unique(names(haplo.fa))
haplo.aln <- AlignSeqs(haplo.fa)

dna <- as.DNAbin(haplo.aln)

dm <- dist.dna(dna, model = "raw", pairwise.deletion = TRUE)
tree <- nj(dm)

plot(tree, cex = 0.4)
write.tree(tree, file = "outputs/phylo/hopper_haplo_tree.nwk")
writeXStringSet(haplo.aln, file = "outputs/phylo/hopper_haplo_msa.fa")
