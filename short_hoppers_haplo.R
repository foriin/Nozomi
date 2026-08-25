library(Biostrings)
library(DECIPHER)
library(dplyr)
library(tidyr)
library(ggplot2)

#### Inputs ####

load("outputs/RData/short_hoppers_seq_dm6loc_group.RData", verbose = TRUE)

## shops.fa = DNAStringSet with short Hopper copies
## mcols(shops.fa)$group = genomic coordinate group

longhop <- readDNAStringSet("~/Work/projects/DrosoTE/popTE/earlgrey/a6_full_hopper_putat.fa")

stopifnot("group" %in% colnames(mcols(shops.fa)))
stopifnot(length(longhop) >= 1)

## Use one long Hopper as reference
longhop <- longhop[1]

if (is.null(names(shops.fa)) || any(names(shops.fa) == "")) {
  names(shops.fa) <- paste0("short_", seq_along(shops.fa))
}

if (is.null(names(longhop)) || names(longhop)[1] == "") {
  names(longhop) <- "long_hopper_ref"
}

short_names <- names(shops.fa)
long_name <- names(longhop)[1]

stopifnot(!anyDuplicated(c(short_names, long_name)))


#### Align short Hoppers + one long reference ####


all.seq <- c(shops.fa, longhop)

all.aln <- AlignSeqs(all.seq)

aln_mat <- as.matrix(all.aln)

short_mat <- aln_mat[short_names, , drop = FALSE]
long_ref <- aln_mat[long_name, ]

stopifnot(nrow(short_mat) == length(shops.fa))
stopifnot(length(long_ref) == ncol(short_mat))

bases <- c("A", "C", "G", "T", "a", "c", "g", "t")


#### Detect universal large central deletion ####


## Long reference has base, short copies mostly have gaps
short_gap_frac <- colMeans(short_mat == "-")
long_has_base <- long_ref %in% bases

universal_del_cols <- long_has_base & short_gap_frac >= 0.90

## Optional: inspect the universal deletion in alignment coordinates
universal_del_runs <- rle(universal_del_cols)

universal_del_df <- data.frame(
  run_id = seq_along(universal_del_runs$lengths),
  value = universal_del_runs$values,
  aln_start = cumsum(c(1, head(universal_del_runs$lengths, -1))),
  aln_end = cumsum(universal_del_runs$lengths),
  length = universal_del_runs$lengths
) %>%
  filter(value) %>%
  arrange(desc(length))

universal_del_df


#### Build haplotypes ####


## STRICT haplotype:
## SNPs + small deletions + small insertions, but excluding universal central deletion
strict_keep_cols <- !universal_del_cols & colSums(short_mat != "-") > 0

short_strict_mat <- short_mat[, strict_keep_cols, drop = FALSE]

if (ncol(short_strict_mat) > 0) {
  strict_informative_cols <- apply(short_strict_mat, 2, function(x) {
    length(unique(x)) > 1
  })
} else {
  strict_informative_cols <- logical(0)
}

short_strict_var_mat <- short_strict_mat[, strict_informative_cols, drop = FALSE]

if (ncol(short_strict_var_mat) == 0) {
  strict_hap_string <- rep("NO_VARIATION", nrow(short_strict_var_mat))
  names(strict_hap_string) <- rownames(short_strict_var_mat)
} else {
  strict_hap_string <- apply(short_strict_var_mat, 1, paste0, collapse = "")
}

strict_hap_id <- paste0(
  "HAP",
  sprintf("%04d", as.integer(factor(strict_hap_string)))
)

names(strict_hap_id) <- names(strict_hap_string)

## SNP-only haplotype:
## use only long-reference columns outside universal deletion
## gaps are converted to N, so small indels do not define haplotypes
snp_keep_cols <- long_has_base & !universal_del_cols & colSums(short_mat != "-") > 0

short_snp_mat <- short_mat[, snp_keep_cols, drop = FALSE]
short_snp_mat[short_snp_mat == "-"] <- "N"

if (ncol(short_snp_mat) > 0) {
  snp_informative_cols <- apply(short_snp_mat, 2, function(x) {
    length(unique(x)) > 1
  })
} else {
  snp_informative_cols <- logical(0)
}

short_snp_var_mat <- short_snp_mat[, snp_informative_cols, drop = FALSE]

if (ncol(short_snp_var_mat) == 0) {
  snp_hap_string <- rep("NO_VARIATION", nrow(short_snp_var_mat))
  names(snp_hap_string) <- rownames(short_snp_var_mat)
} else {
  snp_hap_string <- apply(short_snp_var_mat, 1, paste0, collapse = "")
}

snp_hap_id <- paste0(
  "SNPHAP",
  sprintf("%04d", as.integer(factor(snp_hap_string)))
)

names(snp_hap_id) <- names(snp_hap_string)

## Add haplotypes back to DNAStringSet metadata
mcols(shops.fa)$strict_haplotype <- strict_hap_id[names(shops.fa)]
mcols(shops.fa)$snp_haplotype <- snp_hap_id[names(shops.fa)]


#### Metadata table ####


haplo_df <- data.frame(
  copy_id = names(shops.fa),
  group = as.character(mcols(shops.fa)$group),
  strict_haplotype = as.character(mcols(shops.fa)$strict_haplotype),
  snp_haplotype = as.character(mcols(shops.fa)$snp_haplotype),
  stringsAsFactors = FALSE
)

genome_col <- intersect(
  c("genome", "strain", "sample", "assembly"),
  colnames(mcols(shops.fa))
)[1]

if (!is.na(genome_col)) {
  haplo_df$genome <- as.character(mcols(shops.fa)[[genome_col]])
} else {
  haplo_df$genome <- NA_character_
}


#### Haplotype-level stats ####

strict_haplo_stats <- haplo_df %>%
  group_by(strict_haplotype) %>%
  summarise(
    n_copies = n(),
    n_groups = n_distinct(group),
    n_genomes = n_distinct(genome),
    groups = paste(sort(unique(group)), collapse = ";"),
    genomes = paste(sort(unique(genome)), collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(desc(n_groups), desc(n_copies))

snp_haplo_stats <- haplo_df %>%
  group_by(snp_haplotype) %>%
  summarise(
    n_copies = n(),
    n_groups = n_distinct(group),
    n_genomes = n_distinct(genome),
    groups = paste(sort(unique(group)), collapse = ";"),
    genomes = paste(sort(unique(genome)), collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(desc(n_groups), desc(n_copies))

strict_haplo_stats
snp_haplo_stats


#### Genomic group-level stats ####


group_haplo_stats <- haplo_df %>%
  group_by(group) %>%
  summarise(
    n_copies = n(),
    n_genomes = n_distinct(genome),
    n_strict_haplotypes = n_distinct(strict_haplotype),
    n_snp_haplotypes = n_distinct(snp_haplotype),
    strict_haplotypes = paste(sort(unique(strict_haplotype)), collapse = ";"),
    snp_haplotypes = paste(sort(unique(snp_haplotype)), collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(desc(n_genomes), desc(n_strict_haplotypes), desc(n_snp_haplotypes))

group_haplo_stats


#### Haplotype x genomic group matrix ####


haplo_group_matrix <- haplo_df %>%
  count(strict_haplotype, group, name = "n") %>%
  pivot_wider(
    names_from = group,
    values_from = n,
    values_fill = 0
  )

haplo_group_matrix

#### Useful subsets ####

###### Haplotypes found at multiple independent genomic positions ######
multi_group_haplotypes <- strict_haplo_stats %>%
  filter(n_groups > 1) %>%
  arrange(desc(n_groups), desc(n_copies))

multi_group_haplotypes

###### Genomic groups containing several strict haplotypes ######
multi_haplo_groups <- group_haplo_stats %>%
  filter(n_strict_haplotypes > 1) %>%
  arrange(desc(n_genomes), desc(n_strict_haplotypes))

multi_haplo_groups


#### Plots ####


ggplot(strict_haplo_stats, aes(x = n_groups, y = n_copies)) +
  geom_point(alpha = 0.8) +
  theme_classic() +
  labs(
    x = "Number of genomic groups carrying strict haplotype",
    y = "Number of copies with strict haplotype"
  )

ggplot(group_haplo_stats, aes(x = n_genomes, y = n_strict_haplotypes)) +
  geom_point(alpha = 0.8) +
  theme_classic() +
  labs(
    x = "Number of genomes carrying genomic group",
    y = "Number of strict haplotypes in group"
  )

ggplot(group_haplo_stats, aes(x = n_genomes, y = n_snp_haplotypes)) +
  geom_point(alpha = 0.8) +
  theme_classic() +
  labs(
    x = "Number of genomes carrying genomic group",
    y = "Number of SNP haplotypes in group"
  )

#### Look at different haplotypes ####

snphap0378 <- AlignSeqs(shops.fa[shops.fa@elementMetadata$snp_haplotype == 'SNPHAP0378'])
snphap0378@elementMetadata <- shops.fa[names(snphap0378)]@elementMetadata
snphap0378@elementMetadata

snphap0421 <- AlignSeqs(shops.fa[shops.fa@elementMetadata$snp_haplotype == 'SNPHAP0421'])
snphap0421@elementMetadata <- shops.fa[names(snphap0421)]@elementMetadata
snphap0421@elementMetadata

snphap0373 <- AlignSeqs(shops.fa[shops.fa@elementMetadata$snp_haplotype == 'SNPHAP0373'])
snphap0373@elementMetadata <- shops.fa[names(snphap0373)]@elementMetadata
snphap0373@elementMetadata

###### Add TSD to short hoppers seq md ######
load("outputs/RData/short_hoppers_origloc_seq.RData", verb = T)
rownames(hopsh.tsd) <- hopsh.tsd$name
mcols(shops.fa) <- cbind(mcols(shops.fa), hopsh.tsd[names(shops.fa),2:3])

shops.md <- mcols(shops.fa) %>% as.data.frame()

#### Save objects ####
save(universal_del_df,
     strict_haplo_stats,
     snp_haplo_stats,
     file = "outputs/RData/short_hop_haplo_stats.RData")
save(shops.fa,
     file = "outputs/RData/short_hoppers_seq_dm6loc_group_haplo_tsd.RData")

#### Optional: save tables ####

# write.csv(haplo_df, "hopper_copy_haplotypes.csv", row.names = FALSE)
# write.csv(strict_haplo_stats, "hopper_strict_haplotype_stats.csv", row.names = FALSE)
# write.csv(snp_haplo_stats, "hopper_snp_haplotype_stats.csv", row.names = FALSE)
# write.csv(group_haplo_stats, "hopper_group_haplotype_stats.csv", row.names = FALSE)
# write.csv(haplo_group_matrix, "hopper_haplotype_by_group_matrix.csv", row.names = FALSE)