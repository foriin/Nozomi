library(Biostrings)
library(DECIPHER)
library(ape)
library(dplyr)
library(ggplot2)

## shops.fa = DNAStringSet
## group is stored as mcols(shops.fa)$group

load("outputs/RData/short_hoppers_seq_dm6loc_group.RData", verbose = TRUE)

longhop <- readDNAStringSet("~/Work/projects/DrosoTE/popTE/earlgrey/a6_full_hopper_putat.fa")

stopifnot("group" %in% colnames(mcols(shops.fa)))
stopifnot(!is.null(names(shops.fa)))
stopifnot(!anyDuplicated(names(shops.fa)))

## make sure long Hopper names cannot collide with short Hopper names
if (is.null(names(longhop)) || any(names(longhop) == "")) {
  names(longhop) <- paste0("longhop_", seq_along(longhop))
}

names(longhop) <- paste0("LONGHOP__", make.unique(names(longhop)))

stopifnot(!anyDuplicated(c(names(shops.fa), names(longhop))))

short_names <- names(shops.fa)
long_names <- names(longhop)

## combine short + long and align once
all.seq <- c(shops.fa, longhop)

all.aln <- AlignSeqs(all.seq)

## split aligned object back into short and long
shops.aln <- all.aln[short_names]
longhop.aln <- all.aln[long_names]

stopifnot(identical(names(shops.aln), short_names))
stopifnot(identical(names(longhop.aln), long_names))

## metadata matched to aligned short sequences
group_vec <- mcols(shops.fa)$group
names(group_vec) <- names(shops.fa)
group_vec <- group_vec[names(shops.aln)]

stopifnot(identical(names(group_vec), names(shops.aln)))

## optional genome/sample column detection
genome_col <- intersect(
  c("genome", "strain", "sample", "assembly"),
  colnames(mcols(shops.fa))
)[1]

if (length(genome_col) == 0 || is.na(genome_col)) {
  genome_vec <- rep(NA_character_, length(shops.fa))
  names(genome_vec) <- names(shops.fa)
} else {
  genome_vec <- as.character(mcols(shops.fa)[[genome_col]])
  names(genome_vec) <- names(shops.fa)
}

genome_vec <- genome_vec[names(shops.aln)]

## convert aligned sequences to DNAbin
## pairwise.deletion = TRUE means gap-only central deletion columns are ignored pairwise
all.bin <- as.DNAbin(as.matrix(all.aln))

k80.mat <- dist.dna(
  all.bin,
  model = "K80",
  pairwise.deletion = TRUE,
  as.matrix = TRUE
)

## short-short and short-long K80 matrices
k80.short <- k80.mat[short_names, short_names, drop = FALSE]
k80.short.long <- k80.mat[short_names, long_names, drop = FALSE]

## helper for within-group pairwise K80
k80_within_stats <- function(ids) {
  if (length(ids) < 2) {
    return(data.frame(
      mean_k80 = NA_real_,
      median_k80 = NA_real_,
      max_k80 = NA_real_
    ))
  }
  
  d <- k80.short[ids, ids, drop = FALSE]
  vals <- d[upper.tri(d)]
  vals <- vals[is.finite(vals)]
  
  if (length(vals) == 0) {
    return(data.frame(
      mean_k80 = NA_real_,
      median_k80 = NA_real_,
      max_k80 = NA_real_
    ))
  }
  
  data.frame(
    mean_k80 = mean(vals),
    median_k80 = median(vals),
    max_k80 = max(vals)
  )
}

## helper for distance from group members to long Hopper copies
k80_longhop_stats <- function(ids) {
  d <- k80.short.long[ids, , drop = FALSE]
  vals <- as.numeric(d)
  vals <- vals[is.finite(vals)]
  
  if (length(vals) == 0) {
    return(data.frame(
      mean_k80_to_longhop = NA_real_,
      median_k80_to_longhop = NA_real_,
      min_k80_to_longhop = NA_real_
    ))
  }
  
  data.frame(
    mean_k80_to_longhop = mean(vals),
    median_k80_to_longhop = median(vals),
    min_k80_to_longhop = min(vals)
  )
}

## nearest long Hopper per short copy
nearest_longhop_per_copy <- apply(k80.short.long, 1, function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_character_)
  names(which.min(x))
})

nearest_longhop_dist_per_copy <- apply(k80.short.long, 1, function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  min(x)
})

## split short sequences by genomic group
ids_by_group <- split(names(shops.aln), group_vec)

## all groups, including singletons
k80_by_group <- lapply(names(ids_by_group), function(g) {
  ids <- ids_by_group[[g]]
  
  nearest_longhop_group <- nearest_longhop_per_copy[ids]
  nearest_longhop_group <- nearest_longhop_group[!is.na(nearest_longhop_group)]
  
  nearest_longhop_majority <- if (length(nearest_longhop_group) == 0) {
    NA_character_
  } else {
    names(sort(table(nearest_longhop_group), decreasing = TRUE))[1]
  }
  
  cbind(
    data.frame(
      group = g,
      n_copies = length(ids),
      n_genomes = if (all(is.na(genome_vec[ids]))) NA_integer_ else length(unique(genome_vec[ids])),
      copies = paste(ids, collapse = ";"),
      nearest_longhop_majority = nearest_longhop_majority,
      stringsAsFactors = FALSE
    ),
    k80_within_stats(ids),
    k80_longhop_stats(ids)
  )
}) |> bind_rows()

k80_by_group <- k80_by_group |>
  mutate(
    freq = n_genomes / 46,
    freq_bin = case_when(
      is.na(n_genomes) ~ NA_character_,
      n_genomes == 1 ~ "1 singleton",
      n_genomes <= 4 ~ "2-4 rare",
      n_genomes <= 19 ~ "5-19 intermediate",
      n_genomes <= 39 ~ "20-39 high",
      n_genomes >= 40 ~ "40-46 near-fixed"
    ),
    freq_bin = factor(
      freq_bin,
      levels = c(
        "1 singleton",
        "2-4 rare",
        "5-19 intermediate",
        "20-39 high",
        "40-46 near-fixed"
      )
    )
  ) |>
  arrange(desc(n_copies), desc(mean_k80))

k80_by_group <- merge(ahred.df, k80_by_group, by = "group")
k80_by_group <- k80_by_group[, -c(7,8,12)]

shops.groups.w.stats.gr <- makeGRangesFromDataFrame(k80_by_group, keep.extra.columns = T)


#### Plots ####

ggplot(k80_by_group, aes(x = n_genomes, y = median_k80)) +
  geom_point(aes(color = freq_bin, size = n_copies), alpha = 0.8) +
  theme_classic() +
  labs(
    x = "Number of genomes carrying insertion group",
    y = "Median within-group K80",
    color = "Frequency class",
    size = "Copies"
  )

ggplot(k80_by_group, aes(x = n_genomes, y = mean_k80_to_longhop)) +
  geom_point(aes(color = freq_bin, size = n_copies), alpha = 0.8) +
  theme_classic() +
  labs(
    x = "Number of genomes carrying insertion group",
    y = "Mean K80 distance to active long Hopper",
    color = "Frequency class",
    size = "Copies"
  )

ggplot(k80_by_group %>% filter(n_genomes > 1),
       aes(x = mean_k80_to_longhop, y = median_k80)) +
  geom_point(aes(color = freq_bin), alpha = 0.8) +
  theme_classic() +
  labs(
    x = "Mean K80 distance to active long Hopper",
    y = "Median within-group K80",
    color = "Frequency class",
    size = "Copies"
  )

ggplot(
  k80_by_group %>%
    filter(
      n_genomes > 1,
      !is.na(mean_k80_to_longhop),
      !is.na(median_k80),
      !is.na(freq_bin)
    ) %>%
    droplevels(),
  aes(x = mean_k80_to_longhop, y = median_k80)
) +
  geom_point(aes(color = freq_bin), alpha = 0.8) +
  geom_smooth(
    aes(color = freq_bin, group = freq_bin),
    method = "lm",
    se = FALSE,
    linewidth = 0.8
  ) +
  theme_classic() +
  labs(
    x = "Mean K80 distance to active long Hopper",
    y = "Median within-group K80",
    color = "Frequency class"
  )


#### Just group stats ####
shops.md <- shops.fa@elementMetadata %>% as_tibble()
shops.group.summary <- shops.md %>% group_by(group) %>% 
  summarize(strains = paste(strain, collapse = "; "),
            loc = first(loc), n_genomes = length(unique(strain)))

#### Save group gr with the stats ####
save(shops.groups.w.stats.gr, k80_by_group, file = "outputs/RData/short_hop_groups_GRanges.RData")


