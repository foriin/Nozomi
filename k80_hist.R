library(Biostrings)
library(DECIPHER)
library(ape)

get_k80_to_ref <- function(x, ref, set_name = NA_character_,
                           include_ref = F) {
  
  ## x = unaligned DNAStringSet
  
  if (is.character(ref)) {
    ref_name <- ref
  } else if (is.numeric(ref)) {
    ref_name <- names(x)[ref]
  } else {
    stop("ref must be sequence name or index")
  }
  
  aln <- AlignSeqs(x)
  
  dna <- ape::as.DNAbin(as.matrix(aln))
  
  k80 <- as.matrix(ape::dist.dna(
    dna,
    model = "K80",
    pairwise.deletion = TRUE,
    as.matrix = TRUE
  ))
  
  nm <- names(aln)
  ref_i <- match(ref_name, nm)
  
  if (is.na(ref_i)) {
    stop("Reference sequence not found after alignment: ", ref_name)
  }
  
  n_insertions <- suppressWarnings(
    as.integer(sub(".*_n_([0-9]+)$", "\\1", nm))
  )
  n_insertions[is.na(n_insertions)] <- 1
  
  df <- data.frame(
    set = set_name,
    sequence = nm,
    ref_sequence = ref_name,
    width = width(x)[match(nm, names(x))],
    n = n_insertions,
    k80 = k80[ref_i, ],
    stringsAsFactors = FALSE
  )
  if (!include_ref){
    df <- df %>% filter(sequence != ref_name)
  }
  
  df
  
}



#### Prepare datasets ####


a6.k80 <- get_k80_to_ref(
  c(reverseComplement(a6hopper.fa),
    longhop),
  ref = "A6_hopper_long",
  set_name = "Dmel_Hopper (A6)"
)

names(dsuz_ht_fa) <- paste("Dsuz_ilop", 1:length(dsuz_ht_fa), sep = "_")


dsuz_interloper.k80 <- get_k80_to_ref(
  dsuz_ht_fa,
  ref = "Dsuz_ilop_1",
  set_name = "Dsuz_Interloper",
  include_ref = T
)

dsuz_hops_aln <- AlignSeqs(dsuz_hopst_fa)
BrowseSeqs(dsuz_hops_aln)

dsuz_hopst_fa2 <- c(
  dsuz_hopst_fa[-69],
  reverseComplement(dsuz_hopst_fa[69])
)

###### Chopper in Dsuz ref genome ######

dsuz_chop.bed <- import.bed("/Users/artemilin/Work/projects/DrosoTE/Hopper_A6/teaid/Dsuz_ref_hopster/dsuz_hopster_ref.bed")
dsuz_chop.bed <- GenomicRanges::reduce(dsuz_chop.bed)
dsuz_ref.fa <- "/Users/artemilin/Work/projects/DrosoTE/Hopper_A6/genomes/Dsuz_ref_IsoJPT.fa"
dsuz_ref.genome <- FaFile(dsuz_ref.fa)
open(dsuz_ref.genome)

dsuz_chopper.ref <- getSeq(dsuz_ref.genome, dsuz_chop.bed)
dsuz_chopper.fa <- readDNAStringSet("/Users/artemilin/Work/projects/DrosoTE/Hopper_A6/DNA/dsuz_hopster.fa")
dsuz_chopper.ref <- c(reverseComplement(dsuz_chopper.ref), dsuz_chopper.fa)


dsuz_chopper.k80 <- get_k80_to_ref(
  dsuz_chopper.ref,
  ref = "Dsuz_hopster",
  set_name = "Dsuz_Chopper"
)

###### Chopper in romanian dsuz ######
coga.fa <- "/Users/artemilin/Work/projects/DrosoTE/Hopper_A6/genomes/Dsuz_GB_ls_coga.fa"
coganome <- FaFile(coga.fa)
open(coganome)

coga_chop.bed <- import.bed("/Users/artemilin/Work/projects/DrosoTE/Hopper_A6/teaid/Dsuz_coga_hopster/dsuz_hopster.bed")
coga_chop.bed <- GenomicRanges::reduce(coga_chop.bed)
coga_chop.fa <- getSeq(coganome, coga_chop.bed)
names(coga_chop.fa)[1] <- "Dsuz_chopper"

dsuz_coga_chopper.k80 <- get_k80_to_ref(
  coga_chop.fa,
  ref = "Dsuz_chopper",
  set_name = "Dsuz_Chopper (Bucharest)",
  include_ref = T
)


hops.dm6.fa2 <- c(reverseComplement(hops.dm6.fa), 
                  longhop)
dm6_hopper.k80 <- get_k80_to_ref(
  hops.dm6.fa2,
  ref = "A6_hopper_long",
  set_name = "Dmel_Hopper (dm6)"
)


ral426.k80 <- get_k80_to_ref(
  c(hopsh.all.fa[grepl("RAL-426", names(hopsh.all.fa))],
    longhop),
  ref = "A6_hopper_long",
  set_name = "RAL-426_Hopper"
)



a5.k80 <- get_k80_to_ref(
  c(shops.fa[mcols(shops.fa)$strain == "A5-VAG1"],
    longhop),
  ref = "A6_hopper_long",
  set_name = "A5-VAG1_Hopper"
)

bar.k80 <- get_k80_to_ref(
  c(reverseComplement(barhop),
    longhop),
  ref = "A6_hopper_long",
  set_name = "Dmel_Hopper (BL2969)"
)


###### Combine datasets ######


k80_df <- rbind(
  dm6_hopper.k80,
  a6.k80,
  bar.k80,
  ral426.k80,
  dsuz_interloper.k80,
  dsuz_chopper.k80,
  dsuz_coga_chopper.k80
)




k80_df <- k80_df[is.finite(k80_df$k80) & k80_df$k80 > 0 & k80_df$k80 < 0.5, ]

binwidth <- 0.005

#### Some old bullshit ####
k80_binned <- k80_df %>%
  filter(is.finite(k80)) %>%
  mutate(
    bin = floor(k80 / binwidth) * binwidth,
    bin_mid = bin + binwidth / 2
  ) %>%
  count(set, bin, bin_mid, name = "count") %>%
  group_by(set) %>%
  mutate(
    percent = 100 * count / sum(count)
  ) %>%
  ungroup()


ggplot(k80_binned, aes(x = bin_mid, y = percent, fill = set, color = set)) +
  geom_col(
    width = binwidth,
    position = "identity",
    alpha = 0.45,
    linewidth = 0.25
  ) +
  theme_classic() +
  labs(
    x = "K80 distance to reference sequence",
    y = "Sequences per bin (%)"
  )

ggplot(k80_df, aes(x = k80, fill = set, color = set)) +
  geom_density(alpha = 0.3, linewidth = 0.4) +
  theme_classic() +
  labs(
    x = "K80 distance to reference sequence",
    y = "Density"
  )

ggplot(k80_df, aes(x = k80, fill = set, color = set)) +
  geom_histogram(
    aes(y = after_stat(count / sum(count))),
    binwidth = 0.01,
    position = "identity",
    alpha = 0.35,
    linewidth = 0.3
  ) +
  theme_classic() +
  labs(
    x = "K80 distance to reference sequence",
    y = "Fraction of sequences"
  )


#### Get counts ####
binwidth <- 0.005

k80_counts <- k80_df %>%
  filter(is.finite(k80)) %>%
  mutate(
    bin = floor(k80 / binwidth) * binwidth,
    bin_mid = bin + binwidth / 2
  ) %>%
  group_by(set, bin_mid) %>%
  summarise(copy_number = sum(n), .groups = "drop")


#### Plot Histos ####

###### Hopper, chopper, interloper ####

pal <- c(
  "Dmel_Hopper (dm6)"    = "#2E7D32",  # same in both plots
  "Dsuz_Chopper"         = "#AB91E8",
  "Dsuz_Chopper (Bucharest)" = "#3B2168",
  "Dsuz_Interloper"      = "#E69F00",
  "Dmel_Hopper (A6)"     = "#4E79A7",
  "Dmel_Hopper (BL2969)" = "#D95F5F"
)

te_breaks <- c(
  "Dmel_Hopper (dm6)",
  "Dsuz_Chopper",
  "Dsuz_Interloper"
)

te_labels <- c(
  "Dmel_Hopper (dm6)" = "Hopper (dm6)",
  "Dsuz_Chopper" = "Chopper (Dsuz)",
  "Dsuz_Interloper" = "Interloper (Dsuz)"
)

p_te <- ggplot(
  k80_counts %>% filter(set %in% te_breaks),
  aes(x = bin_mid, y = copy_number, fill = set, color = set)
) +
  geom_col(
    width = binwidth,
    position = "identity",
    alpha = 0.35,
    linewidth = 0.25
  ) +
  scale_fill_manual(
    name = "TE",
    values = pal,
    breaks = te_breaks,
    labels = te_labels
  ) +
  scale_color_manual(
    name = "TE",
    values = pal,
    breaks = te_breaks,
    labels = te_labels
  ) +
  theme_classic() +
  labs(
    x = "K80 distance to reference sequence",
    y = "Copy number per bin"
  )

p_te


###### Hopper in dm6, A6 and BL2969 ######


genome_breaks <- c(
  "Dmel_Hopper (dm6)",
  "Dmel_Hopper (A6)",
  "Dmel_Hopper (BL2969)"
)

genome_labels <- c(
  "Dmel_Hopper (dm6)" = "dm6",
  "Dmel_Hopper (A6)" = "A6",
  "Dmel_Hopper (BL2969)" = "BL2969"
)

p_genome <- ggplot(
  k80_counts %>% filter(set %in% genome_breaks),
  aes(x = bin_mid, y = copy_number, fill = set, color = set)
) +
  geom_col(
    width = binwidth,
    position = "identity",
    alpha = 0.35,
    linewidth = 0.25
  ) +
  scale_fill_manual(
    name = "Genome",
    values = pal,
    breaks = genome_breaks,
    labels = genome_labels
  ) +
  scale_color_manual(
    name = "Genome",
    values = pal,
    breaks = genome_breaks,
    labels = genome_labels
  ) +
  theme_classic() +
  labs(
    x = "K80 distance to reference sequence",
    y = "Copy number per bin"
  )

p_genome

###### Chopper in ref and Bucharest genomes ######

chopper_breaks <- c(
  "Dsuz_Chopper",
  "Dsuz_Chopper (Bucharest)"
)

chopper_labels <- c(
  "Dsuz_Chopper" = "ref",
  "Dsuz_Chopper (Bucharest)" = "GB-ls-coga4 (Bucharest)"
)

p_chopper <- ggplot(
  k80_counts %>% filter(set %in% chopper_breaks),
  aes(x = bin_mid, y = copy_number, fill = set, color = set)
) +
  geom_col(
    width = binwidth,
    position = "identity",
    alpha = 0.55,
    linewidth = 0.25
  ) +
  scale_fill_manual(
    name = "Genome",
    values = pal,
    breaks = chopper_breaks,
    labels = chopper_labels
  ) +
  scale_color_manual(
    name = "Genome",
    values = pal,
    breaks = chopper_breaks,
    labels = chopper_labels
  ) +
  theme_classic() +
  labs(
    x = "K80 distance to reference sequence",
    y = "Copy number per bin"
  )

p_chopper



pdf("outputs/plots/te_k80_hist_te.pdf", width = 5, height = 3)
p_te
dev.off()

pdf("outputs/plots/te_k80_hist_genome.pdf", width = 5, height = 3)
p_genome
dev.off()
