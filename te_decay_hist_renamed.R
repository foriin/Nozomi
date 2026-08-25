library(Biostrings)
library(DECIPHER)
library(ape)

get_te_decay_to_ref <- function(
    x,
    ref,
    set_name = NA_character_,
    k80_cap = 0.5,
    k80_weight = 1,
    length_weight = 1,
    include_ref = FALSE
) {
  
  if (is.character(ref)) {
    ref_name <- ref
  } else if (is.numeric(ref)) {
    ref_name <- names(x)[ref]
  } else {
    stop("ref must be sequence name or index")
  }
  
  if (!ref_name %in% names(x)) {
    stop("Reference sequence not found in input DNAStringSet: ", ref_name)
  }
  
  input_width <- width(x)
  names(input_width) <- names(x)
  
  consensus_len <- input_width[ref_name]
  
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
  
  out <- data.frame(
    set = set_name,
    sequence = nm,
    ref_sequence = ref_name,
    match_len = as.integer(input_width[nm]),
    consensus_len = as.integer(consensus_len),
    n = n_insertions,
    k80 = k80[ref_i, ],
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      coverage = pmin(match_len / consensus_len, 1),
      missing_fraction = 1 - coverage,
      k80_norm = pmin(k80 / k80_cap, 1),
      te_decay = 1 - ((1 - k80_norm)^k80_weight * coverage^length_weight)
    )
  
  if (!include_ref) {
    out <- out %>% filter(sequence != ref_sequence)
  }
  
  out
}


dsuz_haya.decay <- get_te_decay_to_ref(
  dsuz_ht_fa, include_ref = T,
  ref = "Dsuz_ilop_1",
  set_name = "Dsuz_Hayabusa"
)

a6.decay <- get_te_decay_to_ref(
  c(reverseComplement(a6hopper.fa),
    longhop),
  ref = "A6_hopper_long",
  set_name = "Dmel_Nozomi (A6)",
  k80_cap = 0.5
)

bar.decay <- get_te_decay_to_ref(
  c(reverseComplement(barhop),
    longhop),
  ref = "A6_hopper_long",
  set_name = "Dmel_Nozomi (BL2969)"
)

dsuz_darj.decay <- get_te_decay_to_ref(
  dsuz_chopper.ref,
  ref = "Dsuz_hopster",
  set_name = "Dsuz_Darjeeling",
  include_ref = F,
  k80_weight = 1
)

# dsuz_coga_chop.decay <- get_te_decay_to_ref(
#   coga_chop.fa, include_ref = T,
#   ref = "Dsuz_chopper",
#   set_name = "Dsuz_Chopper (Bucharest)",
#   k80_weight = 1
# )

dm6_nozomi.decay <- get_te_decay_to_ref(
  hops.dm6.fa2,
  ref = "A6_hopper_long",
  set_name = "Dmel_Nozomi (dm6)"
)

ral426.decay <- get_te_decay_to_ref(
  c(hopsh.all.fa[grepl("RAL-426", names(hopsh.all.fa))],
    longhop),
  ref = "A6_hopper_long",
  set_name = "RAL-426_Hopper"
)

decay_df.1 <- rbind(
  dsuz_coga_chop.decay,
  # dsuz_interloper.decay,
  # a6.decay,
  dsuz_chopper.decay

)



ggplot(decay_df.1, aes(x = te_decay, weight = n, fill = set, color = set)) +
  geom_histogram(
    binwidth = 0.025,
    position = "identity",
    alpha = 0.35,
    linewidth = 0.25,
    boundary = 0
  ) +
  theme_classic() +
  labs(
    x = "TE decay",
    y = "Copy number per bin",
    fill = "TE",
    color = "TE"
  )+
  scale_y_continuous(
    breaks = scales::breaks_width(5)
  )

#### Decay for different TEs ####


decay_df.2 <- rbind(
  dsuz_haya.decay,
  dm6_nozomi.decay,
  dsuz_darj.decay
  
)

pal <- c(
  "Dmel_Nozomi (dm6)"    = "#2E7D32",  # same in both plots
  "Dsuz_Darjeeling"         = "#AB91E8",
  "Dsuz_Darjeeling (Bucharest)" = "#3B2168",
  "Dsuz_Hayabusa"      = "#E69F00",
  "Dmel_Nozomi (A6)"     = "#4E79A7",
  "Dmel_Nozomi (BL2969)" = "#D95F5F"
)

te_breaks <- c(
  "Dmel_Nozomi (dm6)",
  "Dsuz_Darjeeling",
  "Dsuz_Hayabusa"
)

te_labels <- c(
  "Dmel_Nozomi (dm6)" = "Nozomi (dm6)",
  "Dsuz_Darjeeling" = "Darjeeling (Dsuz)",
  "Dsuz_Hayabusa" = "Hayabusa (Dsuz)"
)


p_decay_tes <- ggplot(decay_df.2, aes(x = te_decay, weight = n, fill = set, color = set)) +
  geom_histogram(
    binwidth = 0.025,
    position = "identity",
    alpha = 0.5,
    linewidth = 0.25,
    boundary = 0
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
    x = "TE decay index",
    y = "Copy number per bin",
    fill = "TE",
    color = "TE"
  )

pdf("outputs/plots/te_decay_hli.pdf", width = 5, height = 3)
p_decay_tes
dev.off()

#### Decay for nozomi in dm6, a6 and bar ####
decay_df.3 <- rbind(
  dm6_nozomi.decay,
  bar.decay,
  a6.decay
  
)

te_breaks2 <- c(
  "Dmel_Nozomi (dm6)",
  "Dmel_Nozomi (BL2969)",
  "Dmel_Nozomi (A6)"
)

te_labels2 <- c(
  "Dmel_Nozomi (dm6)" = "dm6",
  "Dmel_Nozomi (BL2969)" = "BL2969",
  "Dmel_Nozomi (A6)" = "A6"
)


p_decay_tes2 <- ggplot(decay_df.3, aes(x = te_decay, weight = n, fill = set, color = set)) +
  geom_histogram(
    binwidth = 0.025,
    position = "identity",
    alpha = 0.5,
    linewidth = 0.25,
    boundary = 0
  ) +
  scale_fill_manual(
    name = "Genome",
    values = pal,
    breaks = te_breaks2,
    labels = te_labels2
  ) +
  scale_color_manual(
    name = "Genome",
    values = pal,
    breaks = te_breaks2,
    labels = te_labels2
  ) +
  theme_classic() +
  labs(
    x = "TE decay index",
    y = "Copy number per bin",
    fill = "TE",
    color = "TE"
  )+
  ggtitle("Nozomi decay in dm6, A6 and BL2969 genomes")

pdf("outputs/plots/te_decay_nozomi.pdf", width = 5, height = 3)
p_decay_tes2
dev.off()

