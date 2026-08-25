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
