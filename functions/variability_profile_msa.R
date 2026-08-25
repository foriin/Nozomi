# variability_profile_msa.R
#
# Usage (from another R script):
#   source("variability_profile_msa.R")
#   p <- plot_msa_variability(
#     genome_fa = "fasta/WGLKD.ONT.genome.flye_ragtag.fasta",
#     gr = copia.ins.intact,
#     name_col = "Name",
#     msa_method = "ClustalOmega"
#   )
#   print(p)
#
# Requirements:
#   Bioc: Rsamtools, Biostrings, GenomicRanges, msa
#   CRAN: ggplot2

suppressPackageStartupMessages({
  library(Rsamtools)
  library(Biostrings)
  library(GenomicRanges)
  library(msa)
  library(ggplot2)
})

# -------------------------
# Core helpers
# -------------------------

get_seqs_from_fasta <- function(genome_fa, gr, name_col = NULL) {
  stopifnot(inherits(gr, "GRanges"))
  if (!file.exists(genome_fa)) stop("FASTA not found: ", genome_fa, call. = FALSE)
  
  # Ensure FASTA is indexed
  if (!file.exists(paste0(genome_fa, ".fai"))) {
    indexFa(genome_fa)
  }
  
  fa <- FaFile(genome_fa)
  open(fa)
  on.exit(close(fa), add = TRUE)
  
  # Extract sequences (strand is respected by getSeq)
  seqs <- getSeq(fa, gr)
  
  # Names
  if (!is.null(name_col) && (name_col %in% colnames(mcols(gr)))) {
    nm <- as.character(mcols(gr)[[name_col]])
  } else if (!is.null(names(gr)) && all(names(gr) != "")) {
    nm <- names(gr)
  } else {
    nm <- paste0("seq_", seq_along(gr))
  }
  nm[is.na(nm) | nm == ""] <- paste0("seq_", which(is.na(nm) | nm == ""))
  nm <- make.unique(nm)
  
  names(seqs) <- nm
  seqs
}

aln_to_matrix <- function(aln) {
  if (inherits(aln, "MsaDNAMultipleAlignment")) {
    s <- as.character(unmasked(aln))
  } else if (inherits(aln, "DNAStringSet")) {
    s <- as.character(aln)
  } else {
    stop("Unsupported alignment class: ", paste(class(aln), collapse = ", "), call. = FALSE)
  }
  
  split_chars <- strsplit(s, split = "", fixed = TRUE)
  maxlen <- max(lengths(split_chars))
  
  mat <- do.call(rbind, lapply(split_chars, function(x) {
    if (length(x) < maxlen) x <- c(x, rep("-", maxlen - length(x)))
    x
  }))
  
  rownames(mat) <- names(s)
  mat
}

col_entropy <- function(col, alphabet = c("A", "C", "G", "T"), ignore = c("-", "N")) {
  col <- toupper(col)
  col <- col[!col %in% ignore]
  if (length(col) == 0) return(NA_real_)
  tab <- table(factor(col, levels = alphabet))
  p <- tab / sum(tab)
  p <- p[p > 0]
  -sum(p * log2(p))
}

col_mismatch_to_consensus <- function(col, ignore = c("-", "N")) {
  col <- toupper(col)
  keep <- !col %in% ignore
  x <- col[keep]
  if (length(x) == 0) return(NA_real_)
  cons <- names(which.max(table(x)))
  mean(x != cons)
}

col_gap_fraction <- function(col) {
  mean(col == "-")
}

variability_profile <- function(mat) {
  stopifnot(is.matrix(mat))
  L <- ncol(mat)
  
  entropy <- vapply(seq_len(L), function(i) col_entropy(mat[, i]), numeric(1))
  mismatch <- vapply(seq_len(L), function(i) col_mismatch_to_consensus(mat[, i]), numeric(1))
  gapfrac <- vapply(seq_len(L), function(i) col_gap_fraction(mat[, i]), numeric(1))
  
  data.frame(
    pos = seq_len(L),
    entropy = entropy,
    mismatch = mismatch,
    gap_fraction = gapfrac
  )
}

# ---- replace plot_profile() in variability_profile_msa.R with this version ----
plot_profile <- function(prof, anno = NULL,
                         track_height_frac = 0.20, track_gap_frac = 0.05,
                         alpha = 0.25, label_size = 3,
                         title = "Per-position variability profile",
                         subtitle = "Solid: entropy, dashed: mismatch to consensus, dotted: gap fraction") {
  stopifnot(all(c("pos", "entropy", "mismatch", "gap_fraction") %in% colnames(prof)))
  
  # base plot (keeps your styling)
  p <- ggplot(prof, aes(x = pos)) +
    geom_line(aes(y = entropy), linewidth = 0.6, col = "red", na.rm = TRUE) +
    geom_line(aes(y = mismatch), linewidth = 0.6, linetype = 2, col = "green", na.rm = TRUE) +
    geom_line(aes(y = gap_fraction), linewidth = 0.6, linetype = 3, na.rm = TRUE) +
    labs(x = "Alignment position", y = "Value", title = title, subtitle = subtitle) +
    theme_bw(base_size = 12)
  
  # annotation is optional
  if (is.null(anno) || nrow(anno) == 0) return(p)
  
  req <- c("feature", "start", "end")
  if (!all(req %in% colnames(anno))) {
    stop("anno must have columns: feature, start, end", call. = FALSE)
  }
  
  # y-range of plotted signals
  yr <- range(c(prof$entropy, prof$mismatch, prof$gap_fraction), na.rm = TRUE)
  d <- diff(yr)
  if (!is.finite(d) || d == 0) d <- 1
  
  track_h <- d * track_height_frac
  track_gap <- d * track_gap_frac
  
  y1 <- yr[1] - track_gap
  y0 <- y1 - track_h
  
  anno2 <- anno %>%
    mutate(
      start = as.numeric(start),
      end = as.numeric(end),
      xmid = (start + end) / 2
    )
  
  p +
    geom_rect(
      data = anno2,
      aes(xmin = start, xmax = end, ymin = y0, ymax = y1),
      inherit.aes = FALSE,
      alpha = alpha
    ) +
    geom_text(
      data = anno2,
      aes(x = xmid, y = (y0 + y1) / 2, label = feature),
      inherit.aes = FALSE,
      size = label_size
    ) +
    scale_y_continuous(
      limits = c(y0, yr[2]),
      expand = expansion(mult = c(0, 0.02))
    )
}



# -------------------------
# Main function
# -------------------------

# Returns a ggplot object. Optionally writes a PDF/PNG if out_file is given.
plot_msa_variability <- function(genome_fa,
                                 gr,
                                 name_col = "Name",
                                 msa_method = c("ClustalOmega", "Muscle", "ClustalW"),
                                 anno = NULL,
                                 out_file = NULL,
                                 width = 9,
                                 height = 4) {
  msa_method <- match.arg(msa_method)
  
  # 1) sequences
  seqs <- get_seqs_from_fasta(genome_fa, gr, name_col = name_col)
  
  # 2) MSA
  aln <- msa(seqs, method = msa_method)
  
  # 3) matrix
  mat <- aln_to_matrix(aln)
  
  # keep rownames aligned with sequence names (already set by get_seqs_from_fasta)
  # but enforce once in case msa drops names
  if (!is.null(names(seqs))) rownames(mat) <- names(seqs)
  
  # 4) variability
  prof <- variability_profile(mat)
  
  # 5) plot
  p <- plot_profile(
    prof,
    anno = anno,
    title = paste0("Per-position variability (", msa_method, "), n=", length(seqs))
  )
  
  if (!is.null(out_file)) {
    ggsave(out_file, p, width = width, height = height)
  }
  
  p
}
