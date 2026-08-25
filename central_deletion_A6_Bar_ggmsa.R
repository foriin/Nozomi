library(DECIPHER)
library(Biostrings)
library(ggplot2)
library(dplyr)

#### Load alignments ####
load('outputs/RData/bar_a6_unique_collapsed_aln.RData', verb = T)

## ------------------------------------------------------------
## Find the longest gap in selected sequence
## ------------------------------------------------------------

find_longest_gap <- function(aln, which_seq = length(aln), min_gap = 1000) {
  
  x <- strsplit(as.character(aln[[which_seq]]), "")[[1]]
  
  rr <- base::rle(x == "-")
  ends <- cumsum(rr$lengths)
  starts <- ends - rr$lengths + 1
  
  gap_tbl <- data.frame(
    start = starts,
    end = ends,
    gap_width = rr$lengths,
    is_gap = rr$values
  )
  
  gap_tbl <- gap_tbl[
    gap_tbl$is_gap & gap_tbl$gap_width >= min_gap,
  ]
  
  if (nrow(gap_tbl) == 0) {
    stop("No gap longer than min_gap found in selected sequence.")
  }
  
  gap_tbl <- gap_tbl[order(gap_tbl$gap_width, decreasing = TRUE), ]
  gap_tbl[1, , drop = FALSE]
}


## ------------------------------------------------------------
## Build plotting data frame for collapsed breakpoint view
## ------------------------------------------------------------

make_deletion_zoom_df <- function(
    aln,
    which_seq = length(aln),
    min_gap = 1000,
    flank = 35,
    deleted_edge = 15,
    gap_space = 8
) {
  
  mat <- as.matrix(aln)
  seq_names <- names(aln)
  nseq <- nrow(mat)
  aln_width <- ncol(mat)
  
  gap <- find_longest_gap(
    aln = aln,
    which_seq = which_seq,
    min_gap = min_gap
  )
  
  left_cols <- seq(
    from = max(1, gap$start - flank),
    to   = min(aln_width, gap$start + deleted_edge - 1)
  )
  
  right_cols <- seq(
    from = max(1, gap$end - deleted_edge + 1),
    to   = min(aln_width, gap$end + flank)
  )
  
  left_width <- length(left_cols)
  right_x_offset <- left_width + gap_space
  
  left_df <- expand.grid(
    seq_i = seq_len(nseq),
    local_i = seq_along(left_cols)
  )
  
  left_df$aln_col <- left_cols[left_df$local_i]
  left_df$x <- left_df$local_i
  left_df$side <- "left"
  
  right_df <- expand.grid(
    seq_i = seq_len(nseq),
    local_i = seq_along(right_cols)
  )
  
  right_df$aln_col <- right_cols[right_df$local_i]
  right_df$x <- right_x_offset + right_df$local_i
  right_df$side <- "right"
  
  df <- rbind(left_df, right_df)
  
  df$nt <- mat[cbind(df$seq_i, df$aln_col)]
  df$sequence <- seq_names[df$seq_i]
  
  df$y <- nseq - df$seq_i + 1
  
  df$region <- ifelse(
    df$aln_col >= gap$start & df$aln_col <= gap$end,
    "deleted_region_edge",
    "flank"
  )
  
  df$is_gap <- df$nt == "-"
  
  left_deleted_x <- df$x[
    df$side == "left" &
      df$region == "deleted_region_edge"
  ]
  
  right_deleted_x <- df$x[
    df$side == "right" &
      df$region == "deleted_region_edge"
  ]
  
  list(
    df = df,
    gap = gap,
    nseq = nseq,
    seq_names = seq_names,
    left_width = left_width,
    right_x_start = right_x_offset + 1,
    break_x = mean(c(left_width + 0.5, right_x_offset + 0.5)),
    left_deleted_range = range(left_deleted_x),
    right_deleted_range = range(right_deleted_x),
    x_max = max(df$x)
  )
}


## ------------------------------------------------------------
## Plot collapsed central deletion from real AlignSeqs alignment
## ------------------------------------------------------------

parse_alignment_names <- function(seq_names) {
  
  ok <- grepl("_w_[0-9]+_n_[0-9]+$", seq_names)
  
  id <- seq_names
  width <- rep(NA_character_, length(seq_names))
  n <- rep(NA_character_, length(seq_names))
  
  id[ok] <- sub("^(.*)_w_([0-9]+)_n_([0-9]+)$", "\\1", seq_names[ok])
  width[ok] <- sub("^(.*)_w_([0-9]+)_n_([0-9]+)$", "\\2", seq_names[ok])
  n[ok] <- sub("^(.*)_w_([0-9]+)_n_([0-9]+)$", "\\3", seq_names[ok])
  
  data.frame(
    sequence = seq_names,
    id = id,
    width = width,
    n = n,
    stringsAsFactors = FALSE
  )
}


plot_deletion_zoom <- function(
  aln,
  which_seq = length(aln),
  min_gap = 1000,
  flank = 35,
  deleted_edge = 15,
  gap_space = 8,
  deletion_label = NULL,
  base_size = 11,
  nt_text_size = 3.4,
  table_text_size = 3.2,
  family = "mono"
) {
  
  z <- make_deletion_zoom_df(
    aln = aln,
    which_seq = which_seq,
    min_gap = min_gap,
    flank = flank,
    deleted_edge = deleted_edge,
    gap_space = gap_space
  )
  
  df <- z$df
  df$nt <- toupper(df$nt)
  
  if (is.null(deletion_label)) {
    deletion_label <- paste0(z$gap$gap_width, " bp central deletion")
  }
  
  draw_df <- df[!df$is_gap, ]
  
  ## table on the left
  table_df <- parse_alignment_names(z$seq_names)
  table_df$seq_i <- seq_len(z$nseq)
  table_df$y <- z$nseq - table_df$seq_i + 1
  
  max_id_len <- max(nchar(c(table_df$id, "ID")), na.rm = TRUE)
  
  id_x <- -max_id_len * 0.75 - 7
  width_x <- -4.0
  n_x <- -1.0
  
  table_header_y <- z$nseq + 1.02
  
  p <- ggplot(draw_df, aes(x = x, y = y)) +
    
    ## faint background over shown deleted edges
    annotate(
      "rect",
      xmin = z$left_deleted_range[1] - 0.5,
      xmax = z$left_deleted_range[2] + 0.5,
      ymin = 0.5,
      ymax = z$nseq + 0.5,
      fill = "grey96",
      color = NA
    ) +
    annotate(
      "rect",
      xmin = z$right_deleted_range[1] - 0.5,
      xmax = z$right_deleted_range[2] + 0.5,
      ymin = 0.5,
      ymax = z$nseq + 0.5,
      fill = "grey96",
      color = NA
    ) +
    
    ## nucleotide boxes
    geom_tile(
      aes(fill = nt),
      width = 0.92,
      height = 0.72,
      color = "grey25",
      linewidth = 0.18
    ) +
    geom_text(
      aes(label = nt),
      family = family,
      size = nt_text_size
    ) +
    
    ## visual break
    annotate(
      "text",
      x = z$break_x,
      y = (z$nseq + 1) / 2,
      label = "//",
      size = 6,
      fontface = "bold"
    ) +
    
    ## deletion brackets
    annotate(
      "segment",
      x = z$left_deleted_range[1] - 0.45,
      xend = z$left_deleted_range[2] + 0.45,
      y = z$nseq + 0.78,
      yend = z$nseq + 0.78,
      linewidth = 0.35
    ) +
    annotate(
      "segment",
      x = z$right_deleted_range[1] - 0.45,
      xend = z$right_deleted_range[2] + 0.45,
      y = z$nseq + 0.78,
      yend = z$nseq + 0.78,
      linewidth = 0.35
    ) +
    annotate(
      "text",
      x = z$break_x,
      y = z$nseq + 1.02,
      label = deletion_label,
      size = 3.4
    ) +
    
    ## left-side table headers
    annotate(
      "text",
      x = id_x,
      y = table_header_y,
      label = "ID",
      hjust = 0,
      family = family,
      fontface = "bold",
      size = table_text_size
    ) +
    annotate(
      "text",
      x = width_x,
      y = table_header_y,
      label = "width",
      hjust = 1,
      family = family,
      fontface = "bold",
      size = table_text_size
    ) +
    annotate(
      "text",
      x = n_x,
      y = table_header_y,
      label = "n",
      hjust = 1,
      family = family,
      fontface = "bold",
      size = table_text_size
    ) +
    
    ## left-side table rows
    geom_text(
      data = table_df,
      aes(x = id_x, y = y, label = id),
      inherit.aes = FALSE,
      hjust = 0,
      family = family,
      size = table_text_size
    ) +
    geom_text(
      data = table_df,
      aes(x = width_x, y = y, label = width),
      inherit.aes = FALSE,
      hjust = 1,
      family = family,
      size = table_text_size
    ) +
    geom_text(
      data = table_df,
      aes(x = n_x, y = y, label = n),
      inherit.aes = FALSE,
      hjust = 1,
      family = family,
      size = table_text_size
    ) +
    
    ## subtle separator between table and alignment
    annotate(
      "segment",
      x = 0.15,
      xend = 0.15,
      y = 0.45,
      yend = z$nseq + 0.55,
      linewidth = 0.25,
      color = "grey45"
    ) +
    
    scale_x_continuous(
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.03, 0.12))
    ) +
    scale_fill_manual(
      values = c(
        A = "#6AFF6A",
        C = "#59B9FF",
        G = "#FFB454",
        T = "#FF7A83",
        N = "grey70"
      ),
      na.value = "grey90",
      name = NULL
    ) +
    coord_cartesian(
      xlim = c(id_x - 0.5, z$x_max + 0.5),
      ylim = c(0.5, z$nseq + 1.25),
      clip = "off"
    ) +
    theme_classic(base_size = base_size) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      legend.position = "none",
      plot.margin = margin(15, 20, 15, 10)
    )
  
  p
}

## ------------------------------------------------------------
## Run on your AlignSeqs alignment
## ------------------------------------------------------------

p_a6 <- plot_deletion_zoom(
  a6.uq.aln.sorted,
  which_seq = 1,
  min_gap = 1000,
  flank = 35,
  deleted_edge = 15,
  gap_space = 9,
  deletion_label = "1380 bp central deletion"
)

print(p_a6)

p_bar <- plot_deletion_zoom(
  bar.uq.aln.sorted,
  which_seq = 1,
  min_gap = 1000,
  flank = 35,
  deleted_edge = 15,
  gap_space = 9,
  deletion_label = "1380 bp central deletion"
)

print(p_bar)


## ------------------------------------------------------------
## Save
## ------------------------------------------------------------

pdf("outputs/plots/a6_central_deletion_alignment_zoom_colored.pdf", width = 15, height = 10)
p_a6
dev.off()

pdf("outputs/plots/bar_central_deletion_alignment_zoom_colored.pdf", width = 15, height = 8)
p_bar
dev.off()

ggsave(
  "central_deletion_alignment_zoom_colored.svg",
  p,
  width = 9,
  height = 3.5
)