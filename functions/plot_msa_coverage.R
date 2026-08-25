plot_msa_coverage <- function(
    alignment,
    ref_name,
    min_deletion_width = 20,
    min_deleted_fraction = 0.9,
    max_interruption_width = 5,
    include_terminal_gaps = FALSE,
    x_marks = NULL,
    fill_color = "#4E79A7",
    fill_alpha = 0.75,
    line_color = "black",
    line_width = 0.35,
    title = NULL
) {
  
  # -----------------------------------------------------------------------
  # Validate input
  # -----------------------------------------------------------------------
  
  aln_char <- as.character(alignment)
  
  if (is.null(names(aln_char))) {
    stop("The aligned sequences must have names.")
  }
  
  if (!ref_name %in% names(aln_char)) {
    stop("Reference sequence '", ref_name, "' was not found.")
  }
  
  if (sum(names(aln_char) == ref_name) != 1) {
    stop("ref_name must identify exactly one sequence.")
  }
  
  if (length(aln_char) < 2) {
    stop("The alignment must contain a reference and at least one query.")
  }
  
  if (length(unique(nchar(aln_char))) != 1) {
    stop("All sequences must have the same aligned length.")
  }
  
  if (
    min_deleted_fraction < 0 ||
    min_deleted_fraction > 1
  ) {
    stop("min_deleted_fraction must be between 0 and 1.")
  }
  
  # -----------------------------------------------------------------------
  # Convert alignment to a character matrix
  # -----------------------------------------------------------------------
  
  aln_mat <- do.call(
    rbind,
    strsplit(aln_char, split = "", fixed = TRUE)
  )
  
  rownames(aln_mat) <- names(aln_char)
  
  ref <- aln_mat[ref_name, ]
  
  # -----------------------------------------------------------------------
  # Convert alignment columns to reference coordinates
  #
  # Columns in which the reference has a gap are removed.
  # -----------------------------------------------------------------------
  
  ref_has_base <- ref != "-" & ref != "."
  
  ref_pos <- cumsum(ref_has_base)
  ref_pos[!ref_has_base] <- NA_integer_
  
  keep_cols <- !is.na(ref_pos)
  
  aln_ref <- aln_mat[, keep_cols, drop = FALSE]
  ref_pos_kept <- ref_pos[keep_cols]
  
  query_mat <- aln_ref[
    rownames(aln_ref) != ref_name,
    ,
    drop = FALSE
  ]
  
  n_queries <- nrow(query_mat)
  ref_len <- max(ref_pos_kept)
  
  # -----------------------------------------------------------------------
  # Coverage
  # -----------------------------------------------------------------------
  
  query_has_base <- (
    query_mat != "-" &
      query_mat != "."
  )
  
  coverage <- colSums(query_has_base)
  frac_covered <- coverage / n_queries
  
  coverage_df <- tibble::tibble(
    ref = ref_name,
    pos = ref_pos_kept,
    coverage = coverage,
    n_sequences = n_queries,
    frac_covered = frac_covered,
    deleted_fraction = 1 - frac_covered
  )
  
  # -----------------------------------------------------------------------
  # Helper: convert a logical vector into contiguous runs
  # -----------------------------------------------------------------------
  
  logical_runs <- function(x) {
    
    r <- rle(x)
    
    end_index <- cumsum(r$lengths)
    start_index <- end_index - r$lengths + 1L
    
    tibble::tibble(
      value = r$values,
      start_index = start_index,
      end_index = end_index,
      width = r$lengths
    )
  }
  
  # -----------------------------------------------------------------------
  # Helper: join deletion stretches separated by very short interruptions
  #
  # This prevents one or two isolated covered positions from splitting one
  # biologically continuous deletion into multiple regions.
  # -----------------------------------------------------------------------
  
  bridge_interruptions <- function(mask, max_width) {
    
    if (max_width <= 0 || length(mask) < 3) {
      return(mask)
    }
    
    runs <- logical_runs(mask)
    
    for (i in seq_len(nrow(runs))) {
      
      is_short_false_run <- (
        !runs$value[i] &&
          runs$width[i] <= max_width
      )
      
      has_flanking_runs <- (
        i > 1 &&
          i < nrow(runs)
      )
      
      flanked_by_deletions <- (
        has_flanking_runs &&
          runs$value[i - 1] &&
          runs$value[i + 1]
      )
      
      if (is_short_false_run && flanked_by_deletions) {
        mask[
          runs$start_index[i]:
            runs$end_index[i]
        ] <- TRUE
      }
    }
    
    mask
  }
  
  # -----------------------------------------------------------------------
  # Identify internal consensus deletions
  #
  # A position is considered deleted when at least min_deleted_fraction
  # of query sequences contain a gap at that reference position.
  # -----------------------------------------------------------------------
  
  deletion_mask <- (
    coverage_df$deleted_fraction >=
      min_deleted_fraction
  )
  
  deletion_mask <- bridge_interruptions(
    deletion_mask,
    max_width = max_interruption_width
  )
  
  deletion_runs <- logical_runs(deletion_mask) |>
    dplyr::filter(
      value,
      width >= min_deletion_width
    )
  
  if (!include_terminal_gaps) {
    deletion_runs <- deletion_runs |>
      dplyr::filter(
        start_index > 1,
        end_index < nrow(coverage_df)
      )
  }
  
  if (nrow(deletion_runs) > 0) {
    
    consensus_deletions <- tibble::tibble(
      start_index = deletion_runs$start_index,
      end_index = deletion_runs$end_index,
      start = coverage_df$pos[
        deletion_runs$start_index
      ],
      end = coverage_df$pos[
        deletion_runs$end_index
      ]
    ) |>
      dplyr::rowwise() |>
      dplyr::mutate(
        width = end - start + 1L,
        mean_coverage = mean(
          coverage_df$coverage[
            start_index:end_index
          ]
        ),
        mean_covered_fraction = mean(
          coverage_df$frac_covered[
            start_index:end_index
          ]
        ),
        mean_deleted_fraction = mean(
          coverage_df$deleted_fraction[
            start_index:end_index
          ]
        )
      ) |>
      dplyr::ungroup() |>
      dplyr::select(
        start,
        end,
        width,
        mean_coverage,
        mean_covered_fraction,
        mean_deleted_fraction
      )
    
  } else {
    
    consensus_deletions <- tibble::tibble(
      start = integer(),
      end = integer(),
      width = integer(),
      mean_coverage = numeric(),
      mean_covered_fraction = numeric(),
      mean_deleted_fraction = numeric()
    )
  }
  
  # -----------------------------------------------------------------------
  # X-axis marks
  #
  # By default, use reference boundaries and deletion boundaries.
  # -----------------------------------------------------------------------
  
  if (is.null(x_marks)) {
    
    x_marks <- c(
      1L,
      consensus_deletions$start,
      consensus_deletions$end,
      ref_len
    )
    
  } else {
    
    x_marks <- c(
      1L,
      x_marks,
      ref_len
    )
  }
  
  x_marks <- sort(
    unique(
      as.integer(x_marks)
    )
  )
  
  x_marks <- x_marks[
    x_marks >= 1 &
      x_marks <= ref_len
  ]
  
  # -----------------------------------------------------------------------
  # Plot
  # -----------------------------------------------------------------------
  
  if (is.null(title)) {
    title <- paste0(
      "MSA coverage of query sequences\nalong ",
      ref_name
    )
  }
  
  y_max <- max(
    coverage_df$coverage,
    na.rm = TRUE
  )
  
  p <- ggplot2::ggplot(
    coverage_df,
    ggplot2::aes(
      x = pos,
      y = coverage
    )
  ) +
    
    # Fill only the area represented by sequence coverage.
    # Zero-coverage deletion regions remain white.
    ggplot2::geom_area(
      fill = fill_color,
      alpha = fill_alpha,
      colour = NA
    ) +
    
    ggplot2::geom_line(
      colour = line_color,
      linewidth = line_width
    ) +
    
    ggplot2::scale_x_continuous(
      breaks = x_marks,
      labels = as.character(x_marks),
      minor_breaks = NULL,
      expand = ggplot2::expansion(
        mult = c(0, 0)
      )
    ) +
    
    ggplot2::scale_y_continuous(
      limits = c(0, y_max),
      expand = ggplot2::expansion(
        mult = c(0, 0.04)
      ),
      breaks = scales::pretty_breaks(
        n = 5
      )
    ) +
    
    ggplot2::theme_classic() +
    
    ggplot2::theme(
      axis.line.x = ggplot2::element_line(
        linewidth = 0.35
      ),
      axis.line.y = ggplot2::element_line(
        linewidth = 0.35
      )
    ) +
    
    ggplot2::labs(
      x = "Position along reference, nt",
      y = "Number of sequences with nucleotide present",
      title = title,
      subtitle = paste0(
        "n query sequences = ",
        n_queries
      )
    )
  
  # -----------------------------------------------------------------------
  # Return plot and called deletion coordinates
  # -----------------------------------------------------------------------
  
  list(
    plot = p,
    coverage = coverage_df,
    consensus_deletions = consensus_deletions,
    x_marks = x_marks
  )
}