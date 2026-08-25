library(ggplot2)

# -------------------------------------------------------------------------
# Input/output
# -------------------------------------------------------------------------

input_file <- "~/Work/projects/DrosoTE/Hopper_A6/TE_Blast_runner/pelement/copy_number.tsv"
# output_file <- "outputs/plots/Transib1_CN.pdf"
element_name <- "P-element"

cn <- read.delim(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c("genome", "group", "panel", "class", "n")

if (!all(required_columns %in% names(cn))) {
  stop(
    "Input file must contain these columns: ",
    paste(required_columns, collapse = ", ")
  )
}

cn$n <- as.numeric(cn$n)

# -------------------------------------------------------------------------
# Panel order and labels
# -------------------------------------------------------------------------

group_order <- c("Reference", "DSPR", "DGRP", "DrosEU")

cn$group_plot <- factor(
  cn$group,
  levels = group_order,
  labels = c("ref", "DSPR (1950-1960s)", "DGRP (2002)", "DrosEU (2010s)")
)

cn$panel_plot <- factor(
  cn$panel,
  levels = c(
    "All reduced copies",
    "Reduced full-length only"
  ),
  labels = c(
    "All copies",
    "Full-length copies"
  )
)

# -------------------------------------------------------------------------
# Genome order
# -------------------------------------------------------------------------

# Keep genomes grouped from left to right:
# ref -> DSPR -> DGRP -> DrosEU
genome_order <- unlist(
  lapply(group_order, function(grp) {
    sort(unique(cn$genome[cn$group == grp]))
  }),
  use.names = FALSE
)

cn$genome_plot <- factor(
  cn$genome,
  levels = genome_order
)

# Clean labels without changing the underlying genome identifiers
genome_labels <- genome_order

genome_labels <- sub(
  "\\.(fasta|fa|fna)$",
  "",
  genome_labels,
  ignore.case = TRUE
)

genome_labels <- sub(
  "_pacbio$",
  "",
  genome_labels,
  ignore.case = TRUE
)

genome_labels <- sub(
  "^dmel-all-chromosome-r[0-9.]+$",
  "dm6",
  genome_labels
)

# Force every genome in the reference panel to display as dm6
genome_labels[genome_order %in% cn$genome[cn$group == "Reference"]] <- "dm6"

names(genome_labels) <- genome_order

# -------------------------------------------------------------------------
# Copy classes and colours
# -------------------------------------------------------------------------

cn$class_plot <- factor(
  cn$class,
  levels = c(
    "recent/full",
    "older/full",
    "recent/partial",
    "older/partial"
  ),
  labels = c(
    "Low div, full-length",
    "High div, full-length",
    "Low div, partial",
    "High div, partial"
  )
)

copy_colours <- c(
  "Low div, full-length" = "#2E7D32",
  "High div, full-length"  = "#78A66F",
  "Low div, partial"     = "#D6A443",
  "High div, partial"     = "#B8AAA0"
)



# -------------------------------------------------------------------------
# Plot
# -------------------------------------------------------------------------

p <- ggplot(
  cn,
  aes(
    x = genome_plot,
    y = n,
    fill = class_plot
  )
) +
  geom_col(
    width = 0.82,
    colour = "white",
    linewidth = 0.18
  ) +
  
  facet_grid(
    rows = vars(panel_plot),
    cols = vars(group_plot),
    scales = "free",
    space = "free_x",
    switch = "y"
  ) +
  
  scale_fill_manual(
    values = copy_colours,
    breaks = c(
      "Low div, full-length",
      "High div, full-length",
      "Low div, partial",
      "High div, partial"
    ),
    drop = FALSE
  ) +
  
  scale_x_discrete(
    labels = genome_labels,
    expand = expansion(add = c(0.75, 0.75))
  ) +
  
  scale_y_continuous(
    breaks = scales::breaks_pretty(n = 5),
    expand = expansion(mult = c(0, 0.06))
  ) +
  
  labs(
    title = element_name,
    x = NULL,
    y = "Copy number",
    fill = NULL
  ) +
  
  guides(
    fill = guide_legend(
      ncol = 1,
      byrow = TRUE
    )
  ) +
  
  theme_classic(base_size = 9) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13,
      hjust = 0,
      margin = margin(b = 8)
    ),
    
    axis.text.x = element_text(
      angle = 55,
      hjust = 1,
      vjust = 1,
      size = 6.8,
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 8,
      colour = "black"
    ),
    
    axis.title.y = element_text(
      size = 10,
      margin = margin(r = 8)
    ),
    
    strip.background = element_rect(
      fill = "grey94",
      colour = "black",
      linewidth = 0.4
    ),
    
    strip.text.x = element_text(
      face = "bold",
      size = 9
    ),
    
    strip.text.y.left = element_text(
      face = "bold",
      size = 9,
      angle = 90
    ),
    
    strip.placement = "outside",
    
    # Box around every facet panel
    panel.border = element_rect(
      fill = NA,
      colour = "black",
      linewidth = 0.45
    ),
    
    panel.spacing.x = grid::unit(0.7, "lines"),
    panel.spacing.y = grid::unit(0.7, "lines"),
    
    panel.grid.major.y = element_line(
      colour = "grey88",
      linewidth = 0.3
    ),
    
    panel.grid.minor = element_blank(),
    
    # panel.border now provides the axes/frame
    axis.line = element_blank(),
    
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.3
    ),
    
    legend.position = "right",
    legend.direction = "vertical",
    legend.justification = "center",
    
    legend.key.width = grid::unit(0.45, "cm"),
    legend.key.height = grid::unit(0.45, "cm"),
    
    legend.text = element_text(
      size = 8
    ),
    
    plot.margin = margin(
      t = 8,
      r = 8,
      b = 8,
      l = 8
    )
  )


p

pdf("outputs/plots/p_element_POPCN.pdf", width = 9, height = 4)
p
dev.off()


