library(ggplot2)
library(dplyr)
library(tibble)

.seq_df <- function(chars, x_start, y, fill_group) {
  tibble(
    x = seq_along(chars) + x_start - 1,
    y = y,
    nt = chars,
    fill_group = fill_group
  )
}

plot_mh_deletion <- function(
    left_flank,
    left_mh,
    deleted_middle,
    right_mh,
    right_flank,
    top_label = "Original allele",
    bottom_label = "Deleted allele",
    deletion_label = "central deletion",
    family = "mono"
) {
  
  lf  <- strsplit(left_flank, "")[[1]]
  lmh <- strsplit(left_mh, "")[[1]]
  del <- strsplit(deleted_middle, "")[[1]]
  rmh <- strsplit(right_mh, "")[[1]]
  rf  <- strsplit(right_flank, "")[[1]]
  
  #### top row ####
  x0 <- 1
  
  top_lf  <- .seq_df(lf,  x0, 2, "flank")
  x1 <- max(top_lf$x) + 1
  
  top_lmh <- .seq_df(lmh, x1, 2, "left_mh")
  x2 <- max(top_lmh$x) + 1
  
  top_del <- .seq_df(del, x2, 2, "deleted")
  x3 <- if (length(del) > 0) max(top_del$x) + 1 else x2
  
  top_rmh <- .seq_df(rmh, x3, 2, "right_mh")
  x4 <- max(top_rmh$x) + 1
  
  top_rf  <- .seq_df(rf,  x4, 2, "flank")
  
  top_df <- bind_rows(top_lf, top_lmh, top_del, top_rmh, top_rf)
  
  #### bottom row ####
  bottom_total_width <- length(lf) + length(lmh) + length(rf)
  top_total_width <- max(top_df$x)
  
  bottom_x0 <- floor((top_total_width - bottom_total_width) / 2) + 1
  
  bot_lf <- .seq_df(lf, bottom_x0, 1, "flank")
  xb1 <- max(bot_lf$x) + 1
  
  bot_mh <- .seq_df(lmh, xb1, 1, "left_mh")
  xb2 <- max(bot_mh$x) + 1
  
  bot_rf <- .seq_df(rf, xb2, 1, "flank")
  
  bottom_df <- bind_rows(bot_lf, bot_mh, bot_rf)
  
  plot_df <- bind_rows(top_df, bottom_df)
  
  #### deleted region bracket ####
  if (length(del) > 0) {
    del_start <- min(top_del$x)
    del_end   <- max(top_rmh$x)
  } else {
    del_start <- min(top_rmh$x)
    del_end   <- max(top_rmh$x)
  }
  
  bot_mh_center <- mean(bot_mh$x)
  
  p <- ggplot(plot_df, aes(x = x, y = y)) +
    geom_tile(
      aes(fill = fill_group),
      width = 0.95,
      height = 0.72,
      color = "grey20",
      linewidth = 0.25
    ) +
    geom_text(
      aes(label = nt),
      family = family,
      size = 4
    ) +
    
    annotate(
      "segment",
      x = del_start - 0.45,
      xend = del_end + 0.45,
      y = 2.62,
      yend = 2.62,
      linewidth = 0.4
    ) +
    annotate(
      "segment",
      x = del_start - 0.45,
      xend = del_start - 0.45,
      y = 2.54,
      yend = 2.70,
      linewidth = 0.4
    ) +
    annotate(
      "segment",
      x = del_end + 0.45,
      xend = del_end + 0.45,
      y = 2.54,
      yend = 2.70,
      linewidth = 0.4
    ) +
    annotate(
      "text",
      x = mean(c(del_start, del_end)),
      y = 2.82,
      label = deletion_label,
      size = 3.5
    ) +
    
    annotate(
      "text",
      x = min(plot_df$x) - 1.2,
      y = 2,
      label = top_label,
      hjust = 1,
      size = 3.4
    ) +
    annotate(
      "text",
      x = min(plot_df$x) - 1.2,
      y = 1,
      label = bottom_label,
      hjust = 1,
      size = 3.4
    ) +
    
    annotate(
      "text",
      x = bot_mh_center,
      y = 0.56,
      label = "retained MH",
      size = 3
    ) +
    
    scale_fill_manual(
      values = c(
        flank    = "white",
        left_mh  = "grey70",
        deleted  = "grey92",
        right_mh = "grey55"
      ),
      breaks = c("flank", "left_mh", "deleted", "right_mh"),
      labels = c("Flank", "Left MH", "Deleted sequence", "Right MH"),
      name = NULL
    ) +
    coord_cartesian(clip = "off") +
    theme_void(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.margin = margin(20, 20, 20, 90)
    )
  
  p
}


pmd <- plot_mh_deletion(
  left_flank     = "AAAAACC",
  left_mh        = "ATTAG",
  deleted_middle = "AAAATT...AATCGT",
  right_mh       = "ATTTG",
  right_flank    = "TAGCTTCGC",
  deletion_label = "central deletion"
)



pdf("outputs/plots/hopper_central_del_scheme.pdf", width = 9, height = 3)
pmd
dev.off()
