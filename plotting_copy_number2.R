library(dplyr)
library(ggplot2)
library(ggrepel)

main_tes <- c("HOPPER", "COPIA", "DOC", "ROO", "HOBO")

pal <- c(
  HOPPER = "#D55E99",
  COPIA  = "#0072B2",
  DOC    = "#009E73",
  ROO    = "#FAAF00",
  HOBO   = "#7B61A8",
  Other  = "grey80"
)

plot_df2 <- plot_df %>%
  mutate(
    tier = case_when(
      name %in% main_tes ~ name,
      TRUE ~ "Other"
    ),
    tier = factor(
      tier,
      levels = c("Other", "COPIA", "DOC", "ROO", "HOBO", "HOPPER")
    )
  )

lab_df2 <- lab_df %>%
  mutate(
    tier = case_when(
      name %in% main_tes ~ name,
      TRUE ~ "Other"
    ),
    tier = factor(
      tier,
      levels = c("Other", "COPIA", "DOC", "ROO", "HOBO", "HOPPER")
    )
  )

max_point_df <- plot_df2 %>%
  filter(name %in% c("HOPPER", "DOC", "COPIA", "HOBO")) %>%
  group_by(name) %>%
  slice_max(order_by = value, n = 1, with_ties = FALSE) %>%
  ungroup()

pcn_line <- ggplot() +
  geom_line(
    data = plot_df2 %>% filter(tier == "Other"),
    aes(x = genome, y = value, group = name),
    color = unname(pal["Other"]),
    linewidth = 0.22,
    alpha = 0.35
  ) +
  geom_line(
    data = plot_df2 %>% filter(tier %in% c("COPIA", "DOC", "ROO", "HOBO")),
    aes(x = genome, y = value, color = tier, group = name),
    linewidth = 0.6,
    alpha = 0.95
  ) +
  geom_line(
    data = plot_df2 %>% filter(tier == "HOPPER"),
    aes(x = genome, y = value, color = tier, group = name),
    linewidth = 0.6,
    alpha = 0.95
  ) +
  geom_point(
    data = max_point_df,
    aes(x = genome, y = value, color = tier),
    size = 1.9,
    alpha = 0.95
  ) +
  geom_text_repel(
    data = lab_df2 %>% filter(tier %in% c("COPIA", "DOC", "ROO", "HOBO")),
    aes(x = genome, y = value, label = name, color = tier),
    direction = "y",
    hjust = 0,
    nudge_x = 0.12,
    size = 2.7,
    box.padding = 0.15,
    point.padding = 0.1,
    segment.size = 0.2,
    segment.alpha = 0.45,
    min.segment.length = 0,
    force = 0.6,
    show.legend = FALSE
  ) +
  geom_label_repel(
    data = lab_df2 %>% filter(tier == "HOPPER"),
    aes(x = genome, y = value, label = "Hopper", color = tier),
    direction = "y",
    hjust = 0,
    nudge_x = 0.18,
    size = 3.2,
    fontface = "bold",
    label.size = 0.25,
    fill = "white",
    box.padding = 0.2,
    point.padding = 0.12,
    segment.size = 0.22,
    segment.alpha = 0.6,
    min.segment.length = 0,
    force = 0.7,
    show.legend = FALSE
  ) +
  scale_color_manual(values = pal, guide = "none") +
  scale_x_discrete(expand = expansion(add = c(0.2, 0.9))) +
  coord_cartesian(clip = "off") +
  labs(
    x = NULL,
    y = "Number of isolated insertions"
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.grid.major.x = element_line(
      color = "grey90",
      linewidth = 0.3
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(margin = margin(r = 8)),
    
    axis.line = element_line(
      linewidth = 0.35,
      color = "black"
    ),
    axis.ticks = element_line(
      linewidth = 0.3,
      color = "black"
    ),
    
    plot.margin = margin(5.5, 45, 5.5, 5.5)
  )

pdf("outputs/plots/dspr_TE_CN_lineplot.pdf", width = 5, height = 3)
pcn_line
dev.off()
