library(dplyr)
library(ggplot2)
library(ggrepel)

plot_df <- df_long2 %>%
  # filter(name %in% c("ROO", "HOPPER", "DOC", "COPIA", "F-ELEMENT", "HOBO")) %>%
  mutate(genome = factor(genome, levels = unique(genome)))

lab_df <- plot_df %>%
  group_by(name) %>%
  slice_max(value, n = 1, with_ties = FALSE) %>%
  ungroup()

ggplot(plot_df, aes(x = genome, y = value, col = name, group = name)) +
  geom_line(aes(linewidth = name %in% c("COPIA", "DOC", "HOPPER"),
                alpha = name %in% c("COPIA", "DOC", "HOPPER"))) +
  geom_text_repel(
    data = lab_df,
    aes(label = name),
    direction = "y",
    hjust = 0,
    nudge_x = 0.3,
    size = 2.5,
    show.legend = FALSE
  ) +
  scale_linewidth_manual(values = c(`TRUE` = 1.3, `FALSE` = 0.5), guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.5), guide = "none")+
  coord_cartesian(clip = "off") +
  theme_bw() +
  ylab("Number of isolated insertions")+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(5.5, 40, 5.5, 5.5)
  )
