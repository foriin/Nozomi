library(dplyr)
library(tidyr)
library(ape)
library(ggplot2)
library(DECIPHER)
library(purrr)

## Keep only A6 short hopper insertions
shops.a6 <- shops.fa[grepl("^A6", names(shops.fa))]
shops.a6.gr <- GRanges(mcols(shops.a6)$loc)
shops.a6.gr$strainid <- names(shops.a6.gr)

shops.a6.aln <- AlignSeqs(shops.a6)
hopa6.names <- names(shops.a6.aln)
hopper_names <- names(hopper_a6_aln)

## Pairwise K80 divergence
hopper_dnabin <- as.DNAbin(as.matrix(shops.a6.aln))
k80_mat.a6 <- dist.dna(hopper_dnabin, model = "K80", pairwise.deletion = TRUE, as.matrix = TRUE)

k80_df.a6 <- as.data.frame(as.table(k80_mat.a6), stringsAsFactors = FALSE) %>%
  setNames(c("ins1", "ins2", "k80")) %>%
  filter(ins1 != ins2) %>%
  mutate(pair = map2_chr(ins1, ins2, ~ paste(sort(c(.x, .y)), collapse = "___"))) %>%
  distinct(pair, .keep_all = TRUE) %>%
  dplyr::select(-pair)

## Pull genomic coordinates for the same insertions
## Assumes teorig.te has columns: strainid, seqnames, start, end
shopsa6gr_df <- as.data.frame(shops.a6.gr) %>%
  filter(strainid %in% hopa6.names) %>%
  distinct(strainid, .keep_all = TRUE) %>%
  dplyr::select(strainid, seqnames, start, end)

## Join coordinates to each pair and compute genomic distances
pair_df <- k80_df.a6 %>%
  left_join(shopsa6gr_df, by = c("ins1" = "strainid")) %>%
  rename(chr1 = seqnames, start1 = start, end1 = end) %>%
  left_join(shopsa6gr_df, by = c("ins2" = "strainid")) %>%
  rename(chr2 = seqnames, start2 = start, end2 = end) %>%
  mutate(
    mid1 = (start1 + end1) / 2,
    mid2 = (start2 + end2) / 2,
    same_chr = chr1 == chr2,
    genomic_dist_mid = ifelse(same_chr, abs(mid1 - mid2), NA_real_),
    genomic_dist_edge = ifelse(same_chr, pmax(start1, start2) - pmin(end1, end2), NA_real_),
    genomic_dist_edge = ifelse(!is.na(genomic_dist_edge) & genomic_dist_edge < 0, 0, genomic_dist_edge)
  )

## Use only same-chromosome pairs
corr_df <- pair_df %>%
  filter(same_chr, !is.na(k80), !is.na(genomic_dist_mid), genomic_dist_mid > 0)

## Correlation test
cor_res <- cor.test(
  corr_df$k80,
  corr_df$genomic_dist_mid,
  method = "spearman"
)

print(cor_res)

## Label for plot
cor_lab <- paste0(
  "Spearman rho = ", signif(unname(cor_res$estimate), 3),
  ", P = ", format.pval(cor_res$p.value, digits = 3, eps = 1e-3)
)

## Plot midpoint-to-midpoint distance vs K80
x_lab <- min(corr_df$genomic_dist_mid + 1, na.rm = TRUE)
y_lab <- max(corr_df$k80, na.rm = TRUE)

p <- ggplot(corr_df, aes(x = genomic_dist_mid + 1, y = k80)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  geom_smooth(method = "lm", se = FALSE) +
  annotate(
    "text",
    x = x_lab,
    y = y_lab,
    label = cor_lab,
    hjust = 0,
    vjust = 1.2,
    size = 4
  ) +
  labs(
    x = "Genomic distance between Kodama insertions (bp, log10 scale)",
    y = "Kimura 80 divergence",
    title = "Sequence divergence vs genomic distance\namong A6 Kodama insertions"
  ) +
  theme_bw()

print(p)

pdf("outputs/plots/a6_kodama_k80_vs_dist.pdf", width = 6, height = 4)
p
dev.off()

pair_df %>%
  mutate(chr_group = ifelse(chr1 == chr2, "same_chr", "different_chr")) %>%
  filter(!is.na(k80), !is.na(chr_group)) %>%
  { 
    print(
      . %>%
        group_by(chr_group) %>%
        summarise(
          n = n(),
          mean_k80 = mean(k80),
          median_k80 = median(k80),
          .groups = "drop"
        )
    )
    
    print(
      wilcox.test(
        k80 ~ chr_group,
        data = .,
        alternative = "less"
      )
    )
    
    ggplot(., aes(chr_group, k80)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.15, alpha = 0.25) +
      theme_bw()
  }


pair_df %>%
  mutate(
    chr1_main = sub("[LR]$", "", chr1),
    chr2_main = sub("[LR]$", "", chr2),
    chr_group = ifelse(chr1_main == chr2_main, "same_chr", "different_chr")
  ) %>%
  filter(!is.na(k80), !is.na(chr_group)) %>%
  {
    print(
      . %>%
        group_by(chr_group) %>%
        summarise(
          n = n(),
          mean_k80 = mean(k80),
          median_k80 = median(k80),
          .groups = "drop"
        )
    )
    
    print(
      wilcox.test(
        k80 ~ chr_group,
        data = .,
        alternative = "less"
      )
    )
    
    ggplot(., aes(k80)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.15, alpha = 0.25) +
      theme_bw()
  }
