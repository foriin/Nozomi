library(Biostrings)
library(ggplot2)

fa <- readDNAStringSet("inputs/fasta/nozomi_kodama.fa")
stopifnot(length(fa) == 2)

s1 <- fa[[1]]
s2 <- fa[[2]]

aln <- pairwiseAlignment(
  pattern = s1,
  subject = s2,
  type = "global",
  substitutionMatrix = nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE),
  gapOpening = 8,
  gapExtension = 1
)

a1 <- strsplit(as.character(alignedPattern(aln)), "")[[1]]
a2 <- strsplit(as.character(alignedSubject(aln)), "")[[1]]

p1 <- cumsum(a1 != "-")
p2 <- cumsum(a2 != "-")

df <- data.frame(
  x = p1,
  y = p2,
  b1 = a1,
  b2 = a2
)

df <- subset(df, b1 != "-" & b2 != "-")
df$match <- df$b1 == df$b2

p <- ggplot(df, aes(x, y)) +
  geom_point(data = subset(df, !match), size = 0.15, alpha = 0.15) +
  geom_point(data = subset(df, match), size = 0.2) +
  coord_fixed() +
  theme_classic(base_size = 14) +
  labs(
    x = "Nozomi",
    y = "Kodama"
  )+
  ggtitle("Pairwise alignment Nozomi vs Kodama")


pdf("outputs/plots/nozo_koda_pairwise_dotplot.pdf", width = 5, height = 5)
p
dev.off()
