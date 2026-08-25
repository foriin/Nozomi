#### Load packages ####

library(Biostrings)
library(DECIPHER)
library(ggplot2)


#### Check input and metadata ####

stopifnot(inherits(bwhopdm6, "DNAStringSet"))
stopifnot(!is.null(names(bwhopdm6)))
stopifnot(!anyDuplicated(names(bwhopdm6)))
stopifnot(all(c("genome", "count") %in% colnames(mcols(bwhopdm6))))

bwhopdm6 <- bwhopdm6[(width(bwhopdm6) > 1300) & (width(bwhopdm6) < 1700)]

genome <- as.character(mcols(bwhopdm6)$genome)
count <- as.numeric(mcols(bwhopdm6)$count)

stopifnot(all(genome %in% c("A6", "BL2969", "dm6")))
stopifnot(all(!is.na(count[genome %in% c("A6", "BL2969")])))
stopifnot(all(count[genome %in% c("A6", "BL2969")] > 0))


#### Align sequences and calculate distances ####

bwhopdm6.aln <- AlignSeqs(bwhopdm6)

dist.mat <- stringDist(
  bwhopdm6.aln,
  method = "hamming"
)

dist.mat <- as.dist(
  dist.mat / width(bwhopdm6.aln)[1]
)


#### Perform PCoA ####

pcoa <- cmdscale(
  dist.mat,
  k = 2,
  eig = TRUE,
  add = TRUE
)

positive.eig <- pcoa$eig[pcoa$eig > 0]
var.explained <- positive.eig / sum(positive.eig)


#### Prepare plotting data ####

pcoa.df <- data.frame(
  name = names(bwhopdm6.aln),
  PCoA1 = pcoa$points[, 1],
  PCoA2 = pcoa$points[, 2],
  genome = factor(
    genome,
    levels = c("dm6", "A6", "BL2969")
  ),
  count = count
)

xlab <- paste0(
  "PCoA1 (",
  round(var.explained[1] * 100, 1),
  "%)"
)

ylab <- paste0(
  "PCoA2 (",
  round(var.explained[2] * 100, 1),
  "%)"
)

genome.colours <- c(
  "dm6" = "#2E7D32",
  "A6" = "#4E79A7",
  "BL2969" = "#D95F5F"
)


#### Plot PCoA ####

#### Plot PCoA ####

p_pcoa <- ggplot(
  pcoa.df,
  aes(x = PCoA1, y = PCoA2, colour = genome)
) +
  geom_point(
    data = subset(pcoa.df, genome %in% c("A6", "BL2969")),
    aes(size = count),
    alpha = 0.65
  ) +
  geom_point(
    data = subset(pcoa.df, genome == "dm6"),
    size = 2,
    alpha = 0.9
  ) +
  scale_colour_manual(
    name = "Genome",
    values = genome.colours
  ) +
  scale_size_continuous(
    name = "Copy number",
    range = c(1.5, 7)
  ) +
  theme_classic() +
  labs(
    x = xlab,
    y = ylab
  )


p_pcoa

ggplot2::ggsave(
  "outputs/plots/kodama_bar_a6_dm6_PCoA.pdf",
  width = 7,
  height = 5
)
