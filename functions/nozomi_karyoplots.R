library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
library(ggbio)
library(ape)
library(Rsamtools)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(DECIPHER)

#### Import liftovers of insertions from B3, A6 and BL2969 to dm6 ####
a6_dm6 <- import.bed("~/Work/projects/DrosoTE/Hopper_A6/LIFT/a6_dm6.sorted.bed")
a6_dm6 <- sortSeqlevels(a6_dm6)
a6_dm6 <- GenomicRanges::sort(a6_dm6)
mcols(a6_dm6)$width.orig <- as.integer(
  sub(".*_width_([0-9]+)_hopper_.*", "\\1", mcols(a6_dm6)$name)
)

bl_dm6 <- import.bed("~/Work/projects/DrosoTE/Hopper_A6/LIFT/bl2969_dm6.sorted.bed")
bl_dm6 <- sortSeqlevels(bl_dm6)
bl_dm6 <- sort(bl_dm6)
mcols(bl_dm6)$width.orig <- as.integer(
  sub(".*_width_([0-9]+)_hopper_.*", "\\1", mcols(bl_dm6)$name)
)

b3_dm6 <- import.bed("/Users/artemilin/Work/projects/DrosoTE/Hopper_A6/LIFT/B3_dm6.sorted.bed")



#### Karyoplots for hopper insertions ####

###### Prepare common dm6 seqinfo ######

dm6_chr <- paste0("chr", c("2L", "2R", "3L", "3R", "X", "4"))

dm6_seqinfo <- seqinfo(
  BSgenome.Dmelanogaster.UCSC.dm6
)[dm6_chr]


###### Standardize ranges to dm6 ######

prepare_dm6_ranges <- function(gr, sample) {
  
  # Keep only major chromosomes
  chr_bare <- sub("^chr", "", seqlevels(gr))
  keep <- chr_bare %in% sub("^chr", "", dm6_chr)
  
  gr <- keepSeqlevels(
    gr,
    seqlevels(gr)[keep],
    pruning.mode = "coarse"
  )
  
  # Convert 2L -> chr2L, etc.
  old_levels <- seqlevels(gr)
  new_levels <- paste0("chr", sub("^chr", "", old_levels))
  
  gr <- renameSeqlevels(
    gr,
    setNames(new_levels, old_levels)
  )
  
  # Supply missing chromosome lengths
  seqinfo(gr) <- dm6_seqinfo[seqlevels(gr)]
  
  # Remove incompatible metadata before combining
  mcols(gr) <- S4Vectors::DataFrame(
    sample = rep(sample, length(gr))
  )
  
  gr <- sortSeqlevels(gr)
  gr <- GenomicRanges::sort(gr)
  
  gr
}


#### Long Hopper insertions in BL2969 and A6 ####

###### Prepare BL2969 ######

barlonghops.dm6 <- bl_dm6[
  !is.na(bl_dm6$width.orig) &
    bl_dm6$width.orig > 2700
]

barlonghops.dm6 <- barlonghops.dm6[-c(5, 9, 10, 11)]

barlonghops.dm6 <- prepare_dm6_ranges(
  barlonghops.dm6,
  sample = "BL2969"
)


###### Prepare A6 ######

a6longhops.dm6 <- a6_dm6[
  !is.na(a6_dm6$width.orig) &
    a6_dm6$width.orig > 2700
]

a6longhops.dm6 <- prepare_dm6_ranges(
  a6longhops.dm6,
  sample = "A6"
)


###### Combine BL2969 and A6 ######

longhops.dm6 <- c(
  barlonghops.dm6,
  a6longhops.dm6
)

longhops.dm6$sample <- factor(
  longhops.dm6$sample,
  levels = c("BL2969", "A6")
)


###### Plot BL2969 and A6 ######

pdf(
  "outputs/plots/nozomi_a6_bl2969_karyo.pdf",
  width = 6,
  height = 2.5
)

autoplot(
  longhops.dm6,
  layout = "karyogram",
  aes(color = sample, fill = sample)
) +
  scale_color_manual(
    values = c(
      "BL2969" = alpha("red", 0.35),
      "A6" = alpha("blue", 0.35)
    )
  ) +
  scale_fill_manual(
    values = c(
      "BL2969" = alpha("red", 0.35),
      "A6" = alpha("blue", 0.35)
    )
  ) +
  theme_pack_panels() +
  ggtitle("Nozomi insertions in BL2969 and A6")

dev.off()


#### Long Hopper insertions with piRNA clusters ####

###### Import piRNA clusters ######

# Comes from here: https://github.com/kerogens101/Dmel_piCs/tree/main/piC_annotations

pirna.cl <- import.bed(
  "/Users/artemilin/Work/projects/DrosoTE/Hopper_A6/piRNA/Dmel_proTRAC_piC_collapsed_remapped_frequency.bed"
)

pirnacl.cons <- pirna.cl[
  pirna.cl$name %in% c(as.character(7:8))
]


###### Prepare piRNA clusters ######

pirna.clusters.dm6 <- prepare_dm6_ranges(
  pirnacl.cons,
  sample = "piRNA clusters"
)


###### Combine insertions and piRNA clusters ######

karyo.gr <- c(
  barlonghops.dm6,
  a6longhops.dm6,
  pirna.clusters.dm6
)

karyo.gr$sample <- factor(
  karyo.gr$sample,
  levels = c(
    "BL2969",
    "A6",
    "piRNA clusters"
  )
)


###### Plot insertions and piRNA clusters ######

pdf(
  "outputs/plots/nozomi_a6_bl2969_w_pirna.pdf",
  width = 6,
  height = 3
)

autoplot(
  karyo.gr,
  layout = "karyogram",
  aes(color = sample, fill = sample)
) +
  scale_color_manual(
    values = c(
      "BL2969" = alpha("red", 0.35),
      "A6" = alpha("blue", 0.35),
      "piRNA clusters" = alpha("green", 0.30)
    )
  ) +
  scale_fill_manual(
    values = c(
      "BL2969" = alpha("red", 0.35),
      "A6" = alpha("blue", 0.35),
      "piRNA clusters" = alpha("green", 0.30)
    )
  ) +
  theme_pack_panels() +
  ggtitle("Nozomi insertions and piRNA clusters")

dev.off()


#### B3, dm6 and piRNA clusters ####

###### Prepare B3 insertions ######

b3hops.dm6 <- prepare_dm6_ranges(
  b3_dm6,
  sample = "B3"
)


###### Prepare reference dm6 insertions ######

hops.dm6.plot <- prepare_dm6_ranges(
  hops.dm6,
  sample = "dm6"
)


###### Combine B3, dm6 and piRNA clusters ######

b3.karyo.gr <- c(
  hops.dm6.plot,
  b3hops.dm6,
  pirna.clusters.dm6
)

b3.karyo.gr$sample <- factor(
  b3.karyo.gr$sample,
  levels = c(
    "dm6",
    "B3",
    "piRNA clusters"
  )
)


###### Plot B3, dm6 and piRNA clusters ######

pdf(
  "outputs/plots/hoppers_b3_dm6_pi_karyo.pdf",
  width = 10,
  height = 6
)

autoplot(
  b3.karyo.gr,
  layout = "karyogram",
  aes(color = sample, fill = sample)
) +
  scale_color_manual(
    values = c(
      "dm6" = alpha("gray", 0.50),
      "B3" = alpha("violet", 0.35),
      "piRNA clusters" = alpha("green", 0.35)
    )
  ) +
  scale_fill_manual(
    values = c(
      "dm6" = alpha("gray", 0.50),
      "B3" = alpha("violet", 0.35),
      "piRNA clusters" = alpha("green", 0.30)
    )
  ) +
  theme_pack_panels() +
  ggtitle("Nozomi insertions in B3 and dm6")

dev.off()