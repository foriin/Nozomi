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

#### Find Overlaps ####
ove <- findOverlaps(a6_dm6, bl_dm6)

a6blove <- cbind(
  as.data.frame(a6_dm6)[ove@from,],
  as.data.frame(bl_dm6)[ove@to,]
)

subsetByOverlaps(a6_dm6, bl_dm6, invert = T, ignore.strand = T) %>% as.data.frame %>% View

subsetByOverlaps(bl_dm6,a6_dm6,  invert = T, ignore.strand = T) %>% as.data.frame %>% View
length(unique(ove@to))

#### Import original hopper insertions BEDs ####
a6_hop.bed <- import.bed("~/Work/projects/DrosoTE/Hopper_A6/a6_lhop.bed")
a6_hop.bed <- GenomicRanges::reduce(a6_hop.bed)
names(a6_hop.bed) <- paste("A6", "width", width(a6_hop.bed), "hopper", 1:length(a6_hop.bed),
                           sep = "_")
export.bed(a6_hop.bed, "~/Work/projects/DrosoTE/Hopper_A6/a6_lhop.mod.bed")

bl2969_hop.bed <- import.bed("~/Work/projects/DrosoTE/Hopper_A6/bl2969_lhop.bed")
bl2969_hop.bed <- GenomicRanges::reduce(bl2969_hop.bed)
names(bl2969_hop.bed) <- paste("bl2969", "width", width(bl2969_hop.bed), "hopper", 1:length(bl2969_hop.bed),
                           sep = "_")
export.bed(bl2969_hop.bed, "~/Work/projects/DrosoTE/Hopper_A6/bl2969_lhop.mod.bed")

#### Liftover from BL to A6 ####
bl_a6_hop.bed <- import.bed("~/Work/projects/DrosoTE/Hopper_A6/LIFT/bl2969_lhop.a6.sorted.bed")
bl_a6_hop.bed <- sortSeqlevels(bl_a6_hop.bed)
bl_a6_hop.bed <- sort(bl_a6_hop.bed)
as_data_frame(bl_a6_hop.bed) %>% View

subsetByOverlaps(bl_a6_hop.bed, a6_hop.bed, ignore.strand = T) %>% as.data.frame() %>% View()

#### Align hoppers from 2 genomes to each other ####

###### Get sequences for BL2969 hoppers ####
blg <- "~/Work/projects/DrosoTE/Hopper_A6/genomes/GCA_055768665.1_Dmel_BL2969_genomic.fa"
blg.fa <- FaFile(blg)
open(blg.fa)

barhop <- getSeq(blg.fa, bl2969_hop.bed)
names(barhop) <- names(bl2969_hop.bed)
weirdbhops <- barhop[width(barhop) > 2860]
wbsplit1 <- subseq(weirdbhops, 1, c(1432, 1433, 1431, 1432))
names(wbsplit1) <- paste0(
  mapply(
    function(nm, w) sub("width_[0-9]+", paste0("width_", w), nm),
    names(wbsplit1),
    width(wbsplit1),
    USE.NAMES = FALSE
  ),
  "L"
)
wbsplit2 <- subseq(weirdbhops, c(1432, 1433, 1431, 1432) + 1, width(weirdbhops))
names(wbsplit2) <- paste0(
  mapply(
    function(nm, w) sub("width_[0-9]+", paste0("width_", w), nm),
    names(wbsplit2),
    width(wbsplit2),
    USE.NAMES = FALSE
  ),
  "R"
)
barhop <- c(
  barhop[!(names(barhop) %in% names(weirdbhops))],
  wbsplit1,
  wbsplit2
)

barhop_bothtirs <- barhop[startsWith(as.character(barhop), "CACTAT") & endsWith(as.character(barhop), "ATAGTG")]

bar_seq_counts <- sort(table(as.character(barhop_bothtirs)), decreasing = TRUE)

bar.unique <- DNAStringSet(names(bar_seq_counts))
mcols(bar.unique)$count <- as.integer(bar_seq_counts)
names(bar.unique) <- paste("BL2969", 1:length(bar.unique), "w", width(bar.unique),
                           "n", bar_seq_counts, sep = "_")

bar.unique@elementMetadata
bar.unique@elementMetadata$count %>% hist(breaks = 30)

barun.aln <- AlignSeqs(bar.unique)
dm <- DistanceMatrix(barun.aln, type = "dist")
hc <- hclust(dm, method = 'average')
BrowseSeqs(barun.aln[rev(hc$order)])
barun.aln <- barun.aln[rev(hc$order)]

###### Collapse similar sequences ######
collapse_zero_dist <- function(aln, hc) {
  
  nm <- names(aln)
  
  w <- as.integer(sub(".*_w_([0-9]+)_n_[0-9]+$", "\\1", nm))
  n <- as.integer(sub(".*_w_[0-9]+_n_([0-9]+)$", "\\1", nm))
  
  grp <- cutree(hc, h = 0)
  
  keep <- unlist(lapply(split(seq_along(aln), grp), function(ii) {
    ii[which.max(w[ii])]
  }))
  
  n_sum <- sapply(split(seq_along(aln), grp), function(ii) sum(n[ii]))
  
  out <- aln[keep]
  
  names(out) <- mapply(function(old_name, new_n) {
    sub("_n_[0-9]+$", paste0("_n_", new_n), old_name)
  }, names(out), n_sum)
  
  out
}

bar.uq.aln.collapsed <- collapse_zero_dist(barun.aln, hc)

ref_i <- which.max(width(bar.uq.aln.collapsed))
dm2 <- as.matrix(DistanceMatrix(bar.uq.aln.collapsed, type = "dist"))
ord <- order(dm2[ref_i, ])

bar.uq.aln.sorted <- bar.uq.aln.collapsed[ord]

BrowseSeqs(bar.uq.aln.sorted)


barhop <- reverseComplement(c(unique(barhop_bothtirs[width(barhop_bothtirs) < 2500]),
                              barhop_bothtirs[width(barhop_bothtirs) > 2500]))

barhop.aln <- AlignSeqs(reverseComplement(unique(barhop_bothtirs)))

BrowseSeqs(barhop.aln)

###### Get sequences for A6 hoppers ######
a6fasta <- "/Users/artemilin/Work/projects/DrosoTE/popTE/earlgrey/dsrp_genomes/A6-Wild5B_pacbio.fa"

genomio <- FaFile(a6fasta)
open(genomio)
a6hopper.fa <- getSeq(genomio, a6_hop.bed)
names(a6hopper.fa) <- names(a6_hop.bed)

a6hop_bothtirs <- a6hopper.fa[startsWith(as.character(a6hopper.fa), "CACTAT") & endsWith(as.character(a6hopper.fa), "ATAGTG")]
a6hbt.aln <- AlignSeqs(unique(a6hop_bothtirs))
dm <- DistanceMatrix(a6hbt.aln, type = 'dist')
hc <- hclust(dm)
BrowseSeqs(a6hbt.aln[hc$order])

a6_seq_counts <- sort(table(as.character(a6hop_bothtirs)), decreasing = TRUE)

a6.unique <- DNAStringSet(names(a6_seq_counts))
mcols(a6.unique)$count <- as.integer(a6_seq_counts)

names(a6.unique) <- paste("A6", 1:length(a6.unique), "w", width(a6.unique),
                           "n", a6_seq_counts, sep = "_")

a6.unique@elementMetadata
a6.unique@elementMetadata$count %>% hist(breaks = 30)

a6.uq.aln <- AlignSeqs(a6.unique[width(a6.unique) > 1300 &
                                   names(a6.unique) != "A6_43_w_1436_n_1"])
dm <- DistanceMatrix(a6.uq.aln, type = 'dist')
hc <- hclust(dm, method = 'average')
BrowseSeqs(a6.uq.aln[rev(hc$order)])

###### Collapse similar sequences ######


a6.uq.aln.collapsed <- collapse_zero_dist(a6.uq.aln, hc)

ref_i <- which.max(width(a6.uq.aln.collapsed))
dm2 <- as.matrix(DistanceMatrix(a6.uq.aln.collapsed, type = "dist"))
ord <- order(dm2[ref_i, ])

a6.uq.aln.sorted <- a6.uq.aln.collapsed[ord]

save(bar.uq.aln.sorted,
     a6.uq.aln.sorted,
     file = 'outputs/RData/bar_a6_unique_collapsed_aln.RData')


wild5hop <- unique(a6hopper.fa)
names(wild5hop) <- paste("A6", names(wild5hop), sep = "_")

###### Align unique hopper seqs from both genomes and dm6 ######
names(dm6hop_bothtirs)
bwhopdm6 <- c(bar.unique[bar.unique@elementMetadata$count > 2],
              a6.unique[a6.unique@elementMetadata$count > 2],
              dm6hop_bothtirs[width(dm6hop_bothtirs) > 1300])
mcols(bwhopdm6)$genome <- c(
  rep("BL2969", length(bar.unique)),
  rep("A6", length(a6.unique)),
  rep("dm6", length(dm6hop_bothtirs[width(dm6hop_bothtirs) > 1300]))
)
bwhopdm6.aln <- AlignSeqs(bwhopdm6)
dm <- DistanceMatrix(bwhopdm6.aln, type = 'dist')
hc <- hclust(dm, method = 'average')
BrowseSeqs(bwhopdm6.aln[hc$order])

# Plot tree
library(ape)
tr <- bionj(dm)
plot(tr, cex = 0.6, no.margin = TRUE)

#### Karyoplots for hopper insertions ####
new_seqinfo <- Seqinfo(
  seqnames = c("2L", "2R", "3L", "3R", "X", "4"),
  seqlengths = c(23513712, 25286936, 28110227, 32079331, 23542271, 1348131),
  genome = "dm6"
)

## prepare BAR
barlonghops.dm6 <- bl_dm6[bl_dm6$width.orig > 2700]
barlonghops.dm6$sample <- "BL2969"

seqlevels(barlonghops.dm6) <- seqlevels(new_seqinfo)
seqinfo(barlonghops.dm6) <- new_seqinfo
barlonghops.dm6 <- sortSeqlevels(barlonghops.dm6)
barlonghops.dm6 <- sort(barlonghops.dm6)

## prepare A6
a6longhops.dm6 <- a6_dm6[a6_dm6$width.orig > 2700]
a6longhops.dm6$sample <- "A6"

seqlevels(a6longhops.dm6) <- seqlevels(new_seqinfo)
seqinfo(a6longhops.dm6) <- new_seqinfo
a6longhops.dm6 <- sortSeqlevels(a6longhops.dm6)
a6longhops.dm6 <- sort(a6longhops.dm6)

# prepare B3
b3hops.dm6 <- b3_dm6
b3hops.dm6$sample <- "B3"
seqlevels(b3hops.dm6) <- seqlevels(new_seqinfo)
seqinfo(b3hops.dm6) <- new_seqinfo

# prepare reference dm6
hops.dm6$sample <- "dm6"
hops.dm6 <- sortSeqlevels(hops.dm6)
hops.dm6 <- keepStandardChromosomes(hops.dm6, pruning = 'coarse')
seqlevels(hops.dm6) <- seqlevels(new_seqinfo)
seqinfo(hops.dm6) <- new_seqinfo

## combine
longhops.dm6 <- c(barlonghops.dm6, a6longhops.dm6)
# as.data.frame(longhops.dm6) %>% View

###### Add pirna clusters ######
pirna.cl <- import.bed("/Users/artemilin/Work/projects/DrosoTE/Hopper_A6/piRNA/Dmel_proTRAC_piC_collapsed_remapped_frequency.bed")
pirnacl.cons <- pirna.cl[pirna.cl$name == '8']

## check this once if it fails:
## colnames(pirna.tsv)

mcols(pirnacl.cons)$sample <- "piRNA clusters"

pirna.clusters.dm6 <- keepSeqlevels(
  pirnacl.cons,
  value = seqlevels(new_seqinfo),
  pruning.mode = "coarse"
)

seqinfo(pirna.clusters.dm6) <- new_seqinfo
pirna.clusters.dm6 <- sortSeqlevels(pirna.clusters.dm6)
pirna.clusters.dm6 <- sort(pirna.clusters.dm6)

## combine and plot
karyo.gr <- c(
  barlonghops.dm6,
  a6longhops.dm6,
  pirna.clusters.dm6
)

karyo.gr$sample <- factor(
  karyo.gr$sample,
  levels = c("BL2969", "A6", "piRNA clusters")
)

pdf("outputs/plots/longhoppers_a6_bar1_w_pirna.pdf", width = 10, height = 6)
autoplot(
  karyo.gr,
  layout = "karyogram",
  aes(color = sample, fill = sample)
) +
  scale_color_manual(
    values = c(
      "BL2969" = alpha("red", 0.35),
      "A6"     = alpha("blue", 0.35),
      "piRNA clusters" = alpha("green", 0.3)
    )
  ) +
  scale_fill_manual(
    values = c(
      "BL2969" = alpha("red", 0.35),
      "A6"     = alpha("blue", 0.35),
      "piRNA clusters" = alpha("green", 0.3)
    )
  )+
  theme_pack_panels() +
  ggtitle("Long hopper insertions in BL2969 and A6")
dev.off()

pdf("outputs/plots/nozomi_a6_bar1_karyo.pdf", width = 10, height = 6)
autoplot(
  longhops.dm6,
  layout = "karyogram",
  aes(color = sample, fill = sample)
) +
  scale_color_manual(
    values = c(
      "BL2969" = alpha("red", 0.35),
      "A6"     = alpha("blue", 0.35)
    )
  ) +
  scale_fill_manual(
    values = c(
      "BL2969" = alpha("red", 0.35),
      "A6"     = alpha("blue", 0.35),
      "piRNA clusters" = alpha("green", 0.3)
    )
  )+
  theme_pack_panels() +
  ggtitle("Nozomi insertions in BL2969 and A6")
dev.off()

hops.dm6$sample <- "dm6"
seqlevelsStyle(pirna.clusters.dm6) <- "UCSC"
seqlevels(pirna.clusters.dm6) <- paste0('chr', seqlevels(pirna.clusters.dm6))
seqlevels(b3hops.dm6) <- paste0("chr", seqlevels(b3hops.dm6))

seqinfo

pdf("outputs/plots/hoppers_b3_dm6_pi_karyo.pdf", width = 10, height = 6)
autoplot(
  c(hops.dm6,
    b3hops.dm6,
    pirna.clusters.dm6),
  layout = "karyogram",
  aes(color = sample, fill = sample)
) +
  scale_color_manual(
    values = c(
      "dm6" = alpha("gray", 0.5),
      "B3" = alpha("violet", 0.35),
      "piRNA clusters"     = alpha("green", 0.35)
    )
  ) +
  scale_fill_manual(
    values = c(
      "dm6" = alpha("gray", 0.5),
      "B3" = alpha("violet", 0.35),
      "piRNA clusters" = alpha("green", 0.3)
    )
  )+
  theme_pack_panels() +
  ggtitle("Hopper insertions in B3 and dm6")
dev.off()

