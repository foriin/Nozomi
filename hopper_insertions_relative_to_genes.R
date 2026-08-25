library(GenomicRanges)
library(GenomeInfoDb)
library(ChIPseeker)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(org.Dm.eg.db)
library(ggplot2)

txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene

## choose hopper ranges
# barlonghops.dm6 <- bl_dm6[bl_dm6$name %in% names(barhop_bothtirs)]
barhops.dm6 <- bl_dm6
# a6longhops.dm6  <- a6_dm6[a6_dm6$name %in% names(a6hop_bothtirs)]
a6hops.dm6  <- a6_dm6
b3hops.dm6 <- b3_dm6


#### Get dm6 hopper insertions, pick those that have TIRs ####
hops.dm6 <- import.bed("~/Work/projects/DrosoTE/Hopper_A6/dm6_lhop.bed")
hops.dm6 <- GenomicRanges::reduce(hops.dm6)
names(hops.dm6) <- paste("dm6", "w", width(hops.dm6), "id", 1:length(hops.dm6),
                         sep = "_")

dm6genome <- FaFile("~/Work/db/dmel-all-chromosome-r6.67.fasta")
open(dm6genome)
hops.dm6.fa <- getSeq(dm6genome, hops.dm6)
names(hops.dm6.fa) <- names(hops.dm6)
close(dm6genome)

hist(width(hops.dm6.fa))
unique(hops.dm6.fa)
dm6hop_bothtirs <- hops.dm6.fa[startsWith(as.character(hops.dm6.fa), "CACTAT") & endsWith(as.character(hops.dm6.fa), "ATAGTG")]


## make seqlevels compatible with dm6 TxDb: chr2L, chr2R, chr3L, chr3R, chrX, chr4
seqlevelsStyle(barhops.dm6) <- "UCSC"
seqlevelsStyle(a6hops.dm6)  <- "UCSC"
seqlevelsStyle(b3hops.dm6)  <- "UCSC"
seqlevelsStyle(hops.dm6)  <- "UCSC"

barhops.dm6 <- keepStandardChromosomes(barhops.dm6, pruning.mode = "coarse")
a6hops.dm6  <- keepStandardChromosomes(a6hops.dm6, pruning.mode = "coarse")
b3hops.dm6  <- keepStandardChromosomes(b3hops.dm6, pruning.mode = "coarse")
hops.dm6  <- keepStandardChromosomes(hops.dm6, pruning.mode = "coarse")

## for TE insertion sites, annotate midpoint rather than full TE span
bar.ins <- resize(barhops.dm6, width = 1, fix = "center")
a6.ins  <- resize(a6hops.dm6,  width = 1, fix = "center")
dm6.ins <- resize(hops.dm6, width = 1, fix = "center")

## annotate relative to genes
bar.anno <- annotatePeak(
  bar.ins,
  TxDb = txdb,
  tssRegion = c(-1000, 1000),
  annoDb = "org.Dm.eg.db"
)

a6.anno <- annotatePeak(
  a6.ins,
  TxDb = txdb,
  tssRegion = c(-1000, 1000),
  annoDb = "org.Dm.eg.db"
)

dm6.anno <- annotatePeak(
  dm6.ins,
  TxDb = txdb,
  tssRegion = c(-1000, 1000),
  annoDb = "org.Dm.eg.db"
)

## basic annotation pie/bar plots
plotAnnoPie(bar.anno) + ggtitle("BL2969 hopper insertion annotation")
plotAnnoPie(a6.anno)  + ggtitle("A6 hopper insertion annotation")
plotAnnoPie(dm6.anno)  + ggtitle("dm6 hopper insertion annotation")

plotAnnoBar(bar.anno) + ggtitle("BL2969 hopper insertion annotation")
plotAnnoBar(a6.anno)  + ggtitle("A6 hopper insertion annotation")

## compare both samples in one plot
hopper.list <- list(
  BL2969 = bar.ins,
  A6 = a6.ins,
  dm6 = dm6.ins
)

hopper.anno.list <- lapply(
  hopper.list,
  annotatePeak,
  TxDb = txdb,
  tssRegion = c(-1000, 1000),
  annoDb = "org.Dm.eg.db"
)


p_bla6 <- plotAnnoBar(hopper.anno.list) +coord_cartesian()


feature_order <- c(
  "Promoter", "1st Intron", "Other Intron", "Distal Intergenic",
  "Downstream (<=300)", "5' UTR", "3' UTR", "Other Exon", "1st Exon"
)

p_bla6$data$Feature <- factor(p_bla6$data$Feature, levels = rev(feature_order))

p_bla6 <- p_bla6 +
  scale_fill_manual(
    values = c(
      "Promoter" = "#11CF11",
      "1st Intron" = "#3377Fd",
      "Other Intron" = "lightblue",
      "Distal Intergenic" = "lightblue4",
      "Downstream (<=300)" = "yellow3",
      "5' UTR" = "purple",
      "3' UTR" = "violet",
      "Other Exon" = "orange",
      "1st Exon" = "red3"
    ),
    limits = feature_order,
    breaks = feature_order,
    drop = FALSE
  )

ggplot2::ggsave(
  "outputs/plots/bl_a6_dm6_kodama_feature_distr.pdf",
  p_bla6,
  width = 4,
  height = 4
)

## distance to nearest TSS
plotDistToTSS(bar.anno, title = "BL2969 hopper insertions relative to TSS")
plotDistToTSS(a6.anno,  title = "A6 hopper insertions relative to TSS")

plotDistToTSS(
  hopper.anno.list,
  title = "Hopper insertions relative to TSS"
)

## extract annotation table
bar.anno.df <- as.data.frame(bar.anno)
a6.anno.df  <- as.data.frame(a6.anno)

write.table(
  bar.anno.df,
  "BL2969_hopper_ChIPseeker_annotation.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  a6.anno.df,
  "A6_hopper_ChIPseeker_annotation.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)