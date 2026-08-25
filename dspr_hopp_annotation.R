library(GenomicRanges)
library(GenomeInfoDb)
library(ChIPseeker)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(org.Dm.eg.db)
library(ggplot2)

txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene


#### Load short hopper groups data ####
load("outputs/RData/short_hop_groups_GRanges.RData", verb = T)
# shops.groups.w.stats.gr - GRanges with group info
# k80_by_group - group info

## make seqlevels compatible with dm6 TxDb: chr2L, chr2R, chr3L, chr3R, chrX, chr4
seqlevelsStyle(shops.groups.w.stats.gr) <- "UCSC"

shops.groups.w.stats.gr  <- keepStandardChromosomes(shops.groups.w.stats.gr, pruning.mode = "coarse")

## for TE insertion sites, annotate midpoint rather than full TE span
shops.groups.inserts <- resize(shops.groups.w.stats.gr, width = 1, fix = "center")

## Separate into A6 insertions and others
sgi.a6 <- shops.groups.inserts[grepl("A6", shops.groups.inserts$copies)]
sgi.b3 <- shops.groups.inserts[grepl("B3", shops.groups.inserts$copies)]
sgi.a2 <- shops.groups.inserts[grepl("A2", shops.groups.inserts$copies)]
sgi.noa6 <- shops.groups.inserts[!(grepl("A6", shops.groups.inserts$copies) & shops.groups.inserts$freq_bin == '1 singleton')]


## annotate relative to genes
shop.group.A6.anno <- annotatePeak(
  sgi.a6,
  TxDb = txdb,
  tssRegion = c(-1000, 1000),
  annoDb = "org.Dm.eg.db"
)

shop.group.A2.anno <- annotatePeak(
  sgi.a2,
  TxDb = txdb,
  tssRegion = c(-1000, 1000),
  annoDb = "org.Dm.eg.db"
)

shop.group.B3.anno <- annotatePeak(
  sgi.b3,
  TxDb = txdb,
  tssRegion = c(-1000, 1000),
  annoDb = "org.Dm.eg.db"
)

shop.group.anno.noA6 <- annotatePeak(
  sgi.noa6,
  TxDb = txdb,
  tssRegion = c(-1000, 1000),
  annoDb = "org.Dm.eg.db"
)



## basic annotation pie/bar plots
plotAnnoPie(shop.group.A6.anno) + ggtitle("A6 short hopper insertion annotation")
plotAnnoPie(shop.group.A2.anno) + ggtitle("A2 short hopper insertion annotation")
plotAnnoPie(shop.group.B3.anno) + ggtitle("B3 short hopper insertion annotation")
plotAnnoPie(shop.group.anno.noA6) + ggtitle("All short hopper insertion annotation")


plotAnnoBar(shop.group.A6.anno) + ggtitle("All short hopper insertion locations annotation")
plotAnnoBar(shop.group.anno.noA6) + ggtitle("Short hopper insertion locations (except for A6) annotation")

hops1.4.dm6.ins <- resize(hops.dm6[width(hops.dm6) > 1350], width = 1, fix = "center")

## compare both samples in one plot
shopgr.list.dspr <- list(
  All_ins_noA6 = sgi.noa6,
  A1 = shops.groups.inserts[grepl("A1", shops.groups.inserts$copies)],
  A2 = shops.groups.inserts[grepl("A2", shops.groups.inserts$copies)],
  A3 = shops.groups.inserts[grepl("A3", shops.groups.inserts$copies)],
  A4 = shops.groups.inserts[grepl("A4", shops.groups.inserts$copies)],
  A5 = shops.groups.inserts[grepl("A5", shops.groups.inserts$copies)],
  A6 = sgi.a6,
  A7 = shops.groups.inserts[grepl("A7", shops.groups.inserts$copies)],
  AB8 = shops.groups.inserts[grepl("AB8", shops.groups.inserts$copies)],
  B1 = shops.groups.inserts[grepl("B1", shops.groups.inserts$copies)],
  B2 = shops.groups.inserts[grepl("B2", shops.groups.inserts$copies)],
  B3 = shops.groups.inserts[grepl("B3", shops.groups.inserts$copies)],
  B4 = shops.groups.inserts[grepl("B4", shops.groups.inserts$copies)],
  B6 = shops.groups.inserts[grepl("B6", shops.groups.inserts$copies)],
  ORE = shops.groups.inserts[grepl("ORE", shops.groups.inserts$copies)],
  dm6 = hops1.4.dm6.ins
)

shopgr.list <- list(
  All_ins_noA6 = sgi.noa6,
  A6 = sgi.a6,
  dm6 = hops1.4.dm6.ins
)

shopgr.anno.list <- lapply(
  shopgr.list,
  annotatePeak,
  TxDb = txdb,
  tssRegion = c(-1000, 1000),
  annoDb = "org.Dm.eg.db"
)

plotAnnoBar(shopgr.anno.list)

shopgr.anno.list.dspr <- lapply(
  shopgr.list.dspr,
  annotatePeak,
  TxDb = txdb,
  tssRegion = c(-1000, 1000),
  annoDb = "org.Dm.eg.db"
)

p_dspr <- plotAnnoBar(
  shopgr.anno.list.dspr[c(14, 8, 5, 7, 3, 11, 12, 2, 13, 4, 10, 9, 6, 15, 16)]
)

feature_order <- c(
  "Promoter", "1st Intron", "Other Intron", "Distal Intergenic",
  "Downstream (<=300)", "5' UTR", "3' UTR", "Other Exon", "1st Exon"
)

p_dspr$data$Feature <- factor(p_dspr$data$Feature, levels = rev(feature_order))

p_dspr <- p_dspr +
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
  "outputs/plots/dspr_kodama_feature_distr2.pdf",
  p_dspr+ggtitle("Feature distribution relative to Kodama insertions"),
  width = 9,
  height = 6.5
)

## distance to nearest TSS
plotDistToTSS(shop.group.A6.anno, title = "short hopper insertions (no A6) relative to TSS")
plotDistToTSS(shop.group.anno.noA6,  title = "A6 hopper insertions relative to TSS")

plotDistToTSS(
  shopgr.anno.list,
  title = "Hopper insertions relative to TSS"
)

plotDistToTSS(
  shopgr.anno.list.dspr,
  title = "Hopper 1.4 kb insertions relative to TSS"
)

## extract annotation table
shop.group.A6.anno.df <- as.data.frame(shop.group.A6.anno)
shop.group.anno.noA6.df  <- as.data.frame(shop.group.anno.noA6)
table(shop.group.anno.noA6.df$freq_bin)

ggplot(shop.group.anno.noA6.df, aes (x = freq_bin, y = abs(distanceToTSS)))+
  geom_boxplot()+
  geom_point(alpha = 0.3)+
  theme_bw()


ggplot(shop.group.anno.noA6.df%>% filter(median_k80_to_longhop < 0.11), 
       aes (x = median_k80_to_longhop, y = abs(distanceToTSS)))+
  geom_point(size = 0.4)+
  geom_smooth(method = 'lm', se = F)+
  theme_bw()

ggplot(shop.group.anno.noA6.df, aes (x = median_k80, y = abs(distanceToTSS)))+
  geom_point(size = 0.4)+
  theme_bw()

# write.table(
#   bar.anno.df,
#   "BL2969_hopper_ChIPseeker_annotation.tsv",
#   sep = "\t",
#   quote = FALSE,
#   row.names = FALSE
# )
# 
# write.table(
#   a6.anno.df,
#   "A6_hopper_ChIPseeker_annotation.tsv",
#   sep = "\t",
#   quote = FALSE,
#   row.names = FALSE
# )