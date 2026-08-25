library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
library(ggbio)

#### Load data ####

load(file = "../TE_TF/RData/te.hierarchy.Artem.2026.RData", verb = T)
te.hier$ID <- toupper(te.hier$name)

gffs <- dir("/Users/artemilin/Work/projects/DrosoTE/popTE/earlgrey/dm6_liftover_earlgrey", pattern ="*filtered.gff",
            full.names = T)
gna <- sub("_pacbio.dm6.filteredRepeats.filtered.gff", "", basename(gffs))

names(gffs) <- gna

gff_orig <- dir("/Users/artemilin/Work/projects/DrosoTE/popTE/earlgrey/dm6_liftover_earlgrey/eg_final/", pattern =".gff",
            full.names = T)
gna2 <- sub("_pacbio.filteredRepeats.gff", "", basename(gff_orig))

names(gff_orig) <- gna2

genomefolder <- dir("/Users/artemilin/Work/projects/DrosoTE/popTE/earlgrey/dsrp_genomes/", pattern =".fa",
                    full.names = T)
names(genomefolder) <- gna2 # they match, but better to check

tefa <- readDNAStringSet("~/Work/db/Dmel_TE_cons_Artem.srt.fa")
#### Prepare reference, non-reference and isolated insertions ####

refreps <- lapply(gna, function(nnnnn){
  reps <- import.gff(gffs[[nnnnn]])
  reps <- reps[width(reps) > 100]
  reps$strain <- nnnnn
  reps
})

refreps <- unlist(GRangesList(refreps))
table(seqnames(refreps))
table(refreps$NAME) %>% sort() %>% tail(n = 20)
table(refreps$strain)
table(refreps$type) %>% sort

nonrepreps <- lapply(gna, function(nnnnn){
  reps <- import.gff(gffs[[nnnnn]])
  reps <- reps[width(reps) == 1]
  reps$strain <- nnnnn
  reps
})

nonrepreps <- unlist(GRangesList(nonrepreps))
dtn <- distanceToNearest(nonrepreps)

isolatedreps <- nonrepreps[dtn@from][dtn@elementMetadata$distance > 5]

table(isolatedreps$NAME) %>% sort
table(seqnames(isolatedreps))[1:6]/seqlengths(BSgenome.Dmelanogaster.UCSC.dm6)[1:6]
table(isolatedreps$type) %>% sort()
mapply(function(x,y) x/y,
       table(seqnames(isolatedreps))[1:6],
      seqlengths(BSgenome.Dmelanogaster.UCSC.dm6)[1:6])
table(isolatedreps$strain)
names(isolatedreps) <- gna


newins <- lapply(gna, function(x) {
  grs <- isolatedreps[isolatedreps$strain == x]
  print(grs)
  table(grs$NAME) %>% enframe
})
names(newins) <- gna

df <- bind_rows(newins, .id = 'strain')
dfw <- pivot_wider(df, names_from = strain, values_from = value)
dfwf <- dfw[rowSums(dfw[,-1], na.rm = T) > 50,]


df_long <- dfwf %>%
  pivot_longer(
    cols = -name,
    names_to = "genome",
    values_to = "value"
  ) %>%
  mutate(
    value = as.numeric(value),
    genome = factor(genome, levels = setdiff(names(dfwf), "name"))
  )

df_long$value[is.na(df_long$value)] <- 0

library(ggrepel)

ggplot(df_long, aes(x = genome, y = value, col = name,group = name)) +
  geom_line(alpha = 0.4) +
  geom_point(size = 1.8) +
  geom_text_repel(aes(label = name), size = 2.5) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(df_long %>% filter(name %in% c("ROO", "HOPPER", "DOC", "COPIA", "POGO", "F-ELEMENT", "HOBO")), aes(x = genome, y = value, col = name,group = name)) +
  geom_line(alpha = 0.4) +
  geom_point(size = 1.8) +
  geom_text_repel(aes(label = name), size = 2.5) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


#### Karyo plots for new insertions ####


b2_ca1_newdoc <- isolatedreps[isolatedreps$strain == "B2-CA1" & isolatedreps$NAME == "DOC"]

new_seqinfo <- Seqinfo(
  seqnames = c("X", "2L", "2R", "3L", "3R", "4"),
  seqlengths = c(23542271, 23513712, 25286936, 28110227, 32079331, 1348131),
  genome = "dm6"
)

seqlevels(b2_ca1_newdoc) <- seqlevels(new_seqinfo)
seqinfo(b2_ca1_newdoc) <- new_seqinfo

pdf("plots/b2_ca1_iso_Doc_karyo.pdf", width = 6, height = 4)
autoplot(b2_ca1_newdoc, layout = "karyogram")+
  theme_pack_panels()+
  ggtitle("Non-reference isolated Doc insertions in B2-CA1")
dev.off()

a6wild5b_newhopper <- isolatedreps[isolatedreps$strain == "A6-Wild5B" & isolatedreps$NAME == "HOPPER"]
seqlevels(a6wild5b_newhopper) <- seqlevels(new_seqinfo)
seqinfo(a6wild5b_newhopper) <- new_seqinfo

pdf("plots/a6_iso_hopper.pdf", width = 6, height = 4)
autoplot(a6wild5b_newhopper, layout = "karyogram")+
  ggtitle("A6-Wild5b new hopper insertions")+
  theme_pack_panels()+
  ggtitle("Non-reference isolated hopper insertions in A6-Wild5B")
dev.off()

table(seqnames(a6wild5b_newhopper))
table(seqnames(b2_ca1_newdoc))

a5vag1_newcopia <- isolatedreps[isolatedreps$strain == "A5-VAG1" & isolatedreps$NAME == "COPIA"]
seqlevels(a5vag1_newcopia) <- seqlevels(new_seqinfo)
seqinfo(a5vag1_newcopia) <- new_seqinfo

pdf("plots/a5_iso_copia.pdf", width = 6, height = 4)
autoplot(a5vag1_newcopia, layout = "karyogram")+
  ggtitle("A5-VAG1 new copia insertions")+
  theme_pack_panels()+
  ggtitle("Non-reference isolated copia insertions in A6-Wild5B")
dev.off()
table(seqnames(a5vag1_newcopia))

a1canton_newf<- isolatedreps[isolatedreps$strain == "A1-CantonS" & isolatedreps$NAME == "F-ELEMENT"]
seqlevels(a1canton_newf) <- seqlevels(new_seqinfo)
seqinfo(a1canton_newf) <- new_seqinfo

autoplot(a1canton_newf, layout = "karyogram")+
  ggtitle("A1-CantonS new F-element insertions")
table(seqnames(a1canton_newf))

#### Isolate long insertions, look at their distributions ####
teorig <- lapply(gna2, function(nnnnn){
  reps <- import.gff(gff_orig[[nnnnn]])
  reps <- reps[width(reps) > 500]
  reps$strain <- nnnnn
  reps$strainid <- paste(reps$strain,
                         reps$ID,
                         sep = "_")
  reps
})

teorig <- unlist(GRangesList(teorig))

names(isolatedreps) <- isolatedreps$strain
isolatedreps$strainid <- paste(isolatedreps$strain,
                               isolatedreps$ID,
                               sep = "_")


teorig.te <- teorig[teorig$NAME %in% te.hier$ID]
teorig.te$conslen <- te.hier$length[match(teorig.te$NAME, te.hier$ID)]
# Long insertions
teorig.te.truec <- teorig.te[abs(width(teorig.te) - teorig.te$conslen) < 0.1*teorig.te$conslen]

isolated.long <- isolatedreps[isolatedreps$strainid %in% teorig.te.truec$strainid]

###### Do family distribution plot just for long insertions ####

newins2 <- lapply(gna2, function(x) {
  grs <- isolated.long[isolated.long$strain == x]
  print(grs)
  table(grs$NAME) %>% enframe
})
names(newins2) <- gna2

df2 <- bind_rows(newins2, .id = 'strain')
dfw2 <- pivot_wider(df2, names_from = strain, values_from = value)
dfwf2 <- dfw2[rowSums(dfw2[,-1], na.rm = T) > 5,]


df_long2 <- dfwf2 %>%
  pivot_longer(
    cols = -name,
    names_to = "genome",
    values_to = "value"
  ) %>%
  mutate(
    value = as.numeric(value),
    genome = factor(genome, levels = setdiff(names(dfwf2), "name"))
  )

df_long2$value[is.na(df_long2$value)] <- 0

ggplot(df_long2, aes(x = genome, y = value, col = name,group = name)) +
  geom_line(alpha = 0.4) +
  geom_point(size = 1.8) +
  geom_text_repel(aes(label = name), size = 2.5) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(df_long2 %>% filter(name %in% c("ROO", "HOPPER", "DOC", "COPIA", "POGO", "F-ELEMENT", "HOBO")), aes(x = genome, y = value, col = name,group = name)) +
  geom_line(alpha = 0.4) +
  geom_point(size = 1.8) +
  geom_text_repel(aes(label = name), size = 2.5) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(
  df_long2 %>% filter(name %in% c("ROO", "HOPPER", "DOC", "COPIA", "POGO", "F-ELEMENT", "HOBO")),
  aes(x = genome, y = value, col = name, group = name)
) +
  geom_line(aes(linewidth = name == "ROO")) +
  geom_point(size = 1.8) +
  # geom_text_repel(aes(label = name), size = 2.5) +
  scale_linewidth_manual(values = c(`TRUE` = 1.3, `FALSE` = 0.5), guide = "none") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

###### MSA + Distance tree for hyperactive insertions ######

teorig.te.truec.f <- teorig.te.truec[teorig.te.truec$strainid %in% isolated.long$strainid]
toi <- c("ROO", "HOPPER", "DOC",
         "COPIA", "POGO",
         "F-ELEMENT", "HOBO")
te.fa <- lapply(gna2, function(strain){
  print(strain)
  te.gr <- teorig.te.truec.f[teorig.te.truec.f$strain == strain &
                           teorig.te.truec.f$NAME %in% toi]
  genomio <- FaFile(genomefolder[[strain]])
  open(genomio)
  indexFa(genomio)
  te.fa <- getSeq(genomio, te.gr)
  names(te.fa) <- te.gr$strainid
  te.fa
})

te.fa.c <- do.call(c, te.fa)

te.seq.by.id <- lapply(toi, function(te){
  te.fa.c[grepl(te, names(te.fa.c))]
})
names(te.seq.by.id) <- toi
te.fa.c[grepl('COPIA', names(te.fa.c))]

library(DECIPHER)
library(Biostrings)
library(dplyr)
library(ggplot2)
library(ape)

te.pca <- lapply(te.seq.by.id, function(seqset){
  aln <- AlignSeqs(seqset)
  dm  <- DistanceMatrix(aln, type = "dist")
  
  meta <- data.frame(
    seqname = names(seqset),
    strain  = sub("_.*$", "", names(seqset)),
    te      = sub("^[^_]+_([^_]+)_.*$", "\\1", names(seqset))
  )
  
  pcoa <- cmdscale(as.dist(dm), k = 2, eig = TRUE)
  
  plot_df <- meta %>%
    mutate(
      PC1 = pcoa$points[, 1],
      PC2 = pcoa$points[, 2]
    )
})

ggplot(te.pca$`F-ELEMENT`, aes(PC1, PC2, color = strain)) +
  geom_point(size = 2) +
  theme_bw()

ggplot(te.pca$`F-ELEMENT`, aes(PC1, PC2, col = strain)) +
  geom_density_2d() +
  theme_bw()

ggplot(te.pca$`F-ELEMENT` %>% filter(strain == 'A1-CantonS'), aes(PC1, PC2)) +
  geom_density_2d_filled() +
  theme_bw()
ggplot(te.pca$`F-ELEMENT` %>% filter(strain == 'B3-QI2'), aes(PC1, PC2)) +
  geom_density_2d_filled() +
  theme_bw()

ggplot(te.pca$COPIA, aes(PC1, PC2, color = strain)) +
  geom_point(size = 2) +
  facet_wrap(~strain)+
  theme_bw()

ggplot(te.pca$COPIA, aes(PC1, PC2)) +
  geom_density_2d_filled(
    contour_var = "ndensity"
  ) +
  facet_wrap(~strain)+
  theme_bw()

ggplot(te.pca$DOC, aes(PC1, PC2)) +
  geom_density_2d_filled(
    contour_var = "ndensity"
  ) +
  facet_wrap(~strain)+
  theme_bw()

ggplot(te.pca$COPIA, aes(PC1, PC2)) +
  geom_density_2d_filled() +
  theme_bw()

#### Find new elements in Doc ####
ggplot(te.pca$DOC, aes(PC1, PC2)) +
  geom_density_2d_filled(
    contour_var = "ndensity"
  ) +
  facet_wrap(~strain)+
  theme_bw()

doc_cands <- te.pca$DOC %>% 
  # filter(strain %in% c("A1-CantonS", "B2-CA1"),
  #                                  PC1 < 0) %>% 
  pull(seqname)

docsfa <- te.fa.c[names(te.fa.c) %in% doc_cands]
docref <- tefa[names(tefa) == "Doc#LINE/I-Jockey"]
writeXStringSet(docref, "fasta/Doc_cons.fa")
docsfa <- c(docsfa, docref)

alndocs <- AlignSeqs(docsfa)
BrowseSeqs(alndocs, highlight = 218)

b2_ca1_doc_81895 <- docsfa[names(docsfa) == "B2-CA1_DOC_81895"]
writeXStringSet(b2_ca1_doc_81895, "fasta/B2-CA1_DOC_81895.fa")

#### Find new elements in COPIA ####
ggplot(te.pca$COPIA, aes(PC1, PC2)) +
  geom_density_2d_filled(
    contour_var = "ndensity"
  ) +
  facet_wrap(~strain)+
  theme_bw()

COPIA_cands <- te.pca$COPIA %>%
  # filter(strain %in% c("A5-VAG1", "A6-Wild5B"),
  #                                  PC1 < 1000) %>% 
  pull(seqname)

COPIAsfa <- te.fa.c[names(te.fa.c) %in% COPIA_cands]
COPIAref <- tefa[names(tefa) == "copia#LTR/Copia"]
COPIAsfa <- c(COPIAsfa, COPIAref)
a5_vag1_copia_95440 <- te.fa.c[names(te.fa.c) == "A5-VAG1_COPIA_95440"]
writeXStringSet(a5_vag1_copia_95440, "fasta/A5-VAG1_COPIA_95440.fa")

alnCOPIAs <- AlignSeqs(COPIAsfa)
BrowseSeqs(alnCOPIAs, highlight = 237)

tr <- bionj(dm)

plot(tr, type = "fan", cex = 0.3, no.margin = TRUE)



