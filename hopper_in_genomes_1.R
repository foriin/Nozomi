library(DECIPHER)
library(Biostrings)
library(dplyr)
library(ggplot2)
library(ape)
library(Rsamtools)
library(pheatmap)


a6te <- import.gff("/Users/artemilin/Work/projects/DrosoTE/popTE/earlgrey/dm6_liftover_earlgrey/eg_final/A6-Wild5B_pacbio.filteredRepeats.gff")

a6hobo <- a6te[a6te$NAME == "HOBO"]
a6hopper <- a6te[a6te$NAME == "HOPPER"]


hist(width(a6hopper)[width(a6hopper) > 1400], breaks = 50)
hist(width(a6hopper), breaks = 50)

cantonte <- import.gff("/Users/artemilin/Work/projects/DrosoTE/popTE/earlgrey/dm6_liftover_earlgrey/eg_final/A1-CantonS_pacbio.filteredRepeats.gff")

cahopper <- cantonte[cantonte$NAME == "HOPPER"]

hist(width(cahopper), breaks = 50)
hist(width(a6hobo))

a6fasta <- "/Users/artemilin/Work/projects/DrosoTE/popTE/earlgrey/dsrp_genomes/A6-Wild5B_pacbio.fa"
cantonfa <- "/Users/artemilin/Work/projects/DrosoTE/popTE/earlgrey/dsrp_genomes/A1-CantonS_pacbio.fa"

#### Sequence of A6 hoppers ####

genomio <- FaFile(a6fasta)
open(genomio)
a6hopper.fa <- getSeq(genomio, a6hopper)
names(a6hopper.fa) <- a6hopper$ID
a6hopper.fa <- a6hopper.fa[width(a6hopper.fa) > 1350]

longhop <- readDNAStringSet("~/Work/projects/DrosoTE/popTE/earlgrey/a6_full_hopper_putat.fa")
longhop_1 <- subseq(longhop, 1, 432)
longhop_2 <- subseq(longhop, 1813, 1879)
longhop_3 <- subseq(longhop, 1886, width(longhop))
pseudoshorthop <- DNAStringSet(paste0(longhop_1, longhop_2, longhop_3))
names(pseudoshorthop) <- "A6_wild5b_pseudoshort_hopper"
middlehop <- subseq(longhop, 432, 1813)
names(middlehop) <- "Long_hopper_centre"

writeXStringSet(pseudoshorthop, "~/Work/projects/drosote/Hopper_A6/a6_pseudoshorthopper.fa")
writeXStringSet(middlehop, "~/Work/projects/drosote/Hopper_A6/hopper_middle.fa")
a6hopper.fa <- c(a6hopper.fa, longhop)

# make sure long copy has a recognizable unique name
names(a6hopper.fa)[length(a6hopper.fa)] <- "long_hopper"

# align 
hopa6.aln <- AlignSeqs(a6hopper.fa)

# distance matrix
dm <- DistanceMatrix(hopa6.aln, type = "dist")

# order by distance to long copy 

long_name <- "long_hopper"

dist_to_long <- as.matrix(dm)[, long_name]
dist_to_long <- sort(dist_to_long, decreasing = FALSE)

dist_to_long

hopa6_aln_ord <- hopa6.aln[names(dist_to_long)]
hopa6_aln_ord.rev <- reverseComplement(hopa6_aln_ord)

BrowseSeqs(hopa6_aln_ord)

subseq(hopa6_aln_ord.rev, 1039, 2420)

#### Get coverage of short copies per long copy ####

smooth_window <- 50

###### alignment matrix ######

aln_mat <- do.call(
  rbind,
  strsplit(as.character(hopa6_aln_ord.rev), split = "")
)

rownames(aln_mat) <- names(hopa6_aln_ord.rev)

ref <- aln_mat["long_hopper", ]

###### reference-coordinate conversion ######
# alignment columns where reference has a gap are removed 

ref_is_base <- ref != "-"
ref_pos <- cumsum(ref_is_base)
ref_pos[!ref_is_base] <- NA_integer_

keep_cols <- !is.na(ref_pos)

aln_ref <- aln_mat[, keep_cols, drop = FALSE]
ref_pos_kept <- ref_pos[keep_cols]

# query sequences only

query_mat <- aln_ref[rownames(aln_ref) != ref_name, , drop = FALSE]

###### coverage along reference ######

coverage_df <- tibble(
  ref = ref_name,
  pos = ref_pos_kept,
  coverage = colSums(query_mat != "-"),
  n_sequences = nrow(query_mat),
  frac_covered = coverage / n_sequences
) |>
  arrange(pos) |>
  mutate(
    coverage_smooth = zoo::rollmean(
      coverage,
      k = 50,
      fill = NA,
      align = "center"
    ),
    frac_covered_smooth = zoo::rollmean(
      frac_covered,
      k = 50,
      fill = NA,
      align = "center"
    )
  )

###### plot coverage ######


ref_len <- max(coverage_df$pos)
# x_marks <- c(1, 438, 1819, ref_len)
x_marks <- c(1, 1040, 2419, ref_len)

p_cov <- ggplot(coverage_df, aes(x = pos, y = coverage)) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = 2) +
  
  geom_line(linewidth = 0.35) +
  
  # vertical lines at reference start and end
  geom_segment(
    data = tibble(x = c(1, ref_len)),
    aes(
      x = x,
      xend = x,
      y = 0,
      yend = max(coverage_df$coverage, na.rm = TRUE)
    ),
    inherit.aes = FALSE,
    linewidth = 0.4
  ) +
  
  scale_x_continuous(
    breaks = x_marks,
    labels = x_marks,
    limits = c(1, ref_len),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  
  theme_classic() +
  labs(
    x = "Position along reference, nt",
    y = "Number of sequences with nucleotide present",
    title = "MSA coverage of non-autonomours Hopper insertions\nalong autonomous Hopper sequence in A6 genome",
    subtitle = paste0("n query sequences = ", nrow(query_mat))
  )

p_cov

ggsave(
  "outputs/plots/a6_shopp_msa_cov.pdf",
  p_cov,
  width = 6,
  height = 3.5
)

source("functions/plot_msa_coverage.R")
a6.shopp.msa <- plot_msa_coverage(
  alignment = hopa6_aln_ord.rev,
  ref_name = "long_hopper",
  fill_color = "#4E79A7",
  title = paste0(
    "MSA coverage of Kodama insertions\n",
    "along the Nozomi sequence in the A6 genome"
  )
)

a6.shopp.msa$plot

ggplot2::ggsave(
  "outputs/plots/a6_kodanozo_msa_cov.pdf",
  a6.shopp.msa$plot+xlab("Position along Nozomi reference"),
  width = 6,
  height = 3.5
)

hopper_cov$consensus_deletions
hopper_cov$x_marks

#### Show the identity of long hoppers ####

# load long hoppers
load("outputs/RData/A6_long_hoppers.RData")
# a6.longhop.new <- a6lh.seq[-(2:4)]
# align
hopa6.aln <- AlignSeqs(a6lh.seq)


###### alignment matrix ######

aln_mat <- do.call(
  rbind,
  strsplit(as.character(hopa6.aln), split = "")
)

rownames(aln_mat) <- as.character(a6_lh_true)

# Remove alignment columns consisting entirely of gaps
aln_mat <- aln_mat[
  ,
  colSums(aln_mat != "-") > 0,
  drop = FALSE
]

copy_names <- rownames(aln_mat)
n_copies <- nrow(aln_mat)

###### pairwise nucleotide identity ######

identity_mat <- matrix(
  NA_real_,
  nrow = n_copies,
  ncol = n_copies,
  dimnames = list(copy_names, copy_names)
)

for (i in seq_len(n_copies)) {
  for (j in seq_len(n_copies)) {
    
    identity_mat[i, j] <- mean(
      aln_mat[i, ] == aln_mat[j, ]
    ) * 100
  }
}

###### clustering based on sequence divergence ######

identity_dist <- as.dist(100 - identity_mat)

row_clustering <- hclust(
  identity_dist,
  method = "average"
)

col_clustering <- hclust(
  identity_dist,
  method = "average"
)

###### heatmap ######


library(pheatmap)
library(grid)
library(gtable)

rownames(identity_mat) <- c(
  "chr3L_Nozomi_1",
  "chr3L_fixed_1",
  "chr3L_fixed_2",
  "chr3L_fixed_3",
  "chr3R_Nozomi_1",
  "chr3R_Nozomi_2",
  "chr2L_Nozomi_1",
  "chr2L_Nozomi_2",
  "chr2L_Nozomi_3",
  "chr2R_Nozomi_1",
  "chr2R_Nozomi_2",
  "chrX_Nozomi_1"
)

colnames(identity_mat) <- rownames(identity_mat)

ph <- pheatmap(
  identity_mat,
  
  cluster_rows = TRUE,
  cluster_cols = row_clustering,
  treeheight_row = 0,
  treeheight_col = 22,
  
  display_numbers = matrix(
    sprintf("%.2f%%", identity_mat),
    nrow = nrow(identity_mat),
    dimnames = dimnames(identity_mat)
  ),
  
  number_color = "black",
  fontsize_number = 5.5,
  
  color = colorRampPalette(
    c("white", "mistyrose", "red")
  )(100),
  breaks = seq(95, 100, length.out = 101),
  
  border_color = "white",
  legend = FALSE,
  
  angle_col = 45,
  fontsize = 9,
  fontsize_row = 7,
  fontsize_col = 7,
  
  main = "Pairwise nucleotide identity of Nozomi insertions\nin A6 genome",
  fontsize_main = 16,
  
  silent = TRUE
)

pdf(
  "outputs/plots/A6_long_hoppers_heat_id.pdf",
  width = 6,
  height = 6,
  useDingbats = FALSE
)

grid.newpage()

# Real 12-mm margins on both left and right
pushViewport(
  viewport(
    x = 0.5,
    y = 0.5,
    width = unit(1, "npc") - unit(18, "mm"),
    height = unit(1, "npc") - unit(6, "mm")
  )
)

grid.draw(ph$gtable)

popViewport()
dev.off()


#### optional presence/absence heatmap ####

presence_df <- as_tibble(query_mat, rownames = "sequence") |>
  pivot_longer(
    cols = -sequence,
    names_to = "aln_col",
    values_to = "base"
  ) |>
  group_by(sequence) |>
  mutate(pos = ref_pos_kept) |>
  ungroup() |>
  mutate(covered = base != "-") |>
  select(sequence, pos, base, covered)

p_heat <- ggplot(presence_df, aes(x = pos, y = sequence, fill = covered)) +
  geom_tile() +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 6),
    legend.position = "none"
  ) +
  labs(
    x = paste0("Position along reference: ", ref_name),
    y = NULL,
    title = "Presence/absence of aligned sequence along reference"
  )

ggsave(
  paste0(out_prefix, ".presence_heatmap.pdf"),
  p_heat,
  width = 12,
  height = max(4, 0.18 * length(unique(presence_df$sequence)))
)

coverage_df

###### TSDs of A6 hoppers ######

a6hopper_tsds <-  get_tsd_motif(
  gr = a6hopper,
  fasta = a6fasta,
  tsd_width = 5,
  use = "consensus",
  min_identity = 1
)

a6hopper_tsds$tsd_table %>% group_by(left_tsd) %>% 
  summarize(n = n()) %>% View
a6hopper_tsds$consensus
a6hopper_tsds$frequency_matrix
a6hopper_tsds$pwm

#### Sequence and alignment of A6 hobos ####
a6hobo.bed <- GenomicRanges::reduce(a6hobo, min.gapwidth = 10)
hist(width(a6hobo.bed), breaks = 50)
a6genomio <- FaFile(a6fasta)
open(a6genomio)
a6hobo.fa <- getSeq(a6genomio, a6hobo.bed)
a6hobo.fa <- a6hobo.fa[startsWith(as.character(a6hobo.fa), "CAGAGAA")]
a6hobo.fa <- a6hobo.fa[endsWith(as.character(a6hobo.fa), "CTCTG")]
hist(width(a6hobo.fa), breaks = 60)
names(a6hobo.fa) <- a6hobo$ID


hoboa6.aln <- AlignSeqs(a6hobo.fa)
BrowseSeqs(hoboa6.aln)

dm <- DistanceMatrix(hoboa6.aln, type = "dist")
hc <- hclust(dm, method = "average")

# reordered alignment
hobo_aln_ord <- hoboa6.aln[hc$order]

BrowseSeqs(hobo_aln_ord)

###### Get TSDs of hobos ######

a6hobo.bed <- a6hobo.bed[width(a6hobo.bed) != 2721]

a6hobo.bed2 <- resize(a6hobo.bed, width = width(a6hobo.bed) + 10, fix = "end")
a6hobo.bed2 <- resize(a6hobo.bed2, width = width(a6hobo.bed2) + 10, fix = "start")

a6hobo_tsd <- get_tsd_motif(
  gr = a6hobo.bed,
  fasta = a6fasta,
  tsd_width = 8,
  use = "consensus",
  min_identity = 1
)

View(a6hobo_tsd$tsd_table)

#### Hoppers in A1-CantonS ####

genomio <- FaFile(cantonfa)
open(genomio)
a1hopper.fa <- getSeq(genomio, cahopper)
names(a1hopper.fa) <- cahopper$ID
a1hopper.fa <- a1hopper.fa[width(a1hopper.fa) > 1300]
a1hopper.fa <- c(a1hopper.fa, longhop)

hopa1.aln <- AlignSeqs(a1hopper.fa)
BrowseSeqs(hopa1.aln)

#### Bactrocera dorsalis hopper ####

bdhop.bed <- import.bed("~/Work/projects/DrosoTE/Hopper_A6/teaid/bacdor_hopper/bdor_hopper.bed")
bdhop.bed <- GenomicRanges::reduce(bdhop.bed)
hist(width(bdhop.bed))
bdorfa <- "~/Work/projects/DrosoTE/Hopper_A6/genomes/Bactrocera_dorsalis_pacbio.fa"
bdorgenome <- FaFile(bdorfa)
open(bdorgenome)

###### TSDs of Bdor hopper ######
bdhop.bed2 <- GenomicRanges::reduce(bdhop.bed, min.gapwidth = 300)
bdhop.bed2 <- resize(bdhop.bed2, width = width(bdhop.bed2) - 2, fix = "end")
bdhop.bed2 <- resize(bdhop.bed2, width = width(bdhop.bed2) - 2, fix = "start")
hist(width(bdhop.bed2))

bdorhop.fa <- getSeq(bdorgenome, bdhop.bed2)
bdorhop.aln <- AlignSeqs(bdorhop.fa)
bdorhop.aln <- AlignSeqs(bdorhop.fa[c(1:8, 11, 12, 17, 20, 22, 24, 26, 29, 31, 34, 36, 40, 44, 45)])
BrowseSeqs(bdorhop.aln)


bdhop_tsd <- get_tsd_motif(
  gr = bdhop.bed2,
  fasta = bdorfa,
  tsd_width = 2,
  use = "consensus",
  min_identity = 0
)

bdhop_tsd$tsd_table %>% View()

#### Drosophila subpulchrella hopper ####

dsub.bed <- import.bed("~/Work/projects/DrosoTE/Hopper_A6/dsub_nakanokori.bed")
dsub.bed <- GenomicRanges::reduce(dsub.bed)
dsub.bed <- dsub.bed[width(dsub.bed) > 300]
hist(width(dsub.bed))
dsubfa <- "~/Work/projects/DrosoTE/Hopper_A6/genomes/Dsub_genome_v2.fa"
dsubgenome <- FaFile(dsubfa)
open(dsubgenome)

dsubhop.fa <- getSeq(dsubgenome, dsub.bed)
dsubhop.aln <- AlignSeqs(dsubhop.fa)

BrowseSeqs(dsubhop.aln)
nchar("tagtgttgggaactatcga")

###### TSDs of Dsub hopper ######
dsub.bed2 <- resize(dsub.bed, width = width(dsub.bed) + 10, fix = "end")
dsub.bed2 <- resize(dsub.bed2, width = width(dsub.bed2) + 10, fix = "start")
dsubhop2.fa <- getSeq(dsubgenome, dsub.bed2)
dsubhop2.aln <- AlignSeqs(dsubhop2.fa)

dsub_tsd <- get_tsd_motif(
  gr = dsub.bed,
  fasta = dsubfa,
  tsd_width = 5,
  use = "consensus",
  min_identity = 1
)

dsub_tsd$tsd_table %>% View()

#### Drosophila bipectinata hopper - transib6 ####

dbip.bed <- import.bed("~/Work/projects/DrosoTE/Hopper_A6/dbip_transib6.bed")
dbip.bed <- GenomicRanges::reduce(dbip.bed, min.gapwidth = 300)
dbip.bed$id <- paste("transib6", 1:length(dbip.bed), sep = "_")
names(dbip.bed) <- dbip.bed$id
dbip.bed <- dbip.bed[width(dbip.bed) > 300]
hist(width(dbip.bed), breaks = 50)
dbipfa <- "~/Work/projects/DrosoTE/fifteen_genomes/Dbip/Dbip.15g.fasta"
dbipgenome <- FaFile(dbipfa)
open(dbipgenome)

dbiphop.fa <- getSeq(dbipgenome, dbip.bed)
names(dbiphop.fa) <- names(dbip.bed)
dbiphop.fa <- dbiphop.fa[startsWith(as.character(dbiphop.fa), "GCACTAT")]
dbiphop.fa <- dbiphop.fa[endsWith(as.character(dbiphop.fa), "ATAGTGC")]
dbiphop.aln <- AlignSeqs(dbiphop.fa, gapExtension = -4, gapOpening = -7,
                         refinements = 30)

BrowseSeqs(dbiphop.aln)
nchar("tagtgttgggaactatcga")

###### TSDs of dbip hopper ######
dbip.bed2 <- dbip.bed[names(dbiphop.fa)]
dbip.bed3 <- resize(dbip.bed2, width = width(dbip.bed2) + 10, fix = "end")
dbip.bed3 <- resize(dbip.bed3, width = width(dbip.bed3) + 10, fix = "start")
dbiphop2.fa <- getSeq(dbipgenome, dbip.bed3[-c(87,94)])
# dbiphop2.aln <- AlignSeqs(dbiphop2.fa)



dbip_tsd <- get_tsd_motif(
  gr = dbip.bed,
  fasta = dbipfa,
  tsd_width = 3,
  use = "consensus",
  min_identity = 1
)

dbip_tsd$tsd_table %>% View()

###### Dbip Transib-6 MSA Profile, central deletion ######
source("functions/plot_msa_coverage.R")
dbiphop3.fa <- getSeq(dbipgenome, dbip_tsd$insertions)
names(dbiphop3.fa) <- paste("t6", 1:length(dbiphop3.fa), sep = "_")
summary(width(dbiphop3.fa))
hist(width(dbiphop3.fa), breaks = 50)


dbpts6.aln <- AlignSeqs(dbiphop3.fa[width(dbiphop3.fa) > 800][-26])
dbpst6.dm <- DistanceMatrix(dbpts6.aln, type = 'dist')
dbpst6.hc <- hclust(dbpst6.dm, method = 'average')
BrowseSeqs(dbpts6.aln[dbpst6.hc$order])

# dbiphop.msa <- plot_msa_coverage(dbpts6.aln,
#                                  ref_name = "t6_1",
#                                  title = "MSA coverage of Transib-6 insertions\nalong autonomous Transib-6 sequence in D. bipectinata genome",
#                                  min_deletion_width = 50,
#                                  shade_deletions = F)
dbiphop.msa <- plot_msa_coverage(
  alignment = dbpts6.aln,
  ref_name = "t6_1",
  fill_color = "#8E44AD",
  title = paste0(
    "MSA coverage of Transib-6 insertions\n",
    "along reference in D. bipectinata genome"
  )
)
dbiphop.msa$plot+xlab("Position along Transib-6 reference")
dbiphop.msa$individual_deletions

ggsave(
  "outputs/plots/dbip_t6_msa_cov.pdf",
  dbiphop.msa$plot+xlab("Position along Transib-6 reference"),
  width = 6,
  height = 3.5
)

#### All dspr lines ####

egdir <- "/Users/artemilin/Work/projects/DrosoTE/popTE/earlgrey/dm6_liftover_earlgrey/eg_final/"
alltes <- lapply(dir(egdir, full.names = T), function(fil){
  import.gff(fil)
})
allhop <- lapply(alltes, function(gff){
  gff[gff$NAME == "HOPPER"]
})
names(allhop) <- sub("_pacbio.filteredRepeats.gff", "", dir(egdir))
hop_widths <- bind_rows(lapply(names(allhop), function(x) {
  data.frame(
    element = x,
    width = width(allhop[[x]])
  )
}))

ggplot(hop_widths, aes(x = width)) +
  geom_histogram(bins = 50, color = "black", fill = "grey70") +
  facet_wrap(~ element, ncol = 4, nrow = 4, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Insertion width",
    y = "Count"
  )

#### DrosEU hoppers length ####
droseueg <- "/Users/artemilin/Work/projects/DrosoTE/popTE/DrosEU/earlgrey"
deu.alltes <- lapply(dir(droseueg, full.names = T), function(fil){
  import.gff(fil)
})
deu.allhop <- lapply(deu.alltes, function(gff){
  gff[gff$NAME == "HOPPER"]
})

names(deu.allhop) <- sub("_EarlGrey.filteredRepeats.gff", "", dir(droseueg))
deu.hop_widths <- bind_rows(lapply(names(deu.allhop), function(x) {
  data.frame(
    element = x,
    width = width(deu.allhop[[x]])
  )
}))

ggplot(deu.hop_widths, aes(x = width)) +
  geom_histogram(bins = 50, color = "black", fill = "grey70") +
  facet_wrap(~ element, ncol = 4, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Insertion width",
    y = "Count"
  )

###### Hopper2 ######
deu.allhop2 <- lapply(deu.alltes, function(gff){
  gff[gff$NAME == "HOPPER2"]
})

names(deu.allhop2) <- sub("_EarlGrey.filteredRepeats.gff", "", dir(droseueg))
deu.hop2_widths <- bind_rows(lapply(names(deu.allhop2), function(x) {
  data.frame(
    element = x,
    width = width(deu.allhop2[[x]])
  )
}))

ggplot(deu.hop2_widths, aes(x = width)) +
  geom_histogram(bins = 50, color = "black", fill = "grey70") +
  facet_wrap(~ element, ncol = 4, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Insertion width",
    y = "Count"
  )

#### Scat ster hopper ####

schop <- readDNAStringSet("~/Work/projects/DrosoTE/Hopper_A6/scatster_hopper.fa")

schopac <- substring(schop, 950, 2266) %>% DNAString()
