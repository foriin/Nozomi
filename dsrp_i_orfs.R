# BiocManager::install("ORFik")
# remotes::install_github("qsbase/qs")
library(ORFik)
library(Biostrings)
library(GenomicRanges)
library(GenomicFeatures)
library(ape)
library(phangorn)
library(patchwork)
library(gt)

# source("fasta/variability_profile_msa.R")
source("fasta/variability_profile_msa_ref.R")

#### Deal with ORFs in Copia ####
copia_orfs <- findORFs(COPIAsfa, minimumLength = 400)
names(copia_orfs) <- names(COPIAsfa)
corf_long <- unlist(copia_orfs) %>% as.data.frame() %>% 
  filter(width > 4000) %>% 
  mutate(strain = sub("_.*", "", names))

table(corf_long$strain)

# keep only ORFs from sequences that passed the width filter
copia_orfs_long <- copia_orfs[names(copia_orfs) %in% corf_long$names]

# extract ORF DNA sequences
copia_orf_seqs <- DNAStringSet(
  unlist(
    lapply(seq_along(copia_orfs_long), function(i) {
      ir <- copia_orfs_long[[i]]
      s  <- COPIAsfa[[ names(copia_orfs_long)[i] ]]
      
      if (length(ir) == 0) return(NULL)
      
      subseq(s, start = start(ir), end = end(ir))
    }),
    use.names = FALSE
  )
)

names(copia_orf_seqs) <- rep(names(copia_orfs_long), lengths(copia_orfs_long))

# translate to protein
copia_orf_pep <- translate(copia_orf_seqs)
names(copia_orf_pep) <- names(copia_orf_seqs)

copia_pep_aln <- AlignSeqs(copia_orf_pep[names(copia_orf_pep) != "A5-VAG1_COPIA_28222"])
BrowseSeqs(copia_pep_aln, highlight = 100)

copia_pep_a5 <- copia_orf_pep[grepl("A5|LTR", names(copia_orf_pep))]
copia_pep_a5.aln <- AlignSeqs(copia_pep_a5[names(copia_pep_a5) != "A5-VAG1_COPIA_28222"])
BrowseSeqs(copia_pep_a5.aln, highlight = 28)

copias_a5_long <- COPIAsfa[names(COPIAsfa) %in% names(copia_pep_a5)]
cla5.aln <- AlignSeqs(copias_a5_long[names(copias_a5_long) != "A5-VAG1_COPIA_28222"])
BrowseSeqs(cla5.aln, highlight = 28)

# ape wants matrix/alignment-like input; this usually works:
pep_phy <- as.phyDat(as.matrix(copia_pep_a5.aln), type = "AA")
dm <- dist.ml(pep_phy)
tr <- bionj(dm)

plot(tr, cex = 0.3, no.margin = TRUE)

# PCoA

pcoa <- cmdscale(dm, k = 2, eig = TRUE)

plot_df <- tibble(
  name = rownames(pcoa$points),
  PCo1 = pcoa$points[, 1],
  PCo2 = pcoa$points[, 2]
) %>%
  mutate(
    strain = sub("_.*", "", name),
    te = sub("^[^_]+_([^_]+)_.*$", "\\1", name),
    id = sub("^.*_([0-9]+)$", "\\1", name)
  )

var_exp <- 100 * pcoa$eig / sum(pcoa$eig[pcoa$eig > 0])

ggplot(plot_df, aes(PCo1, PCo2, color = strain, label = strain)) +
  geom_point(alpha = 0.7, size = 2) +
  ggrepel::geom_text_repel(size = 2.5) +
  theme_bw()

ggplot(plot_df, aes(PCo1, PCo2)) +
  geom_density_2d() +
  geom_point(aes(color = strain), size = 1, alpha = 0.5) +
  facet_wrap(~strain) +
  theme_bw()

##### Get back to the sequence level for "functional" copias ######
cla.seq <- COPIAsfa[names(COPIAsfa) %in% names(copia_orf_pep)]
cla.aln <- AlignSeqs(cla.seq)
BrowseSeqs(cla.aln, highlight = 101)

###### variability profile ######
cla.aln.mat <- aln_to_matrix(cla.aln)
cla.vp <- variability_profile(cla.aln.mat)
plot_profile(cla.vp)
#### Hopper ####
hopper_cands <- te.pca$HOPPER %>% 
  pull(seqname)

hoppersfa <- te.fa.c[names(te.fa.c) %in% hopper_cands]
hopperref <- tefa[names(tefa) == "hopper#DNA/CMC-Transib"]
hoppersfa <- c(hoppersfa, hopperref)
hoppersfa <- reverseComplement(hoppersfa)
writeXStringSet(hoppersfa[grepl("A6", names(hoppersfa))], "fasta/A6_hoppers.fa")
c("TAA", "TAG", "TGA")
hopper_orfs <- findORFs(hoppersfa[131], longestORF = T, minimumLength = 50,
                        stopCodon = "TGA")
names(hopper_orfs) <- names(hoppersfa)[as.integer(names(hopper_orfs))]
horf <- unlist(hopper_orfs) %>% as.data.frame() %>% 
  mutate(strain = sub("_.*", "", names))

table(corf_long$strain)

###### Align sequences ####
# Reorder based on similarity
hopper_aln <- AlignSeqs(hoppersfa)

# keep the last sequence fixed
last_seq <- hopper_aln[131]
rest_seqs <- hopper_aln[-131]
# Let's look only at A6
rest_seqs <- rest_seqs[grepl("A6", names(rest_seqs))]

# cluster the rest by similarity
dm <- DistanceMatrix(rest_seqs, type = "dist")
hc <- hclust(dm, method = "average")

# reordered alignment
hopper_aln_ord <- c(rest_seqs[hc$order], last_seq)

BrowseSeqs(hopper_aln_ord, highlight = 80)

# Plot tree
tr <- bionj(dm)
plot(tr, cex = 0.6, no.margin = TRUE)


###### Get long hopperinos ####

a6_long_hop <- import.bed("/Users/artemilin/Work/projects/DrosoTE/popTE/earlgrey/a6_long_hopper.bed")
a6_lh_true <- a6_long_hop[width(a6_long_hop) > 2500 & width(a6_long_hop) < 3500]
a6fasta <- "~/Work/projects/DrosoTE/popTE/earlgrey/dsrp_genomes/A6-Wild5B_pacbio.fa"
a6genome <- FaFile(a6fasta)
open(a6genome)
a6lh.seq <- getSeq(a6genome, a6_lh_true)
save(a6lh.seq, a6_lh_true, file = "outputs/RData/A6_long_hoppers.RData")


a6genome <- FaFile(genomefolder[grepl("A6", genomefolder)])
open(a6genome)

start(a6_lh_true)[2:4] <- start(a6_lh_true)[2:4] - 8

start(a6_lh_true) <- start(a6_lh_true) - 5
end(a6_lh_true) <- end(a6_lh_true) + 5

a6_lh_true.seq <- getSeq(a6genome, a6_lh_true)
names(a6_lh_true.seq) <- paste0("hopperL_", 1:12)

a6_lh_orfs <- findORFs(a6_lh_true.seq, minimumLength = 400)

names(a6_lh_orfs) <- names(a6_lh_true.seq)[as.integer(names(a6_lh_orfs))]
horf_long <- unlist(a6_lh_orfs) %>% as.data.frame() %>% 
  filter(width > 1500) %>% 
  mutate(strain = sub("_.*", "", names))

horfl_aln <- AlignSeqs(a6_lh_true.seq[horf_long$names])
BrowseSeqs(horfl_aln, highlight = 9)

p_msa <- plot_msa_variability(
  genome_fa = "~/Work/projects/DrosoTE/popTE/earlgrey/dsrp_genomes/A6-Wild5B_pacbio.fa",
  gr = a6_lh_true[c(1,5:12)],
  name_col = "name",
  msa_method = "ClustalOmega"
)

p_msa|wrap_table(as.data.frame(a6_lh_true[c(1,5:12)]) %>% dplyr::select(-name, -score),
                 panel = 'full', space = 'fixed')


###### Combine short and long into bed #####
a6hopperS <- teorig.te[teorig.te$strain == "A6-Wild5B" & teorig.te$NAME == "HOPPER"]
a6hopper_all <- c(a6hopperS, a6_lh_true)
export.bed(sort(a6hopper_all), "bed/A6_hoppers.bed")


#### Doc ####
doc_cands <- te.pca$DOC %>% 
  pull(seqname)

docsfa <- te.fa.c[names(te.fa.c) %in% doc_cands]
docref <- tefa[names(tefa) == "Doc#LINE/I-Jockey"]
docsfa <- c(docsfa, docref)

doc_orfs <- findORFs(docsfa, minimumLength = 400)
names(doc_orfs) <- names(docsfa)
dorf_long <- unlist(doc_orfs) %>% as.data.frame() %>% 
  filter(width > 1500) %>% 
  mutate(strain = sub("_.*", "", names))

table(dorf_long$strain)
# This thresholds were set up after manual expection
dorf_gag <- dorf_long %>% filter(width <= 1750)
dorf_pol <- dorf_long %>% filter(width > 2000)

###### GAG ######

# keep only ORFs from sequences that passed the gag width filter
doc_gags <- doc_orfs[names(doc_orfs) %in% dorf_gag$names]

# within those, keep only ORFs starting near the beginning
doc_gags <- endoapply(doc_gags, function(ir) ir[start(ir) < 500])

# drop empty elements
doc_gags <- doc_gags[lengths(doc_gags) > 0]

# extract ORF DNA sequences
doc_gag_seqs <- DNAStringSet(
  unlist(
    lapply(seq_along(doc_gags), function(i) {
      ir <- doc_gags[[i]]
      s <- docsfa[[names(doc_gags)[i]]]
      
      if (length(ir) == 0) return(NULL)
      
      subseq(s, start = start(ir), end = end(ir))
    }),
    use.names = FALSE
  )
)

names(doc_gag_seqs) <- rep(names(doc_gags), lengths(doc_gags))

# translate to protein
doc_orf_pep <- translate(doc_gag_seqs)
writeXStringSet(doc_orf_pep[names(doc_orf_pep) == "B2-CA1_DOC_119029"],
                "fasta/B2-CA1_DOC_119029_gag.fa")
writeXStringSet(doc_orf_pep[204],
                "fasta/Doc_ref_gag.fa")
names(doc_orf_pep) <- names(doc_gag_seqs)

doc_pep_aln <- AlignSeqs(doc_orf_pep)
BrowseSeqs(doc_pep_aln, highlight = 204)

# Reorder based on similarity

# keep the last sequence fixed
last_seq <- doc_pep_aln[204]
rest_seqs <- doc_pep_aln[-204]

# cluster the rest by similarity
dm <- DistanceMatrix(rest_seqs, type = "dist")
hc <- hclust(dm, method = "average")

# reordered alignment
doc_pep_aln_ord <- c(rest_seqs[hc$order], last_seq)

BrowseSeqs(doc_pep_aln_ord, highlight = 204)

#### Roo ####
names(te.fa) <- gna2
sapply(te.fa, function(x){
  roos <- x[grepl("ROO", names(x))]
  summary(width(roos))
})

grep("copia", names(tefa), ignore.case = T, value = T)

b2roos <- te.fa.c[grepl("_ROO", names(te.fa.c))]
rooref <- tefa[names(tefa) == "roo#LTR/Pao"]
b2roos <- c(b2roos, rooref)

roombas <- findORFs(b2roos, minimumLength = 400)
names(roombas) <- names(b2roos)[as.integer(names(roombas))]
roorf <- as.data.frame(unlist(roombas))
roo_orf <- roorf %>% filter(width > 7000) %>% pull(names)


b2roos <- b2roos[roo_orf]
b2roo.aln <- AlignSeqs(b2roos)

# keep the last sequence fixed
last_seq <- b2roo.aln[311]
rest_seqs <- b2roo.aln[-311]

# cluster the rest by similarity
dm <- DistanceMatrix(rest_seqs, type = "dist")
hc <- hclust(dm, method = "average")

# reordered alignment
b2roo_aln_ord <- c(rest_seqs[hc$order], last_seq)
BrowseSeqs(b2roo.aln, highlight = 311)

###### variability profile ######
b2roo.aln.mat <- aln_to_matrix(b2roo_aln_ord)
b2roo.vp <- variability_profile(b2roo.aln.mat)
plot_profile(b2roo.vp)

writeXStringSet(b2roos[311], "fasta/roo.fa")
anno_roo <- data.frame(
  feature = c("LTR5", "gag-pol-env", "LTR3"),
  start = c(1, 1224, 8665),
  end   = c(429, 8354, 9092)
)

anno_roo5utr <- data.frame(
  feature = c("LTR5", "DLS", "A-tract", "mu"),
  start = c(1, 815, 981, 1063),
  end   = c(429, 905, 996, 1158)
)

roo_res <- plot_msa_variability(
  seqs = b2roos,
  ref_name = "roo#LTR/Pao",
  anno = anno_roo,
  iterations = 10,
  processors = 8
)
print(roo_res$plot)
roo_res$profile

b2roos <- b2roos[names(b2roos) != "B6_ROO_5508"]
roos5utr <- subseq(b2roos, start = 1, end = 1500)
rooutr_res <- plot_msa_variability(
  seqs = roos5utr,
  ref_name = "roo#LTR/Pao",
  anno = anno_roo5utr,
  iterations = 10,
  processors = 8
)
print(rooutr_res$plot)


BrowseSeqs(b2roos[311])

#### 1731 ####
te.fa.all <- lapply(gna2, function(strain){
  print(strain)
  te.gr <- teorig.te.truec.f[teorig.te.truec.f$strain == strain]
  genomio <- FaFile(genomefolder[[strain]])
  open(genomio)
  indexFa(genomio)
  te.fa <- getSeq(genomio, te.gr)
  names(te.fa) <- te.gr$strainid
  te.fa
})

names(te.fa.all) <- gna2


sapply(te.fa.all, function(x){
  dm1731 <- x[grepl("1731", names(x))]
  length(dm1731)
})

te.fa.all.c <- do.call(c, te.fa.all)

sub("^.*_(.*)_.*$", "\\1", names(te.fa.all.c))  %>% sort() %>% table()

