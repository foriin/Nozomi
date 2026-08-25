#### ~/R/projects/HOPPER/hopper_short_lift2dm6_consloc.R ####
# Script for checking location conservation of short hoppers
library(GenomicRanges)
library(rtracklayer)
library(Biostrings)
library(Rsamtools)
library(ggplot2)

beddir <- "/Users/artemilin/Work/projects/DrosoTE/Hopper_A6/LIFT/hopper_bed_dm6"
hopldm6 <- dir(beddir, pattern = "sorted.bed", full.names = T)
names(hopldm6) <- sub("_hopper_1.4kb_dm6.sorted.bed", "", dir(beddir, pattern = 'sorted.bed'))

allh1.4.dm6 <- lapply(hopldm6, import.bed)

allh1.4.dm6 <- unlist(GRangesList(allh1.4.dm6)) %>% sort()
names(allh1.4.dm6) <- allh1.4.dm6$name
as.data.frame(allh1.4.dm6) %>% View()
allh1.4.dm6$id <- sub("_t[35]", "", allh1.4.dm6$name)
allh1.4.dm6$strain <- sub("_[0-9]+$", "", allh1.4.dm6$id)

allh1.4.red <- GenomicRanges::reduce(allh1.4.dm6, with.revmap = T)

allh1.4.red$revmap_l <- sapply(allh1.4.red$revmap, length)
ahred.df <- as.data.frame(sort(allh1.4.red))
ahred.df$group <- 1:nrow(ahred.df)

kekw <- apply(ahred.df, 1, function(x){
  grr <- allh1.4.dm6[x$revmap]
  grr$group <- x$group
  grr$group_len <- x$revmap_l
  grr
  
})

allshops.dm6 <- unlist(GRangesList(kekw))
names(allshops.dm6) <- sub("_t[35]", "", names(allshops.dm6))


#### Load sequence data ####
load("outputs/RData/short_hoppers_origloc_seq.RData", verb = T)

shops.fa <- hopsh.fa[names(allshops.dm6)]
grep("_t", names(allshops.dm6), value = T)
grep("_t3", names(hopsh.fa), value = T)


longhop <- readDNAStringSet("~/Work/projects/DrosoTE/popTE/earlgrey/a6_full_hopper_putat.fa")
shops.fa <- c(shops.fa, longhop)

shop.aln <- AlignSeqs(shops.fa)

#### distance matrix ####

dm <- DistanceMatrix(shop.aln, type = "dist")

#### order by distance to long copy ####

long_name <- "A6_hopper_long"

dist_to_long <- as.matrix(dm)[, long_name]
dist_to_long <- sort(dist_to_long, decreasing = FALSE)

dist_to_long

shop_aln_ord <- shop.aln[names(dist_to_long)]

BrowseSeqs(shop_aln_ord)

###### Add group data to sequences ######
shops.fa@elementMetadata <- DataFrame(group = allshops.dm6[names(shops.fa)]$group,
                                      strain = allshops.dm6[names(shops.fa)]$strain,
                                      loc = as.character(allshops.dm6[names(shops.fa)]))

###### Look at MSAs for groups ######

enframe(shops.fa@elementMetadata$group) %>% 
  group_by(value) %>% summarize(n = n()) %>% 
  arrange(-n)

shop268.aln <- AlignSeqs(shops.fa[shops.fa@elementMetadata$group == 268])
BrowseSeqs(shop268.aln)

shop412.aln <- AlignSeqs(shops.fa[shops.fa@elementMetadata$group == 412])
BrowseSeqs(shop412.aln)

shop250.aln <- AlignSeqs(shops.fa[shops.fa@elementMetadata$group == 250])
mcols(shop250.aln) <- mcols(shops.fa[shops.fa@elementMetadata$group == 250])
names(shop250.aln) <- paste(names(shop250.aln), mcols(shop250.aln)$group,
                            mcols(shop250.aln)$strict_haplotype)
BrowseSeqs(shop250.aln)

shop146.aln <- AlignSeqs(shops.fa[shops.fa@elementMetadata$group == 146])
BrowseSeqs(shop146.aln)

shop134.aln <- AlignSeqs(shops.fa[shops.fa@elementMetadata$group == 134])
mcols(shop134.aln) <- mcols(shops.fa[shops.fa@elementMetadata$group == 134])
names(shop134.aln) <- paste(names(shop134.aln), mcols(shop134.aln)$group,
                            mcols(shop134.aln)$strict_haplotype)
dm134 <- DistanceMatrix(shop134.aln, type = 'dist')
hc134 <- hclust(dm134, method = 'average')
BrowseSeqs(shop134.aln[hc134$order])

shop91.aln <- AlignSeqs(shops.fa[shops.fa@elementMetadata$group == 91])
mcols(shop91.aln) <- mcols(shops.fa[shops.fa@elementMetadata$group == 91])
names(shop91.aln) <- paste(names(shop91.aln), mcols(shop91.aln)$group,
                            mcols(shop91.aln)$strict_haplotype)
BrowseSeqs(shop91.aln[order(mcols(shop91.aln)$strict_haplotype)])

shop90.aln <- AlignSeqs(shops.fa[shops.fa@elementMetadata$group == 90])
mcols(shop90.aln) <- mcols(shops.fa[shops.fa@elementMetadata$group == 90])
names(shop90.aln) <- paste(names(shop90.aln), mcols(shop90.aln)$group,
                           mcols(shop90.aln)$strict_haplotype)
BrowseSeqs(shop90.aln[order(mcols(shop90.aln)$strict_haplotype)])

save(shops.fa, file = 'outputs/RData/short_hoppers_seq_dm6loc_group.RData')

#### Back to A6 - some stats ####

###### First, testing the hotspots model ######
N <- 450   # total unique short Hopper locations
A <- 149   # locations present in A6
B <- 306   # locations present in at least one non-A6 genome
x <- 5     # A6 locations also present in at least one non-A6 genome

expected_overlap <- A * B / N

pval <- phyper(
  q = x,
  m = B,
  n = N - B,
  k = A,
  lower.tail = TRUE
)

expected_overlap
pval

x <- 5
n <- 149

bt <- binom.test(x, n)

data.frame(
  n_A6_insertions = n,
  n_shared_with_other_genomes = x,
  shared_fraction = x / n,
  private_fraction = 1 - x / n,
  ci95_low = bt$conf.int[1],
  ci95_high = bt$conf.int[2]
)

binom.test(5, 149, p = 0.10, alternative = "less")
binom.test(5, 149, p = 0.20, alternative = "less")

