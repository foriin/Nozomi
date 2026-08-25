#### hop_dspr_droseu_insertions_prep.R ####
# Script for preparing the bed files with hopper insertions #

library(GenomicRanges)
library(rtracklayer)
library(Biostrings)
library(Rsamtools)
library(ggplot2)

#### Prepare all 1.4 kb Hopper insertions ####

beddir <- "/Users/artemilin/Work/projects/DrosoTE/Hopper_A6/LIFT/hopper_bed"

shortbeds <- dir(beddir, pattern = "short.bed", full.names = T)
names(shortbeds) <- sub("\\.hopper_short\\.bed", "",
                        dir(beddir, pattern = "short.bed"))

hopsh <- lapply(names(shortbeds), function(bedd){
  gra <- import.bed(shortbeds[[bedd]])
  gra <- GenomicRanges::reduce(gra)
  gra <- gra[width(gra) > 1400]
  gra$origin <- paste(bedd, 1:length(gra), sep = "_")
  names(gra) <- gra$origin
  gra
})

names(hopsh) <- names(shortbeds)

hopsh.all <- lapply(names(shortbeds), function(bedd){
  gra <- import.bed(shortbeds[[bedd]])
  gra <- GenomicRanges::reduce(gra)
  gra$origin <- paste(bedd, 1:length(gra), sep = "_")
  names(gra) <- gra$origin
  gra
})

names(hopsh.all) <- names(shortbeds)



beddir_out <- "/Users/artemilin/Work/projects/DrosoTE/Hopper_A6/LIFT/hopper_bed_prepped/"
dir.create(beddir_out)
for (i in names(hopsh)){
  export.bed(hopsh[[i]],
             paste0(beddir_out, i, "_hopper_1.4kb.bed"))
}

###### Get all the sequences ######

# Prepare genomes
dspr_fa_dir <- "/Users/artemilin/Work/projects/DrosoTE/popTE/earlgrey/dsrp_genomes"
dspr_fa <- dir(dspr_fa_dir, pattern = '.fa$', full.names = T)
names(dspr_fa) <- sub("_pacbio.fa", "", dir(dspr_fa_dir, pattern = ".fa$"))

droseu_fa_dir <- "/Users/artemilin/Work/projects/DrosoTE/popTE/earlgrey/droseu_genomes"
droseu_fa <- dir(droseu_fa_dir, pattern = '.fasta$', full.names = T) 
names(droseu_fa) <- sub("\\.fasta", "", dir(droseu_fa_dir, pattern = '.fasta$'))

genomes_all <- c(dspr_fa, droseu_fa)
all(names(hopsh) %in% names(genomes_all))

hopsh.fa <- lapply(names(hopsh), function(strain){
  gr <- hopsh[[strain]]
  genome <- FaFile(genomes_all[[strain]])
  open(genome)
  dna <- getSeq(genome, gr)
  close(genome)
  names(dna) <- gr$origin
  dna
})

hopsh.fa <- DNAStringSetList(hopsh.fa) %>% unlist

# Without size filtering
hopsh.all.fa <- lapply(names(hopsh.all), function(strain){
  gr <- hopsh.all[[strain]]
  genome <- FaFile(genomes_all[[strain]])
  open(genome)
  dna <- getSeq(genome, gr)
  close(genome)
  names(dna) <- gr$origin
  dna
})

hopsh.all.fa <- DNAStringSetList(hopsh.all.fa) %>% unlist
save(hopsh.all.fa, file = "outputs/RData/shorthopper_all_dspr_droseu.RData")

unique(hopsh.fa)
hopsh.aln <- AlignSeqs(hopsh.fa)
dm <- DistanceMatrix(hopsh.aln, type = 'dist')
hc <- hclust(dm)

BrowseSeqs(hopsh.aln[hc$order])


###### get TSDs ######
source("get_TSD.R")
get_tsd_motif()
hopsh.tsd <- lapply(names(hopsh), function(strain){
  gr <- hopsh[[strain]]
  print(length(gr))
  genome <- FaFile(genomes_all[[strain]])
  open(genome)
  tsds <- get_tsd_both_sides(gr, genome)
  close(genome)
  tsds
})
hopsh.tsd <- do.call(rbind, hopsh.tsd) 

save(hopsh, hopsh.fa, hopsh.tsd, file = "outputs/RData/short_hoppers_origloc_seq.RData")

#### Prepare all 2.8 kb Hopper insertions ####

longbeds <- dir(beddir, pattern = "long.bed", full.names = T)
names(longbeds) <- sub("\\.hopper_long\\.bed", "",
                        dir(beddir, pattern = "long.bed"))

hopl <- lapply(names(longbeds), function(bedd){
  gra <- import.bed(longbeds[[bedd]])
  gra <- GenomicRanges::reduce(gra)
  gra <- gra[width(gra) > 2400]
  if (length(gra) == 0){
    return(NULL)
  }
  gra$origin <- paste(bedd, "longhopper", 1:length(gra), sep = "_")
  names(gra) <- gra$origin
  gra
})

names(hopl) <- names(longbeds)

hopl <- hopl[sapply(hopl, length) > 0]

for (i in names(hopl)){
  export.bed(hopl[[i]],
             paste0(beddir_out, i, "_hopper_long.bed"))
}


###### Get sequences ######
all(names(hopl) %in% names(genomes_all))

hopl.fa <- lapply(names(hopl), function(strain){
  gr <- hopl[[strain]]
  genome <- FaFile(genomes_all[[strain]])
  open(genome)
  dna <- getSeq(genome, gr)
  close(genome)
  names(dna) <- gr$origin
  dna
})

hopl.fa <- DNAStringSetList(hopl.fa) %>% unlist
unique(hopl.fa)
hopl.aln <- AlignSeqs(hopl.fa)

hopl.aln <- hopl.aln[names(allhopsl.dm6)]
dm <- DistanceMatrix(hopl.aln, type = 'dist')
hc <- hclust(dm, method = 'average')
BrowseSeqs(hopl.aln[hc$order])

names(hopl.aln) <- paste(names(hopl.aln), allhopsl.dm6[names(hopl.aln)]$group, sep = '_')

dna <- as.DNAbin(hopl.aln)

dm <- dist.dna(dna, model = "raw", pairwise.deletion = TRUE)
tree <- nj(dm)

plot(tree, cex = 0.4)

dem1 <- subseq(hopl.fa['DE_Mun_15_20_longhopper_3'], 1, 1437)
dem2 <- subseq(hopl.fa['DE_Mun_15_20_longhopper_3'], 1438, 2873)
names(dem2) <- "DE_Mun_15_20_longhopper_4"
hopl.fa <- c(hopl.fa[names(hopl.fa) != 'DE_Mun_15_20_longhopper_3'], dem1, dem2)
hopl.aln <- AlignSeqs(hopl.fa)
dm <- DistanceMatrix(hopl.aln, type = 'dist')
hc <- hclust(dm, method = 'average')
BrowseSeqs(hopl.aln[hc$order])


hopall.fa <- c(hopsh.fa, dem1, dem2, longhop)

hopal.aln <- AlignSeqs(hopall.fa)
dna <- as.DNAbin(hopal.aln)

dm <- dist.dna(dna, model = "raw", pairwise.deletion = TRUE)
tree <- nj(dm)

plot(tree, cex = 0.4)
dir.create("outputs/phylo")
write.tree(tree, file = "outputs/phylo/hopper_tree.nwk")



