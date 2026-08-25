
#### Start with DSPR ####
genome.dir <- "/Users/artemilin/Work/projects/DrosoTE/popTE/earlgrey/dsrp_genomes"
genomes <- list.files(genome.dir, pattern = ".fa$", full.names = T)

teaid_dir <- "/Users/artemilin/Work/projects/DrosoTE/Hopper_A6/teaid/DSRP_long_hopper"
teaid <- dir(teaid_dir, full.names = T)
names(teaid) <- dir(teaid_dir)

names(genomes) <- names(teaid)

hopps <- sapply(names(teaid), function(x) list.files(teaid[[x]], pattern = ".bed$", full.names = T))

hopps.gr <- lapply(names(hopps)[-6], function(strain){
  fp <- hopps[[strain]]
  geg <- import.bed(fp)
  geg$strain <- strain
  geg <- geg[width(geg) > 2500]
  geg
})
names(hopps.gr) <- names(hopps)[-6]

hopps.gr <- hopps.gr[sapply(hopps.gr, function(x) length(x) > 0)]

hopps.fa <- lapply(names(hopps.gr)[-6], function(strain){
  gr <- hopps.gr[[strain]]
  genome <- FaFile(genomes[[strain]])
  open(genome)
  seqs <- getSeq(genome, gr)
  names(seqs) <- paste(strain, as.character(gr), 1:length(seqs), sep = "_")
  close(genome)
  seqs
})

hopps.fa <- do.call(c, hopps.fa)
hopps.fix.fa <- hopps.fa[!grepl("chrX", names(hopps.fa))]
a6lhfix <- a6lh.seq[2:4]
names(a6lhfix) <- paste("A6_chr3L", 1:3, sep = "_")
a6lhfix <- reverseComplement(a6lhfix)
hopps.fix.fa <- c(hopps.fix.fa, a6lhfix)
save(a6lhfix, file = "outputs/RData/A6_Nozomi_chr3L_fix_seq.RData")

unique(sub("_.*", "", names(hopps.fix.fa)))

# hopps.fa <- reverseComplement(hopps.fa)

hopps.fix.aln <- AlignSeqs(hopps.fix.fa)

BrowseSeqs(hopps.fix.aln)

writeXStringSet(hopps.fa[1], "fasta/hopper_shorter.fa")

dm <- DistanceMatrix(hopps.fix.aln, type = 'dist')
median(dm)
mean(dm)
hc <- hclust(dm, method = "average")

hopps.aln <- hopps.aln[hc$order]

hopps.gr$B3
hopps.gr$A6
as.character(hopps.gr$A6)


#### Continue with DrosEU ####

deu.genome.dir <- "/Users/artemilin/Work/projects/DrosoTE/popTE/earlgrey/droseu_genomes"
deu.genomes <- list.files(deu.genome.dir, pattern = ".fasta$", full.names = T)

teaid_dir <- "/Users/artemilin/Work/projects/DrosoTE/Hopper_A6/teaid/DrosEU_longhop"
teaid <- dir(teaid_dir, full.names = T)
names(teaid) <- dir(teaid_dir)

names(deu.genomes) <- names(teaid)

hopps <- sapply(names(teaid), function(x) list.files(teaid[[x]], pattern = ".bed$", full.names = T))

hopps.gr <- lapply(names(hopps), function(strain){
  fp <- hopps[[strain]]
  geg <- import.bed(fp)
  # geg <- GenomicRanges::reduce(geg)
  geg$strain <- strain
  geg <- geg[width(geg) > 2500]
  geg
})
names(hopps.gr) <- names(hopps)

hopps.gr <- hopps.gr[sapply(hopps.gr, function(x) length(x) > 0)]

deu.hopps.fa <- lapply(names(hopps.gr), function(strain){
  gr <- hopps.gr[[strain]]
  genome <- FaFile(deu.genomes[[strain]])
  open(genome)
  seqs <- getSeq(genome, gr)
  names(seqs) <- paste(strain, as.character(gr), 1:length(seqs), sep = "_")
  close(genome)
  seqs
})

deu.hopps.fa <- do.call(c, deu.hopps.fa)
length(names(hopps.gr))

#### Combine dspr and droseu ####

lhops.all.fa <- c(hopps.fix.fa, deu.hopps.fa)
lhops.all.aln <- AlignSeqs(lhops.all.fa)
lha.dm <- DistanceMatrix(lhops.all.aln, type = 'dist')
lha.hc <- hclust(lha.dm, method = 'average')
BrowseSeqs(lhops.all.aln[lha.hc$order])
lhops.fix.fa <- lhops.all.fa[-c(47,44,38)]

lhops.fix.aln <- AlignSeqs(lhops.fix.fa)
lhf.dm <- DistanceMatrix(lhops.fix.aln, type = 'dist')
median(lhf.dm)
mean(lhf.dm)


#### group 134 ham dist ####

gr134 <- shops.fa[mcols(shops.fa)$group == 134]
gr134.aln <- AlignSeqs(gr134)
gr134.dm <- DistanceMatrix(gr134.aln, type = 'dist')
median(gr134.dm)
mean(gr134.dm)

#### group 255 ham dist ####

gr255 <- shops.fa[mcols(shops.fa)$group == 255]
gr255.aln <- AlignSeqs(gr255)
gr255.dm <- DistanceMatrix(gr255.aln, type = 'dist')
hc.gr255 <- hclust(gr255.dm, method = 'average')

mcols(gr255)[hc.gr255$order,] %>% as_tibble %>% View
median(gr255.dm)
mean(gr255.dm)

#### group 93 ham dist ####

gr93 <- shops.fa[mcols(shops.fa)$group == 93]
gr93.aln <- AlignSeqs(gr93)
gr93.dm <- DistanceMatrix(gr93.aln, type = 'dist')
hc.gr93 <- hclust(gr93.dm, method = 'average')

mcols(gr93)[hc.gr93$order,] %>% as_tibble %>% View
median(gr93.dm)
mean(gr93.dm)

#### Cluster of A6 long fixed hoppers ####
a6.sf.fa <- readDNAStringSet("inputs/fasta/a6_chr3l_short_fix_hopper.fa")
a6.lsf <- c(a6.sf.fa, reverseComplement(a6lhfix))
a6.lsf.aln <- AlignSeqs(a6.lsf)
BrowseSeqs(a6.lsf.aln)


