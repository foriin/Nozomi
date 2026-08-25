#### ~/R/projects/HOPPER/hopper_short_lift2dm6_consloc.R ####
# Script for checking location conservation of short hoppers
library(GenomicRanges)
library(rtracklayer)
library(Biostrings)
library(Rsamtools)
library(ggplot2)

beddir <- "/Users/artemilin/Work/projects/DrosoTE/Hopper_A6/LIFT/hopper_bed_dm6"
hopldm6 <- dir(beddir, pattern = "long_dm6.sorted.bed", full.names = T)
names(hopldm6) <- sub("_hopper_long_dm6.sorted.bed", "",
                      dir(beddir, pattern = 'long_dm6.sorted.bed'))

allhl.dm6 <- lapply(hopldm6, import.bed)

allhl.dm6 <- unlist(GRangesList(allhl.dm6)) %>% sort()
names(allhl.dm6) <- allhl.dm6$name
as.data.frame(allhl.dm6) %>% View()
allhl.dm6$id <- sub("_t[35]", "", allhl.dm6$name)
allhl.dm6$strain <- sub("_longhopper_[0-9]+$", "", allhl.dm6$id)

allhl.red <- GenomicRanges::reduce(allhl.dm6, with.revmap = T)

allhl.red$revmap_l <- sapply(allhl.red$revmap, length)
ahlred.df <- as.data.frame(sort(allhl.red))
ahlred.df$group <- 1:nrow(ahlred.df)

kekw <- apply(ahlred.df, 1, function(x){
  grr <- allhl.dm6[x$revmap]
  grr$group <- x$group
  grr$group_len <- x$revmap_l
  grr
  
})

allhopsl.dm6 <- unlist(GRangesList(kekw))
as.data.frame(allhopsl.dm6) %>% View
names(allhopsl.dm6) <- sub("_t[35]", "", names(allhopsl.dm6))


