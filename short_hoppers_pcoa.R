library(GenomicRanges)
library(Biostrings)
library(DECIPHER)
library(ggplot2)

## objects:
## allshops.dm6 = GRanges
## shops.fa     = DNAStringSet

stopifnot(!is.null(names(allshops.dm6)))
stopifnot(!is.null(names(shops.fa)))
stopifnot(!anyDuplicated(names(allshops.dm6)))
stopifnot(!anyDuplicated(names(shops.fa)))
stopifnot("group" %in% colnames(mcols(allshops.dm6)))

## keep only entries present in both objects
common_names <- intersect(names(allshops.dm6), names(shops.fa))

allshops.dm6.sub <- allshops.dm6[common_names]
shops.fa.sub <- shops.fa[common_names]

## reorder GRanges to match DNAStringSet exactly
allshops.dm6.sub <- allshops.dm6.sub[names(shops.fa.sub)]

stopifnot(identical(names(allshops.dm6.sub), names(shops.fa.sub)))

## align sequences
shops.aln <- AlignSeqs(shops.fa.sub)

## sequence distance matrix
dist.mat <- stringDist(shops.aln, method = "hamming")

## normalize by alignment length
dist.mat <- as.dist(dist.mat / width(shops.aln)[1])

## PCoA
pcoa <- cmdscale(dist.mat, k = 2, eig = TRUE)

## variance explained
eig <- pcoa$eig[pcoa$eig > 0]
var.explained <- eig / sum(eig)

pcoa.df <- data.frame(
  name = names(shops.aln),
  PCoA1 = pcoa$points[, 1],
  PCoA2 = pcoa$points[, 2],
  group = as.factor(mcols(allshops.dm6.sub)$group)
)

xlab <- paste0("PCoA1 (", round(var.explained[1] * 100, 1), "%)")
ylab <- paste0("PCoA2 (", round(var.explained[2] * 100, 1), "%)")

ggplot(pcoa.df %>% filter(PCoA1 < 0.4), aes(PCoA1, PCoA2, color = group, label = group)) +
  geom_point(size = 0.2, alpha = 0.45) +
  geom_text(vjust = -0.7, size = 3, check_overlap = TRUE) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(
    x = xlab,
    y = ylab
  )


###### Connect insertions from the same group to see the groups more clearly ######
make_mst_edges <- function(df) {
  split_df <- split(df, df$group)
  
  edge_list <- lapply(split_df, function(d) {
    if (nrow(d) < 2) return(NULL)
    
    coords <- as.matrix(d[, c("PCoA1", "PCoA2")])
    dm <- as.matrix(dist(coords))
    
    n <- nrow(d)
    selected <- rep(FALSE, n)
    selected[1] <- TRUE
    
    edges <- list()
    
    while (sum(selected) < n) {
      candidates <- which(selected, arr.ind = TRUE)
      nonselected <- which(!selected)
      
      best_i <- NA
      best_j <- NA
      best_dist <- Inf
      
      for (i in which(selected)) {
        for (j in nonselected) {
          if (dm[i, j] < best_dist) {
            best_dist <- dm[i, j]
            best_i <- i
            best_j <- j
          }
        }
      }
      
      edges[[length(edges) + 1]] <- data.frame(
        group = d$group[best_i],
        x = d$PCoA1[best_i],
        y = d$PCoA2[best_i],
        xend = d$PCoA1[best_j],
        yend = d$PCoA2[best_j],
        from = d$name[best_i],
        to = d$name[best_j]
      )
      
      selected[best_j] <- TRUE
    }
    
    do.call(rbind, edges)
  })
  
  do.call(rbind, edge_list)
}

mst.edges <- make_mst_edges(pcoa.df)

ggplot(pcoa.df %>% filter(PCoA1 < 0.4), aes(PCoA1, PCoA2, color = group, label = group)) +
  geom_segment(
    data = mst.edges,
    aes(x = x, y = y, xend = xend, yend = yend, color = group),
    inherit.aes = FALSE,
    alpha = 0.5,
    linewidth = 0.4
  ) +
  geom_point(size = 0.2, alpha = 0.55) +
  geom_text(vjust = -0.7, size = 3, check_overlap = TRUE) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(
    x = xlab,
    y = ylab
  )
