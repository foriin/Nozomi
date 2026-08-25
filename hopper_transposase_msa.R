
library(Biostrings)
library(ggmsa)

#### input files ####

transposases <- "~/Work/projects/DrosoTE/Hopper_A6/PROTEIN/MSA/muscle5_hopper_transposase.fa"

orig_fa <- "~/Work/projects/DrosoTE/Hopper_A6/PROTEIN/transposases.fa"   # original unaligned FASTA, desired order
msa_fa  <- "~/Work/projects/DrosoTE/Hopper_A6/PROTEIN/MSA/muscle5_hopper_transposase.fa"   # MUSCLE output alignment

#### read original and aligned sequences ####

orig <- readAAStringSet(orig_fa)
aln  <- readAAStringSet(msa_fa)

#### check names ####

orig_names <- names(orig)
aln_names  <- names(aln)

missing_in_aln <- setdiff(orig_names, aln_names)
extra_in_aln   <- setdiff(aln_names, orig_names)

if (length(missing_in_aln) > 0) {
  stop("These original sequences are missing from alignment:\n",
       paste(missing_in_aln, collapse = "\n"))
}

if (length(extra_in_aln) > 0) {
  warning("These aligned sequences are not in original FASTA:\n",
          paste(extra_in_aln, collapse = "\n"))
}

#### reorder alignment ####

aln_reordered <- aln[c(6,7,5,8,9,4,1,3,2)]
names(aln_reordered)[5] <- "DBp_Transib-6"
names(aln_reordered)[3] <- "Dsuz_chopper"
names(aln_reordered)[4] <- "Dsuz_interloper"
names(aln_reordered)[6] <- "Bdor_blooper"

#### save reordered alignment ####

writeXStringSet(aln_reordered, "~/Work/projects/DrosoTE/Hopper_A6/PROTEIN/muscle5_transpo_reordered.fa")


#### Supplementary Figure 3A ####

ggmsa(aln_reordered[-7],start = 292, end = 437, char_width = 0.5,
      seq_name = TRUE, disagreement = F, consensus_views = T,
      use_dot = T, border = 'white') + geom_msaBar()

ggmsa("~/Work/projects/DrosoTE/Hopper_A6/PROTEIN/muscle5_transpo_reordered.fa",start = 110, end = 270, char_width = 0.5,
      seq_name = TRUE, disagreement = F, consensus_views = T,
      use_dot = T, border = 'white') + geom_msaBar()

# DNA binding second coil domain
ggmsa("~/Work/projects/DrosoTE/Hopper_A6/PROTEIN/muscle5_transpo_reordered.fa",start = 70, end = 140, char_width = 0.5,
      seq_name = TRUE, disagreement = F, consensus_views = T,
      use_dot = T, border = 'white') + geom_msaBar()

ggmsa("~/Work/projects/DrosoTE/Hopper_A6/PROTEIN/muscle5_transpo_reordered.fa",start = 1, end = 250, char_width = 0.5,
      seq_name = TRUE, disagreement = F, consensus_views = T,
      use_dot = F, border = 'white') + geom_msaBar()



ggmsa(transposases,start = 289, end = 450, char_width = 0.5,
      seq_name = TRUE, ref = "Dmel_hopper", color = "Hydrophobicity",
      disagreement = F, consensus_views = T, border = 'white',
      by_conservation = T) + geom_msaBar()

dbip_hopper <- readDNAStringSet("~/Work/projects/DrosoTE/Hopper_A6/transib6_cand.fa")
dbip_hopper <- subseq(dbip_hopper, 276, 3285)
writeXStringSet(dbip_hopper, "~/Work/projects/DrosoTE/Hopper_A6/transib6_dbip.fa")

hztransib <- readDNAStringSet("~/Work/projects/DrosoTE/Hopper_A6/hztransib.fa")
hztransib <- subseq(hztransib, 1565, 3935)
