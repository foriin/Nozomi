suppressPackageStartupMessages({
  library(Biostrings)
  library(tidyverse)
  library(ggmsa)
})

#### input ####

aln_file <- "~/Work/projects/DrosoTE/Hopper_A6/PROTEIN/MSA/muscle5_hopper_transposase.fa"
ref_name <- "Dmel_Nozomi"

#### read aligned protein FASTA ####

aln <- readAAStringSet(aln_file)

if (!ref_name %in% names(aln)) {
  stop(
    "Reference sequence not found: ", ref_name, "\n",
    "Available names:\n",
    paste(names(aln), collapse = "\n")
  )
}

aln <- aln[-1]
names(aln)[1] <- "HzTransib (Hzea)"
names(aln)[2] <- "Transib1 (Dmel)"
names(aln)[3] <- "Oriex (Bdor)"
names(aln)[4] <- "Darjeeling (Dsuz)"
names(aln)[5] <- "Nozomi (Dmel)"
names(aln)[6] <- "Nozomi (Dmau)"
names(aln)[7] <- "Hayabusa (Dsuz)"
names(aln)[8] <- "Transib-6 (Dbip)"

id_tbl_rnase$sequence[1] <- "Nozomi (Dmel)"
id_tbl_rnase$sequence[2] <- "Nozomi (Dmau)"
id_tbl_rnase$sequence[3] <- "Darjeeling (Dsuz)"
id_tbl_rnase$sequence[4] <- "Hayabusa (Dsuz)"
id_tbl_rnase$sequence[5] <- "Transib-6 (Dbip)"
id_tbl_rnase$sequence[6] <- "Oriex (Bdor)"
id_tbl_rnase$sequence[7] <- "Transib1 (Dmel)"
id_tbl_rnase$sequence[8] <- "HzTransib (Hzea)"

#### function ####

aa_identity_to_ref <- function(aln, ref_name, count_gap_mismatch = TRUE) {
  ref <- strsplit(as.character(aln[[ref_name]]), "")[[1]]
  
  purrr::map_dfr(names(aln), function(nm) {
    query <- strsplit(as.character(aln[[nm]]), "")[[1]]
    
    if (length(query) != length(ref)) {
      stop("Sequences are not aligned or have different lengths: ", nm)
    }
    
    both_gap <- ref == "-" & query == "-"
    
    if (count_gap_mismatch) {
      keep <- !both_gap
    } else {
      keep <- ref != "-" & query != "-"
    }
    
    n_compared <- sum(keep)
    n_identical <- sum(ref[keep] == query[keep])
    
    tibble(
      sequence = nm,
      ref = ref_name,
      n_compared = n_compared,
      n_identical = n_identical,
      aa_identity = n_identical / n_compared,
      aa_identity_pct = 100 * aa_identity
    )
  }) |>
    arrange(desc(aa_identity_pct))
}

#### calculate identity ####
#### count_gap_mismatch = TRUE: gaps in one sequence count as mismatches
#### count_gap_mismatch = FALSE: compare only aligned AA-vs-AA positions

id_tbl <- aa_identity_to_ref(
  aln = aln,
  ref_name = ref_name,
  count_gap_mismatch = FALSE
)

#### Sigh... Rename everything ####
id_tbl$ref <- "Nozomi (Dmel)"
id_tbl$sequence[1] <- "Nozomi (Dmel)"
id_tbl$sequence[2] <- "Nozomi (Dmau)"
id_tbl$sequence[3] <- "Darjeeling (Dsuz)"
id_tbl$sequence[4] <- "Hayabusa (Dsuz)"
id_tbl$sequence[5] <- "Transib-6 (Dbip)"
id_tbl$sequence[6] <- "Oriex (Bdor)"
id_tbl$sequence[7] <- "Transib1 (Dmel)"
id_tbl$sequence[8] <- "HzTransib (Hzea)"


#### save table ####

write_tsv(id_tbl, "aa_identity_to_reference.tsv")

#### plot ####

p_transposase <- ggplot(id_tbl %>% filter(sequence != "Nozomi (Dmel)"), aes(x = reorder(sequence, aa_identity_pct), y = aa_identity_pct)) +
  geom_col() +
  geom_hline(
    yintercept = c(50, 75),
    color = "red",
    linetype = "dashed"
  ) +
  coord_flip() +
  theme_classic() +
  labs(
    x = NULL,
    y = "Amino acid identity to reference (%)",
    title = paste("AA identity to Nozomi (Dmel)")
  )

ggsave(
  filename = "outputs/plots/aa_identity_to_Nozomi_transp.pdf",
  plot = p_transposase,
  width = 8,
  height = max(4, 0.25 * nrow(id_tbl))
)

id_tbl


#### Catalytic RNAse-H like Domain ####

rnaseh <- subseq(aln , 292, 437)

id_tbl_rnase <- aa_identity_to_ref(
  aln = rnaseh,
  ref_name = ref_name,
  count_gap_mismatch = FALSE
)

id_tbl_rnase$ref <- "Nozomi (Dmel)"
id_tbl_rnase$sequence[1] <- "Nozomi (Dmel)"
id_tbl_rnase$sequence[2] <- "Nozomi (Dmau)"
id_tbl_rnase$sequence[3] <- "Darjeeling (Dsuz)"
id_tbl_rnase$sequence[4] <- "Hayabusa (Dsuz)"
id_tbl_rnase$sequence[5] <- "Transib-6 (Dbip)"
id_tbl_rnase$sequence[6] <- "Oriex (Bdor)"
id_tbl_rnase$sequence[7] <- "Transib1 (Dmel)"
id_tbl_rnase$sequence[8] <- "HzTransib (Hzea)"


#### plot ####

p_rnase <- ggplot(
  id_tbl_rnase %>% filter(sequence != "Nozomi (Dmel)"),
  aes(x = reorder(sequence, aa_identity_pct), y = aa_identity_pct)
) +
  geom_col() +
  geom_hline(
    yintercept = c(50, 75),
    color = "red",
    linetype = "dashed"
  ) +
  coord_flip() +
  theme_classic() +
  labs(
    x = NULL,
    y = "Amino acid identity to reference (%)",
    title = paste("Catalytic domain AA identity to Nozomi (Dmel)")
  )

ggsave(
  filename = "outputs/plots/aa_identity_to_Nozomi_cat.pdf",
  plot = p_rnase,
  width = 8,
  height = max(4, 0.25 * nrow(id_tbl))
)

#### GGMSA of catalytic domain ####

aln_reord <- aln[id_tbl_rnase$sequence]

ggmsa(aln_reord,start = 292, end = 437, char_width = 0.5,
      seq_name = TRUE, disagreement = F, consensus_views = T,
      use_dot = T, border = 'white') + geom_msaBar()
