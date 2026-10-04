library(tidyverse)
library(biomaRt)
library(ReactomePA)
library(enrichplot)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# ---- Chunk 1: helper to turn a FASTA of transcripts into an Entrez gene list ----

fasta_to_entrez <- function(fasta_path) {
  fa_lines   <- readLines(fasta_path)
  fa_headers <- fa_lines[str_starts(fa_lines, ">")]        # keep only header lines, not sequence lines

  transcripts <- fa_headers %>%
    str_remove("^>") %>%                                   # drop the leading ">"
    str_remove("\\.\\d+$")                                 # strip version suffix, e.g. ".7"

  id_map <- getBM(
    attributes = c("ensembl_transcript_id", "entrezgene_id"),
    filters    = "ensembl_transcript_id",
    values     = transcripts,
    mart       = ensembl
  )

  unique(na.omit(id_map$entrezgene_id))
}

# ---- Chunk 2: build the stable, unstable, and background gene sets ----

stable_genes     <- fasta_to_entrez("stable_train.fa")
unstable_genes   <- fasta_to_entrez("unstable_train.fa")
background_genes <- fasta_to_entrez("all_human_halflife_data.fa")

# ---- Chunk 3: run ORA separately for the stable and unstable gene sets ----

stable_pathways <- enrichPathway(
  gene          = stable_genes,
  universe      = background_genes,
  organism      = "human",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

unstable_pathways <- enrichPathway(
  gene          = unstable_genes,
  universe      = background_genes,
  organism      = "human",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

# ---- Chunk 4: inspect, plot, and save both results ----

head(as.data.frame(stable_pathways))
head(as.data.frame(unstable_pathways))

dotplot(stable_pathways, showCategory = 15)
ggsave("stable_pathways_dotplot.png", width = 8, height = 8, dpi = 300, bg = "white")

dotplot(unstable_pathways, showCategory = 15)
ggsave("unstable_pathways_dotplot.png", width = 8, height = 8, dpi = 300, bg = "white")

write_csv(as.data.frame(stable_pathways),   "stable_pathways.csv")
write_csv(as.data.frame(unstable_pathways), "unstable_pathways.csv")
