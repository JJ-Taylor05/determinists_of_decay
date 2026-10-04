library(tidyverse)
library(biomaRt)   
library(ReactomePA)
library(enrichplot)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# ---- Chunk 1: load hits and flag which ones fall in the 5' CDS peak ----

peak_fraction <- 0.10   # first 10% of CDS length - adjust this to widen/narrow the peak window

hits <- read_tsv("plot_data/stable_motif_mapping_hits.tsv", na = "NA")

peak_hits <- hits %>%
  filter(region == "cds") %>%                                        # only hits inside the CDS
  mutate(
    hit_centre   = (start + end) / 2,                                # midpoint of each motif hit
    cds_fraction = (hit_centre - orf_start) / (orf_end - orf_start)  # position within CDS, 0 = start, 1 = end
  ) %>%
  filter(cds_fraction >= 0, cds_fraction <= peak_fraction)            # keep only hits in the 5' peak window

# ---- Chunk 2: get a clean, unique transcript list ----

peak_transcripts <- peak_hits %>%
  distinct(seq_id) %>%
  mutate(ensembl_transcript = str_remove(seq_id, "\\.\\d+$"))   # strip version suffix, e.g. ".8"

# ---- Chunk 3: convert transcript IDs to Entrez gene IDs ----

gene_map <- getBM(
  attributes = c("ensembl_transcript_id", "entrezgene_id", "hgnc_symbol"),
  filters = "ensembl_transcript_id",
  values = peak_transcripts$ensembl_transcript,
  mart = ensembl
)

gene_list <- unique(na.omit(gene_map$entrezgene_id))   # this is the gene list ReactomePA needs

# ---- Chunk 4: build the background gene set from the FASTA headers ----

fa_lines   <- readLines("plot_data/all_human_halflife_data.fa")
fa_headers <- fa_lines[str_starts(fa_lines, ">")]              # keep only header lines, not sequence lines

background_transcripts <- fa_headers %>%
  str_remove("^>") %>%                                          # drop the leading ">"
  str_remove("\\.\\d+$")                                        # strip version suffix, e.g. ".7"

background_map <- getBM(
  attributes = c("ensembl_transcript_id", "entrezgene_id"),
  filters = "ensembl_transcript_id",
  values = background_transcripts,
  mart = ensembl
)

background_genes <- unique(background_map$entrezgene_id)

# ---- Chunk 5: run the enrichment ----

peak_pathways <- enrichPathway(
  gene          = gene_list,
  universe      = background_genes,
  organism      = "human",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE          # converts Entrez IDs back to gene symbols in the output table
)

head(as.data.frame(peak_pathways))

# ---- Chunk 6: visualize and save ----

dotplot(peak_pathways, showCategory = 15)

write_csv(as.data.frame(peak_pathways), "5prime_cds_peak_pathways.csv")
ggsave("5prime_cds_peak_dotplot.png", width = 8, height = 8, dpi = 300, bg = "white")
