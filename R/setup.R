# Shared setup sourced at the top of every chapter -----------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(limma)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  library(ComplexHeatmap)
  library(circlize)
  library(dendextend)
  library(ggrepel)
  library(patchwork)
  library(cowplot)
  library(scales)
  library(knitr)
  library(kableExtra)
})

# Resolve namespace clashes in favour of the tidyverse verbs used in the book
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
count  <- dplyr::count
slice  <- dplyr::slice

set.seed(111)
options(knitr.kable.NA = "—")

for (d in c("data", "data/raw", "tables", "figures")) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

theme_set(theme_classic(base_size = 12))

# Study design (GEO GSE226609, whole mouse lung) -------------------------------
# Female C57BL/6 mice, wild type or ALOX15 knockout, challenged intranasally
# with saline (PBS) or house dust mite (HDM) extract daily for 14 days.
# The deposited count columns are coded A-D by group; GEO sample titles give
# the mapping used below.

geo_accession <- "GSE226609"
geo_counts_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE226nnn/GSE226609/suppl/",
  "GSE226609_rsem_count.txt.gz"
)
geo_counts_file <- "data/raw/GSE226609_rsem_count.txt.gz"

group_levels <- c("WT_PBS", "WT_HDM", "KO_PBS", "KO_HDM")

group_labels <- c(
  WT_PBS = "WT saline",
  WT_HDM = "WT HDM",
  KO_PBS = "ALOX15 KO saline",
  KO_HDM = "ALOX15 KO HDM"
)

group_codes <- c(A = "WT_PBS", B = "WT_HDM", C = "KO_PBS", D = "KO_HDM")

# Colour-blind-safe palette (Okabe-Ito): blue = saline, vermillion = HDM;
# knockout groups use the lighter tone of each pair.
group_colours <- c(
  WT_PBS = "#0072B2",
  WT_HDM = "#D55E00",
  KO_PBS = "#56B4E9",
  KO_HDM = "#E69F00"
)

treatment_colours <- c(Saline = "#0072B2", HDM = "#D55E00")
de_colours <- c(Down = "#0072B2", NS = "grey75", Up = "#D55E00")

# Thresholds used throughout the book
fdr_cutoff <- 0.05
lfc_cutoff <- 1   # |log2FC| >= 1 (two-fold) for "strong" DE genes

source("R/helpers.R")
source("R/gene_panels.R")
