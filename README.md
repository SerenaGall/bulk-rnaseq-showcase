# House dust mite allergic airway inflammation: bulk RNA-seq workflow

A reproducible Quarto book that walks through a complete bulk RNA-seq analysis of whole mouse lung with edgeR, using public data from [GEO GSE226609](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226609): wild-type and ALOX15 knockout C57BL/6 mice challenged with house dust mite (HDM) or saline, five mice per group.

**Read the book:** https://serenagall.github.io/bulk-rnaseq-showcase/

## Contents

QC and tissue-contamination scoring · Filtering and TMM normalisation · MDS / PCA · Factorial edgeR quasi-likelihood model with contamination covariates · Volcano, MA and heatmaps · KEGG ORA and GSEA · Curated type 2 allergy panels · ALOX15 knockout: confounder flagging, passenger-gene mapping and interaction analysis

## Reproduce

Requires R (≥ 4.4) with edgeR, clusterProfiler, org.Mm.eg.db, ComplexHeatmap, kableExtra, dendextend and the tidyverse, plus Quarto.

```bash
quarto render
```

The first chapter downloads the GEO count matrix (~1 MB) to `data/raw/`. Computed results are cached in `_freeze/`; the GitHub Actions workflow publishes the site from that cache without re-running R.

## Data note

The deposited matrix contains RSEM expected counts keyed by UCSC gene symbols (some are older names, e.g. *Clca3* = gob-5/*Clca1*). Column codes A–D map to WT saline, WT HDM, KO saline and KO HDM according to the GEO sample records.
