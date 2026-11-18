# Curated gene panels for type 2 allergic airway inflammation (mouse symbols) --
# Assembled from the allergic-airway literature for this showcase; each panel
# is a short list of well-characterised marker or effector genes. Symbols follow
# the UCSC gene set of the deposited counts (Clca3 = gob-5, now Clca1; Igj = Jchain).

allergy_panels <- list(
  "Th2 cytokines & alarmins" = c(
    "Il4", "Il5", "Il13", "Il9", "Il10", "Il33", "Tslp", "Il25", "Areg",
    "Il1rl1", "Il4ra", "Il13ra2"
  ),
  "Eosinophil chemokines & recruitment" = c(
    "Ccl11", "Ccl24", "Ccl8", "Ccl17", "Ccl22", "Ccl2", "Ccl7", "Ccl12",
    "Ccr3", "Ccr4", "Ccr8"
  ),
  "Eosinophils" = c(
    "Prg2", "Prg3", "Epx", "Ear1", "Ear2", "Ear6", "Siglecf", "Il5ra", "Alox15"
  ),
  "Mast cells & IgE" = c(
    "Cpa3", "Mcpt1", "Mcpt2", "Mcpt4", "Mcpt8", "Cma1", "Tpsb2", "Fcer1a",
    "Ms4a2", "Kit", "Hdc"
  ),
  "Alternatively activated macrophages" = c(
    "Chil3", "Chil4", "Retnla", "Retnlg", "Arg1", "Mrc1", "Mgl2", "Ear10",
    "Ch25h", "Fn1"
  ),
  "Mucus & goblet cells" = c(
    "Muc5ac", "Muc5b", "Clca3", "Agr2", "Foxa3", "Spdef", "Tff2", "Itln1",
    "Slc26a4", "Scgb3a1"
  ),
  "Type 2 transcription factors" = c(
    "Gata3", "Stat6", "Irf4", "Batf", "Pparg", "Foxp3", "Tbx21", "Rorc", "Bcl6"
  ),
  "Lipid mediators" = c(
    "Alox15", "Alox5", "Alox5ap", "Ltc4s", "Cysltr1", "Ptgs2", "Hpgds",
    "Ptgdr2", "Pla2g4a"
  ),
  # Immunoglobulin heavy-chain loci are not in the UCSC gene set used by the
  # depositors, so plasma-cell and B-cell genes stand in for the antibody response
  "Antibody & B-cell response" = c(
    "Igj", "Mzb1", "Xbp1", "Prdm1", "Tnfrsf17", "Pou2af1", "Cd79a", "Aicda"
  )
)

# Long form: one row per gene; a gene listed in several panels keeps its first
allergy_panel_tbl <- tibble::enframe(allergy_panels, name = "block", value = "SYMBOL") %>%
  tidyr::unnest(SYMBOL) %>%
  dplyr::mutate(block = factor(block, levels = names(allergy_panels)))

# Tissue-of-origin markers used to detect dissection carry-over ----------------
# `core` genes define each per-sample score; `flag` genes (core plus other
# highly tissue-specific transcripts) are marked in DE tables so they are not
# read as lung biology.
contamination_markers <- list(
  "Skeletal muscle" = list(
    core = c("Acta1", "Myh1", "Myh4", "Ckm", "Tnnt3", "Mylpf", "Pvalb", "Atp2a1", "Actn3", "Myl1"),
    flag = c("Acta1", "Myh1", "Myh2", "Myh4", "Myh8", "Ckm", "Tnnt1", "Tnnt3", "Tnni2", "Tnnc2",
             "Mylpf", "Myl1", "Pvalb", "Atp2a1", "Actn2", "Actn3", "Mb", "Casq1", "Ryr1", "Myot",
             "Neb", "Ttn", "Mybpc2", "Cox6a2", "Ampd1", "Kbtbd10", "Smpx", "Myoz1", "Eno3")
  ),
  "Brown fat (mediastinum)" = list(
    core = c("Ucp1", "Cidea", "Adipoq", "Cfd", "Plin1", "Car3", "Pck1", "Thrsp", "Cox8b", "Retn"),
    flag = c("Ucp1", "Cidea", "Cidec", "Adipoq", "Cfd", "Plin1", "Plin4", "Car3", "Pck1", "Thrsp",
             "Cox8b", "Cox7a1", "Retn", "Lep", "Elovl3", "Dio2", "Ppargc1a")
  ),
  "Liver" = list(
    core = c("Alb", "Apoa1", "Apoc3", "Fga", "Fgb", "Ahsg", "Cps1", "Aldob", "Serpina1a", "Ttr"),
    flag = c("Alb", "Apoa1", "Apoa2", "Apoa5", "Apoc1", "Apoc2", "Apoc3", "Apoc4", "Apoh", "Apob",
             "Ahsg", "Fga", "Fgb", "Fgg", "Cps1", "Aldob", "Serpina1a", "Serpina1b", "Serpina1c",
             "Serpina1d", "Serpina1e", "Serpina3k", "Ttr", "Tf", "Hpx", "Kng1", "Gc", "Pzp", "Fabp1",
             "Itih1", "Itih2", "Itih3", "Itih4", "F2", "Plg", "Hrg", "Cyp2e1", "Cyp3a11", "Ass1",
             "Mup3", "Mup20", "Rbp4", "Ambp", "Fetub", "Hamp", "Serpinc1", "C8a", "C8b", "C9",
             "Mat1a", "Tat", "Hpd", "Hgd", "Bhmt", "Gnmt", "Apom", "Ces1c", "Slc10a1",
             "Cyp2c29", "Cyp7a1")
  ),
  "Red blood cells" = list(
    core = c("Hba-a1", "Hbb-bs", "Alas2", "Slc4a1", "Gypa", "Rhd", "Kel", "Epb4.2"),
    flag = c("Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt", "Hbb-b1", "Hbb-b2", "Alas2", "Slc4a1",
             "Gypa", "Rhd", "Rhag", "Kel", "Epb4.2", "Spta1", "Sptb", "Ank1", "Trim10", "Trim58",
             "Bpgm", "Snca", "Ermap", "Tmcc2", "Fech", "Urod", "Add2", "Tspo2", "Slc25a37",
             "Cldn13", "Klf1", "Gata1", "Rhced", "Hemgn")
  )
)

contamination_flag_tbl <- purrr::imap_dfr(
  contamination_markers, ~ tibble::tibble(SYMBOL = .x$flag, tissue_flag = .y)
) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

contamination_scores <- function(logcpm) {
  sapply(contamination_markers, function(m) {
    colMeans(logcpm[intersect(m$core, rownames(logcpm)), , drop = FALSE])
  })
}

# Keyword rules for grouping enriched KEGG pathways into broader themes.
# The first matching rule wins; anything unmatched becomes "Other".
kegg_theme_rules <- tibble::tribble(
  ~pattern,                                                              ~theme,
  "cytokine|chemokine|JAK-STAT|IL-17|TNF|NF-kappa|Toll|NOD|RIG|C-type",   "Cytokine & innate signalling",
  "antigen|Th1|Th2|Th17|T cell|B cell|Fc|IgA|hematopoietic|leukocyte",    "Adaptive immunity & antigen presentation",
  "asthma|allerg|inflammatory bowel|rheumatoid|lupus|graft|autoimmune",   "Immune-mediated disease",
  "infection|virus|viral|bacterial|tuberculosis|leishmania|malaria|toxoplasm|pertussis|influenza|measles|herpes|hepatitis|COVID|Chagas|amoebiasis|shigellosis|salmonella|legionell|yersinia|staphylococcus|leukemia virus|E. coli|Epstein", "Host defence (infection pathways)",
  "cell cycle|DNA replication|mismatch|p53|meiosis|base excision|recombination|Fanconi",        "Cell cycle & DNA replication",
  "ECM|focal adhesion|adhesion|cytoskeleton|tight junction|migration",    "ECM, adhesion & migration",
  "metabolism|biosynthesis|degradation|oxidative|PPAR|fatty|arachidonic|linoleic|glutathione|steroid", "Metabolism & lipid mediators",
  "complement|coagulation|phagosome|lysosome|apoptosis|necroptosis|extracellular trap",      "Complement, phagocytosis & cell death"
)

assign_kegg_theme <- function(description) {
  out <- rep("Other", length(description))
  for (i in rev(seq_len(nrow(kegg_theme_rules)))) {
    hit <- grepl(kegg_theme_rules$pattern[i], description, ignore.case = TRUE)
    out[hit] <- kegg_theme_rules$theme[i]
  }
  out
}
