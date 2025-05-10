# ================================================================
#  quick_enrich.R – one-function wrapper around clusterProfiler
#  Joseph Boktor • May-2025 (GPT assisted scritp)
# ================================================================
suppressPackageStartupMessages({
    library(DESeq2) # differential expression infrastructure
    library(clusterProfiler) # ORA & GSEA framework
    library(enrichplot) # ggplot2 extensions for enrichment objects
    library(org.Mm.eg.db) # ↔ mouse Entrez annotation
    library(AnnotationDbi)
    library(tidyverse)
})

# ---------- main -------------------------------------------------
run_GO_enrich <- function(res,
                          alpha = 0.05, # FDR threshold
                          lfc = 1, # abs(log2FC) threshold
                          ont = "BP", # GO ontology
                          top_n = 20, # terms to show
                          plot_file = "GO_dotplot.pdf") {
    
    ## filter significant DEGs -------------------------------------
    sig <- as.data.frame(res) |>
        tibble::rownames_to_column("ensembl") |>
        dplyr::filter(!is.na(padj), padj < alpha, abs(log2FoldChange) > lfc)

    if (nrow(sig) < 10) {
        stop("Fewer than 10 significant genes – enrichment not meaningful.")
    }

    ## convert Ensembl → Entrez (strip version suffix) -------------
    sig$ensembl <- sub("\\..*$", "", sig$ensembl) # ENSMUSG… → ENSMUSG…
    genes <- clusterProfiler::bitr(sig$ensembl,
        fromType = "ENSEMBL",
        toType   = "ENTREZID",
        OrgDb    = org.Mm.eg.db,
        drop     = TRUE
    )$ENTREZID

    ## perform enrichment ------------------------------------------
    ego <- clusterProfiler::enrichGO(
        gene = genes,
        OrgDb = org.Mm.eg.db,
        keyType = "ENTREZID",
        ont = ont,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
    ) %>%
        clusterProfiler::simplify()

    ## visualise & save --------------------------------------------
    p <- enrichplot::dotplot(ego,
        showCategory = top_n,
        font.size    = 10) +
        ggtitle("GO enrichment (clusterProfiler)") +
        theme_bw()
    ggsave(plot_file, p, width = 6.5, height = 8)

    return(invisible(list(results = ego, plot = p)))
}
