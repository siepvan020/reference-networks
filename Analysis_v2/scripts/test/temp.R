



ranked_genes <- residuals_list[["blood_platelet_vs_lung_2pneumo"]]





# > names(gsea_results)
#  [1] "BP_cd4_blood_vs_lung"              "MF_cd4_blood_vs_lung"
#  [3] "BP_cd8_blood_vs_lung"              "MF_cd8_blood_vs_lung"
#  [5] "BP_monocyte_blood_vs_lung"         "MF_monocyte_blood_vs_lung"
#  [7] "BP_lung_cd4_vs_cd8"                "MF_lung_cd4_vs_cd8"
#  [9] "BP_lung_cd4_vs_monocyte"           "MF_lung_cd4_vs_monocyte"
# [11] "BP_lung_cd8_vs_monocyte"           "MF_lung_cd8_vs_monocyte"
# [13] "BP_blood_platelet_vs_lung_2pneumo" "MF_blood_platelet_vs_lung_2pneumo"


sig_pathways <- gsea_results[[8]] %>%
    filter(padj < 0.05)

topPathwaysUp <- sig_pathways[NES > 0] %>%
    arrange(desc(NES)) %>%
    head(10)
topPathwaysUp$pathway <- factor(topPathwaysUp$pathway, levels = topPathwaysUp$pathway)
topPathwaysUp$pathway <- gsub("_", " ", topPathwaysUp$pathway)
topPathwaysUp$pathway <- str_wrap(topPathwaysUp$pathway, width = 60, whitespace_only = FALSE)
topPathwaysUp$pathway <- gsub(" ", "_", topPathwaysUp$pathway)


topPathwaysDown <- sig_pathways[NES < 0] %>%
    arrange(NES) %>%
    head(10)
topPathwaysDown$pathway <- factor(topPathwaysDown$pathway, levels = topPathwaysDown$pathway)
topPathwaysDown$pathway <- gsub("_", " ", topPathwaysDown$pathway)
topPathwaysDown$pathway <- str_wrap(topPathwaysDown$pathway, width = 60, whitespace_only = FALSE)
topPathwaysDown$pathway <- gsub(" ", "_", topPathwaysDown$pathway)



plotEnrichment(
    pathway = msigdb_BP[["GOBP_GAS_TRANSPORT"]], # Get gene set for first pathway
    stats = ranked_genes
) + labs(title = paste("Enrichment Plot:", "GOBP_GAS_TRANSPORT"))






# temporarily store this here, didnt really seem to work properly
get_top_pathways <- function(gsea_result) {
    sig_pathways <- gsea_result %>%
        filter(padj < 0.05)

    topPathwaysUp <- sig_pathways[NES > 0] %>%
        arrange(desc(NES)) %>%
        head(10)
    # topPathwaysUp$pathway <- str_trunc(as.character(topPathwaysUp$pathway), width = 60)
    topPathwaysUp$pathway <- factor(topPathwaysUp$pathway, levels = topPathwaysUp$pathway)
    topPathwaysUp$pathway <- gsub("_", " ", topPathwaysUp$pathway)
    topPathwaysUp$pathway <- str_wrap(topPathwaysUp$pathway, width = 50, whitespace_only = FALSE)
    topPathwaysUp$pathway <- gsub(" ", "_", topPathwaysUp$pathway)


    topPathwaysDown <- sig_pathways[NES < 0] %>%
        arrange(NES) %>%
        head(10)
    # topPathwaysDown$pathway <- str_trunc(as.character(topPathwaysDown$pathway), width = 60)
    topPathwaysDown$pathway <- factor(topPathwaysDown$pathway, levels = topPathwaysDown$pathway)
    topPathwaysDown$pathway <- gsub("_", " ", topPathwaysDown$pathway)
    topPathwaysDown$pathway <- str_wrap(topPathwaysDown$pathway, width = 50, whitespace_only = FALSE)
    topPathwaysDown$pathway <- gsub(" ", "_", topPathwaysDown$pathway)

    # topPathways <- bind_rows(topPathwaysUp, topPathwaysDown) %>%
    #     arrange(NES)
    # topPathways$pathway <- factor(topPathways$pathway, levels = topPathways$pathway)

    return(list(up = topPathwaysUp, down = topPathwaysDown))
}
