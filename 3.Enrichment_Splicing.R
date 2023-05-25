library(clusterProfiler);library(org.Hs.eg.db);library(enrichplot);library(ggplot2);library(stringr)

res <- read.table("C://Users/Gerard/Desktop/AAA/RNAseq/SUPPA/AAA_Sig_Results.txt")

event <- sub(".*;(..).*","\\1", rownames(res))

ggplot(data.frame(event), aes(x=event, fill = event)) +
  geom_bar()

#ggsave("C://Users/Gerard/Desktop/AAA/RNAseq/Plots/Splicing_Barplot.png")


gene_list <- unique(substr(rownames(res), 1, 15))

entrez_ids <- bitr(gene_list, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# GO enrichment analysis
go_enrich <- enrichGO(gene          = entrez_ids$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "ALL", # Choose from "BP" (Biological Process), "CC" (Cellular Component), or "MF" (Molecular Function)
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2)


gene_id_list <- lapply(go_enrich@result$geneID, function(x) strsplit(x, "/")[[1]])

gene_symbols_df <- bitr(unlist(gene_id_list), fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)

gene_symbol_list <- lapply(gene_id_list, function(x) gene_symbols_df$SYMBOL[match(x, gene_symbols_df$ENTREZID)])
names(gene_symbol_list) <- names(gene_id_list)


med <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/TunicaSpecificResults.xlsx", sheet = 1)
adv <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/TunicaSpecificResults.xlsx", sheet = 2)



intersect(cc$V1, unique(unlist(gene_symbol_list)))


med$Gene.Symbol <- gsub("-", "", med$Gene.Symbol)
adv$Gene.Symbol <- gsub("-", "", adv$Gene.Symbol)

med[med$Gene.Symbol %in% "FHL1",]
adv[adv$Gene.Symbol %in% "FHL1",]

med[med$Gene.Symbol %in% "VCL",]
adv[adv$Gene.Symbol %in% "VCL",]

med[med$Gene.Symbol %in% "CAV2",]
adv[adv$Gene.Symbol %in% "CAV2",]



# GO enrichment dot plot:
png("C://Users/Gerard/Desktop/GO_Enrichment_Splicing.png", width = 4*800, height = 4*600, res = 300)
dotplot(go_enrich, showCategory = 10) + ggtitle("GO Enrichment")
dev.off()



# KEGG enrichment analysis: No results.
#kegg_enrich <- enrichKEGG(gene         = entrez_ids$ENTREZID,
#                          organism     = "hsa",
#                          keyType      = "kegg",
#                          pAdjustMethod = "BH",
#                          pvalueCutoff = 0.05,
#                          qvalueCutoff = 0.2)








