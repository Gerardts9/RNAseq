library(clusterProfiler);library(org.Hs.eg.db);library(enrichplot);library(ggplot2);library(stringr);library(data.table)

res <- read.table("C://Users/Gerard/Desktop/AAA/RNAseq/SUPPA/AAA_Results.dpsi")
res <- res[order(res$AAA_Cases.AAA_Controls_p.val),]
res <- na.omit(res)

res[grep("ENSG00000146648", rownames(res)),]

res$BH <- p.adjust(res$AAA_Cases.AAA_Controls_p.val, method = "fdr")

nrow(res[res$AAA_Cases.AAA_Controls_p.val < 0.05,])
nrow(res[res$BH < 0.05,])

# Nominally significant genes:
res.sig <- res[res$AAA_Cases.AAA_Controls_p.val < 0.05,]

# Significant BH genes:
res.sig.BH <- res[res$BH < 0.05,]

#write.table(res.sig.BH, "C://Users/Gerard/Desktop/AAA/RNAseq/SUPPA/AAA_Sig_Results.txt") # 15

length(table(table(substr(rownames(res),1,15))))

sum(table(table(substr(rownames(res),1,15)))[1:5])
sum(table(table(substr(rownames(res),1,15)))[5:51])

table(table(substr(rownames(res),1,15)))[table(table(substr(rownames(res),1,15))) > 5]

hist(res$AAA_Cases.AAA_Controls_p.val)

event <- sub(".*;(..).*","\\1", rownames(res.sig.BH))

ggplot(data.frame(event), aes(x=event, fill = event)) +
  geom_bar() + scale_y_continuous(breaks = seq(0,max(table(event)), by = 2))

#ggsave("C://Users/Gerard/Desktop/AAA/RNAseq/Plots/Splicing_Barplot.png")


# All genes:
gene_list <- unique(substr(rownames(res), 1, 15))
length(unique(gene_list)) # 9371 unique genes

# Nominally significant genes:
gene_list <- unique(substr(rownames(res.sig), 1, 15))
length(unique(gene_list)) # 261 unique genes

# Significant BH genes:
gene_list <- unique(substr(rownames(res.sig.BH), 1, 15))
length(unique(gene_list)) # 11 unique genes

symbol_ids <- bitr(gene_list, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
symbol_ids

cc <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/SigRes/Cases_Controls.xlsx")
cc$gene_id2 <- substr(cc$Gene.ID, 1, 15)

length(intersect(substr(rownames(res.sig.BH), 1, 15), cc$gene_id2))/length(substr(rownames(res.sig.BH), 1, 15)) # 46.6% Of the significant splicing are DEG.

intersect(substr(rownames(res.sig.BH), 1, 15), cc$gene_id2)


# Perform functional enrichment analysis:

entrez_ids <- bitr(gene_list, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# GO enrichment analysis
go_enrich <- enrichGO(gene          = entrez_ids$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "ALL", # Choose from "BP" (Biological Process), "CC" (Cellular Component), or "MF" (Molecular Function)
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2)

# GO enrichment dot plot:
#png("C://Users/Gerard/Desktop/GO_Enrichment_Splicing_Nominally_Significant.png", width = 4*800, height = 4*600, res = 300)
dotplot(go_enrich, showCategory = 10) + ggtitle("GO Enrichment")
#dev.off()


gene_id_list <- lapply(go_enrich@result$geneID, function(x) strsplit(x, "/")[[1]])

gene_symbols_df <- bitr(unlist(gene_id_list), fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)

gene_symbol_list <- lapply(gene_id_list, function(x) gene_symbols_df$SYMBOL[match(x, gene_symbols_df$ENTREZID)])
names(gene_symbol_list) <- names(gene_id_list)


# KEGG enrichment analysis: No results.
#kegg_enrich <- enrichKEGG(gene         = entrez_ids$ENTREZID,
#                          organism     = "hsa",
#                          keyType      = "kegg",
#                          pAdjustMethod = "BH",
#                          pvalueCutoff = 0.05,
#                          qvalueCutoff = 0.2)








