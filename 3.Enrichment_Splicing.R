library(clusterProfiler);library(org.Hs.eg.db);library(enrichplot);library(ggplot2);library(stringr);library(data.table)

res <- read.table("C://Users/Gerard/Desktop/AAA/RNAseq/SUPPA/AAA_Results.dpsi")
res <- res[order(res$AAA_Cases.AAA_Controls_p.val),]
res <- na.omit(res)

res[res$AAA_Cases.AAA_Controls_p.val == 0,]$AAA_Cases.AAA_Controls_p.val <- 5e-10
res$BH <- p.adjust(res$AAA_Cases.AAA_Controls_p.val, method = "fdr")

nrow(res[res$AAA_Cases.AAA_Controls_p.val < 0.05,])
nrow(res[res$BH < 0.05,])

rownames(res)

# Save supplementary table:
write.xlsx(res, "C://Users/Gerard/Desktop/AAA/RNAseq/Supplementary_Tables/SUPPA.Results.xlsx", row.names = T)

rownames(res)
rownames(res.sig)
rownames(res.sig.BH)

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
table(event)

event_counts <- table(event)

df <- data.frame(event = names(event_counts), count = 100*(event_counts/15))
rownames(df) <- 1:4
df

new_labels <- c("Alternative 5' Splice Site", "Alternative First Exon",  "Mutually Exclusive Exons", "Skipping Exon")

png("C://Users/Gerard/Desktop/AAA/RNAseq/Plots/Splicing_Piechart.png", width = 1500, height = 1500, res = 200)

ggplot(df, aes(x = "", y = count.Freq, fill = count.event)) +
  geom_bar(stat = "identity", color = "black", size = 0.8) +
  labs(title = "Proportion of splicing events types") +
  coord_polar("y", start = 80) +
  geom_text(aes(label = paste0(round(count.Freq, 2), '%')), position = position_stack(vjust =
                                                                                        0.5)) +
  scale_fill_brewer(palette = "Pastel1",
                    name = "",
                    label = new_labels) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    legend.position = "bottom",
    legend.box.spacing = unit(0, "pt")
  )

dev.off()

# SPP1:
grep("ENSG00000118785", rownames(res.sig.BH), value = T)
res.sig.BH[grep("ENSG00000118785",rownames(res.sig.BH)),]


# FHL1:
grep("ENSG00000022267", rownames(res.sig.BH), value = T)


ann <- fread("C://Users/Gerard/Desktop/AAA/RNAseq/Annotation/gencode.v26.annotation.fixed.gtf.gz")
ann$gene_id <- substr(ann$gene_id, 1 , 15)
head(ann)

res.sig.BH$gene_id <- sub("\\..*", "", rownames(res.sig.BH))

res.sig.BH <- merge(res.sig.BH, ann[,c("gene_id","gene_name")], by = "gene_id")


all.genes <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/AllGenes.xlsx")
all.genes[all.genes$Gene.Name %in% res.sig.BH$gene_name,]


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

cc <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/SigRes/Cases_Controls_No_Smoking.xlsx")
cc$gene_id2 <- substr(cc$Gene.ID, 1, 15)

length(intersect(substr(rownames(res.sig.BH), 1, 15), cc$gene_id2))/length(substr(rownames(res.sig.BH), 1, 15)) # 46.6% Of the significant splicing are DEG.

length(intersect(substr(rownames(res.sig.BH), 1, 15), cc$gene_id2))


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
png("C://Users/Gerard/Desktop/AAA/RNAseq/Plots/Enrichment/GO_Enrichment_Splicing_Nominally_Significant.png", width = 4*800, height = 4*600, res = 300)
dotplot(go_enrich, showCategory = 10) + ggtitle("GO Enrichment")
dev.off()


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



### Sashimi plot:

# Prepare file:

file_names <- list.files(path = "/home/gerard/ggsashimi/Results/SPP1", pattern = "Aligned.sortedByCoord.out.bam$", full.names = TRUE)

data <- data.frame(
  Name = basename(gsub(".Aligned.sortedByCoord.out.bam", "", file_names)),
  Path = file_names
)


data$Group <- as.factor(c(rep("Case", 96), rep("Control",44)))


#write.table(data, file = "/home/gerard/ggsashimi/Results/input_bams_SPP1.tsv", sep = "\t", row.names = FALSE, col.names = F)





