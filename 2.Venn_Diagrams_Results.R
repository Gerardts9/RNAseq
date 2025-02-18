library(VennDiagram);library(clusterProfiler);library(org.Hs.eg.db);library(enrichplot);library(readxl);library(data.table);library(openxlsx)

ann <- fread("C://Users/Gerard/Desktop/AAA/RNAseq/Annotation/gencode.v26.annotation.fixed.gtf.gz")
ann$gene_id <- substr(ann$gene_id, 1, 15)
head(ann)


# Load Results:
#################################################
cc <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/SigRes/Cases_Controls_Confounders.xlsx")
cc_no_it <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/SigRes/Cases_Controls_Confounders_No_IT.xlsx")
dia <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/SigRes/Diameter.xlsx")

length(unique(c(dia$Gene.Name, sym$Gene.Name)))
length(intersect(cc$Gene.Name, dia$Gene.Name))
intersect(c(dia$Gene.Name,sym$Gene.Name), cc$Gene.Name)


alt <- read.table("C://Users/Gerard/Desktop/AAA/RNAseq/SUPPA/AAA_Sig_Results.txt")
alt$gene_id <- substr(rownames(alt), 1, 15)
alt <- merge(alt, ann, by = "gene_id")

med <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/TunicaSpecificResults.xlsx", sheet = 1)
med$Gene.Symbol <- gsub("-", "", med$Gene.Symbol)

adv <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/TunicaSpecificResults.xlsx", sheet = 2)
adv$Gene.Symbol <- gsub("-", "", adv$Gene.Symbol)

nrow(cc)
nrow(cc_no_it)

intersect(unique(alt$gene_name), cc$Gene.Name)
intersect(unique(alt$gene_name), cc_no_it$Gene.Name)

intersect(unique(dia$Gene.Name), cc$Gene.Name)
intersect(unique(dia$Gene.Name), cc_no_it$Gene.Name)

#################################################

# CaseControl vs Previous Results (Media and Adventitia):

venn.plot <- venn.diagram(
  x = list(A = cc$Gene.Name, B = med$Gene.Symbol, C = adv$Gene.Symbol),
  filename = "C://Users/Gerard/Desktop/AAA/RNAseq/Venn/CasesControls_Media_Adventitia.png",
  category.names = c("New","Media","Adventitia"),
  fill = c("red", "green", "blue"),
  alpha = 0.5)

# CaseControl vs Media

venn.plot <- venn.diagram(
  x = list(A = cc$Gene.Name, B = med$Gene.Symbol),
  filename = "C://Users/Gerard/Desktop/AAA/RNAseq/Venn/CasesControls_Media.png",
  category.names = c("New","Media"),
  fill = c("red", "green"),
  alpha = 0.5)

# CaseControl vs Adventitia:

venn.plot <- venn.diagram(
  x = list(A = cc$Gene.Name, B = adv$Gene.Symbol),
  filename = "C://Users/Gerard/Desktop/AAA/RNAseq/Venn/CasesControls_Adventitia.png",
  category.names = c("New","Adventitia"),
  fill = c("red", "blue"),
  alpha = 0.5)



# Status vs Diameter:
########################################
v <- venn.diagram(
  x = list(A = dia$Gene.Name, B = cc$Gene.Name),
  disable.logging = T,
  category.names = c("Diameter", "Status"),
  fill = c("yellow", "red"),
  alpha = 0.2,
  filename = NULL,
  scaled = FALSE,
  cat.pos = c(-1, 1), 
  cex = 2,
  cat.cex = 2
)


lapply(v, function(i) i$label)

# Intersection:
v[[7]]$label <- paste(intersect(cc$Gene.Name, dia$Gene.Name), collapse="\n")

v[[7]]$label <- substitute(paste(italic("EXTL3")),
                           paste(italic("EXTL3")))


# Plot:
grid.newpage()
grid.draw(v)
########################################


# Status vs Symptoms:
########################################
v <- venn.diagram(
  x = list(A = sym$Gene.Name, B = cc$Gene.Name),
  category.names = c("Symptoms", "Status"),
  fill = c("green", "red"),
  alpha = 0.5,
  filename = NULL,
  scale = FALSE,
  cat.pos = -5
)

lapply(v, function(i) i$label)

# Intersection:
v[[7]]$label <- paste(intersect(cc$Gene.Name, sym$Gene.Name), collapse="\n")  

# Plot:
grid.newpage()
grid.draw(v)
########################################


# Previous Results vs Diameter:
########################################
v <- venn.diagram(
  x = list(A = dia$Gene.Name, B = med$Gene.Symbol, C = adv$Gene.Symbol),
  category.names = c("Diameter", "Media", "Adventitia"),
  fill = c("yellow", "green", "blue"),
  alpha = 0.5,
  filename = NULL,
  scale = FALSE
)

grid.draw(v)

venn.plot <- venn.diagram(
  x = list(A = dia$Gene.Name, B = med$Gene.Symbol, C = adv$Gene.Symbol),
  filename = NULL,
  category.names = c("New","Media", "Adventitia"),
  fill = c("yellow", "green", "blue"),
  alpha = 0.5)


lapply(v, function(i) i$label)

# Intesection:
v[[7]]$label <- paste(intersect(cc$Gene.Name, dia$Gene.Name), collapse="\n")  

# Plot  
grid.newpage()
grid.draw(v)
########################################


# Status vs Diameter vs Symptoms:
########################################
v <- venn.diagram(
  x = list(A = cc$Gene.Name, B = dia$Gene.Name, C = sym$Gene.Name),
  category.names = c("Status", "Diameter", "Symptoms"),
  fill = c("red", "yellow","green"),
  alpha = 0.5,
  filename = NULL,
  cat.pos = c(-10,10,0)
)


lapply(v, function(i) i$label)

v[[10]]$label <- paste(intersect(intersect(cc$Gene.Name, dia$Gene.Name), sym$Gene.Name), collapse="\n")  
v[[9]]$label <- paste(intersect(cc$Gene.Name, dia$Gene.Name)[!intersect(cc$Gene.Name, dia$Gene.Name) == "LAMA2"], 
                      collapse="\n")
v[[11]]$label <- paste(intersect(cc$Gene.Name, sym$Gene.Name)[!intersect(cc$Gene.Name, sym$Gene.Name) == "LAMA2"],
                      collapse="\n")

grid.newpage()
grid.draw(v)
########################################


# Splicing vs Status:
########################################
v.splicing <- venn.diagram(sub.cex = 10,
  x = list(A = alt$gene_name, B = cc$Gene.Name),
  disable.logging = T,
  category.names = c("Splicing", "Status"),
  fill = c("cornflowerblue", "red"),
  alpha = 0.2,
  filename = NULL,
  cat.pos = c(-1, 1), 
  scaled = F,
  cex = 2,
  cat.cex = 2
)

lapply(v.splicing, function(i) i$label)

intersect(alt$gene_name, cc$Gene.Name)


v.splicing[[6]]$label <- paste(setdiff(alt$gene_name, cc$Gene.Name), collapse="\n")  
v.splicing[[7]]$label <- paste(intersect(alt$gene_name, cc$Gene.Name), collapse="\n")  

grid.newpage()
grid.draw(v.splicing)
########################################














