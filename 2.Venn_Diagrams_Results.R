library(VennDiagram)

ann <- fread("C://Users/Gerard/Desktop/AAA/RNAseq/Annotation/gencode.v26.annotation.fixed.gtf.gz")
ann$gene_id <- substr(ann$gene_id, 1, 15)
head(ann)


# Load Results:
#################################################
cc <- fread("C://Users/Gerard/Desktop/AAA/RNAseq/SigRes/Cases_Controls.txt")
dia <- fread("C://Users/Gerard/Desktop/AAA/RNAseq/SigRes/Diameter.txt")
sym <- fread("C://Users/Gerard/Desktop/AAA/RNAseq/SigRes/Symptoms.txt")

alt <- read.table("C://Users/Gerard/Desktop/AAA/RNAseq/SUPPA/AAA_Sig_Results.txt")
alt$gene_id <- substr(rownames(alt), 1, 15)
alt <- merge(alt, ann, by = "gene_id")


med <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/TunicaSpecificResults.xlsx", sheet = 1)
med$Gene.Symbol <- gsub("-", "", med$Gene.Symbol)

adv <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/TunicaSpecificResults.xlsx", sheet = 2)
adv$Gene.Symbol <- gsub("-", "", adv$Gene.Symbol)
#################################################


# CaseControl vs Previous Results (Media and Adventitia):

venn.plot <- venn.diagram(
  x = list(A = cc$V1, B = med$Gene.Symbol, C = adv$Gene.Symbol),
  filename = "C://Users/Gerard/Desktop/AAA/RNAseq/Venn/CasesControls_Media_Adventitia.png",
  category.names = c("New","Media","Adventitia"),
  fill = c("red", "green", "blue"),
  alpha = 0.5)


# CaseControl vs Media

venn.plot <- venn.diagram(
  x = list(A = cc$V1, B = med$Gene.Symbol),
  filename = "C://Users/Gerard/Desktop/AAA/RNAseq/Venn/CasesControls_Media.png",
  category.names = c("New","Media"),
  fill = c("red", "green"),
  alpha = 0.5)

# CaseControl vs Adventitia:

venn.plot <- venn.diagram(
  x = list(A = cc$V1, B = adv$Gene.Symbol),
  filename = "C://Users/Gerard/Desktop/AAA/RNAseq/Venn/CasesControls_Adventitia.png",
  category.names = c("New","Adventitia"),
  fill = c("red", "blue"),
  alpha = 0.5)



# CasesControls vs Diameter:
########################################
v <- venn.diagram(
  x = list(A = dia$V1, B = cc$V1),
  category.names = c("Diameter", "CaseControl"),
  fill = c("yellow", "red"),
  alpha = 0.5,
  filename = NULL,
  scale = FALSE,
  cat.pos = -5
)

lapply(v, function(i) i$label)

# Intersection:
v[[7]]$label <- paste(intersect(cc$V1, dia$V1), collapse="\n")  

# Plot:
grid.newpage()
grid.draw(v)
########################################

# CasesControls vs Symptoms:
########################################
v <- venn.diagram(
  x = list(A = sym$V1, B = cc$V1),
  category.names = c("Symptoms", "CaseControl"),
  fill = c("green", "red"),
  alpha = 0.5,
  filename = NULL,
  scale = FALSE,
  cat.pos = -5
)

lapply(v, function(i) i$label)

# Intersection:
v[[7]]$label <- paste(intersect(cc$V1, sym$V1), collapse="\n")  

# Plot:
grid.newpage()
grid.draw(v)
########################################


# Previous Results vs Diameter:
########################################
v <- venn.diagram(
  x = list(A = dia$V1, B = med$Gene.Symbol, C = adv$Gene.Symbol),
  category.names = c("Diameter", "Media", "Adventitia"),
  fill = c("yellow", "green", "blue"),
  alpha = 0.5,
  filename = NULL,
  scale = FALSE
)

grid.draw(v)

venn.plot <- venn.diagram(
  x = list(A = dia$V1, B = med$Gene.Symbol, C = adv$Gene.Symbol),
  filename = NULL,
  category.names = c("New","Media", "Adventitia"),
  fill = c("yellow", "green", "blue"),
  alpha = 0.5)


lapply(v, function(i) i$label)

# Intesection:
v[[7]]$label <- paste(intersect(cc$V1, dia$V1), collapse="\n")  

# Plot  
grid.newpage()
grid.draw(v)
########################################


# CaseControl vs Diameter vs Symptoms:
########################################
v <- venn.diagram(
  x = list(A = cc$V1, B = dia$V1, C = sym$V1),
  category.names = c("CaseControl", "Diameter", "Symptoms"),
  fill = c("red", "yellow","green"),
  alpha = 0.5,
  filename = NULL,
  cat.pos = c(-10,10,0)
)


lapply(v, function(i) i$label)

v[[10]]$label <- paste(intersect(intersect(cc$V1, dia$V1), sym$V1), collapse="\n")  
v[[9]]$label <- paste(intersect(cc$V1, dia$V1)[!intersect(cc$V1, dia$V1) == "LAMA2"], 
                      collapse="\n")
v[[11]]$label <- paste(intersect(cc$V1, sym$V1)[!intersect(cc$V1, sym$V1) == "LAMA2"],
                      collapse="\n")

grid.newpage()
grid.draw(v)
########################################


# CaseControl vs Diameter vs Symptoms:
########################################
v <- venn.diagram(
  x = list(A = alt$gene_name, B = cc$V1),
  category.names = c("Splicing", "CaseControl"),
  fill = c("blue", "red"),
  alpha = 0.5,
  filename = "C://Users/Gerard/Desktop/AAA/RNAseq/Venn/CaseControl_Splicing.png",
  cat.pos = c(-30, 1)
)


lapply(v, function(i) i$label)

v[[10]]$label <- paste(intersect(intersect(cc$V1, dia$V1), sym$V1), collapse="\n")  

grid.newpage()
grid.draw(v)
########################################














