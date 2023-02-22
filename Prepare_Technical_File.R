library(openxlsx)
library(tidyr)

tec <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/Tables/Libraries_Unique.xlsx")
tec

cr <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/Tables/Correspondencias.xlsx")
cr

# Add name:
tec2 <- merge(tec, cr, by.x = "Sample.N", by.y = "NIM_ID")
tec2

names(tec2)[names(tec2) == "Sample.N"] <- "SampleID"
names(tec2)[names(tec2) == "EXHEUS_ID"] <- "Seq_ID"

cc <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/Tables/CaseControl.xlsx")
cc

tec3 <- merge(tec2, cc, by.x = "Seq_ID", by.y = "SantPau_ID")
tec3


rin <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/Tables/RIN.xlsx")
rin

tec4 <- merge(tec3, rin[,c("Sample_ID","RIN","DV200","Qubit")], by.x = "Seq_ID", by.y = "Sample_ID")
tec4

gc1 <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/Tables/MultiQC1.xlsx")
gc2 <- read.xlsx("C://Users/Gerard/Desktop/AAA/RNAseq/Tables/MultiQC2.xlsx")

names(gc1)[names(gc1) == "%.GC"] <- "GC_Mean"
names(gc2)[names(gc2) == "%.GC"] <- "GC_Mean"

gc1$Sample.Name <- gsub("_.*", "", gc1$Sample.Name)
gc2$Sample.Name <- gsub("_.*", "", gc2$Sample.Name)

gc1 <- fill(gc1, "GC_Mean", .direction = "up")
gc2 <- fill(gc2, "GC_Mean", .direction = "up")

first_occurrence1 <- !duplicated(gc1$Sample.Name)
first_occurrence2 <- !duplicated(gc2$Sample.Name)

gc1 <- gc1[first_occurrence1, ]
gc2 <- gc2[first_occurrence2, ]

tec5 <- merge(tec4, gc1[,c("Sample.Name","GC_Mean")], all.x = T, by.x = "SampleID", by.y = "Sample.Name")
tec6 <- merge(tec5, gc2[,c("Sample.Name","GC_Mean")], all.x = T, by.x = "SampleID", by.y = "Sample.Name")

tec6$GC_Mean <- ifelse(is.na(tec6$GC_Mean.x), tec6$GC_Mean.y, tec6$GC_Mean.x)

tec7 <- tec6[,c(1:11,14)]
tec7

write.xlsx(tec7, "C://Users/Gerard/Desktop/AAA/RNAseq/Tables/Technical_Complete.xlsx")

