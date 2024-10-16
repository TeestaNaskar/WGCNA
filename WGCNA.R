#this R script is for WGCNA analysis as well as the visualization
setwd('/Users/teestanaskar/Dropbox/Teesta/Placenta/Rat_placenta/RNAseq/MEGENA/')
library(MEGENA)
library(Matrix)
library(openxlsx)
library(matrixStats)
library(WGCNA) ## Eigengene calculation uses WGCNA code
library(ComplexHeatmap)
library(ggplot2)
#loading VST normalized count data from data inventory
vst.rat = read.csv("../DEG_STARalignedcounts/bothsex.rat.VST_counts.csv", header=TRUE,check.names=FALSE)
rownames(vst.rat) = toupper(vst.rat[,1])
vst.rat=vst.rat[2:ncol(vst.rat)]

#IMPORT METADATA
rat.meta=read.table("../data/metadata.txt", header = TRUE, check.names = FALSE)
rownames(rat.meta) = rat.meta$ID
rownames(rat.meta) = gsub("X","",rownames(rat.meta))
rat.meta = rat.meta[2:ncol(rat.meta)]

rat.meta.sub=rat.meta[rat.meta$Group== "THC_CBD",]
#load metadata for ecb rats
rat.ecb = read.xlsx("../../../../byJMF/rat_ecbs.xlsx")
rownames(rat.ecb) = rat.ecb$ID
#remove rows that has no values for any of the measures
rat.ecb = rat.ecb[!(rat.ecb$AEA == 'NA'),]
rat.ecb = rat.ecb[!(rownames(rat.ecb) == 'NA'), ]
vst.rat.sub=vst.rat[,as.character(rat.ecb$ID)]
print(dim(vst.rat.sub))

#remove samples from the vst count data that are not in the rat.ecb
animals_to_remove = setdiff(rownames(rat.ecb), colnames(vst.rat))
vst.rat.sub = vst.rat[!(rownames(rat.ecb) %in% animals_to_remove),]
setwd("Rat.Placenta/THC_CBD/")
#remove genes with stddev=0
sd_rows.rat <- apply(vst.rat.sub, 1, sd)
vst.rat.sub = vst.rat.sub[which(sd_rows>0),]
print(dim(vst.rat.sub))
set.seed(12345)
library(msigdbi)
library(GOtest)
require(WGCNA)
#running WGCNA with megena's sugnificant modules
#load significant modules of MEGENA
sigmod = read.table("output/significant_module_2column_table.txt", header = T)
vst.rat.sub=vst.rat.sub[match(sigmod$values,rownames(vst.rat.sub)),rat.ecb$Genotype=="THC_CBD"]
rat.meta=rat.ecb[rat.ecb$Genotype=="THC_CBD",]

#int= c("AEA", "OEA", "AA", "2-AG+1AG")
ecb = rat.meta[,4:7]
dim(ecb)
sum(complete.cases(ecb))
modE=moduleEigengenes(t(vst.rat.sub),sigmod$ind)
ME0 = modE$eigengenes
MEs = orderMEs(ME0)
moduleTraitCor = cor(MEs, ecb, use = "p",method="spearman") 
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, 20) #20 = number of subjects
View(moduleTraitCor)
View(moduleTraitPvalue)
moduleTraitCor = data.frame(moduleTraitCor)
moduleTraitCor$ID = rownames(moduleTraitCor)

#writing module trait correlation and pvalue
write.xlsx(moduleTraitCor, "moduleTraitcor.xlsx")

moduleTraitPvalue = data.frame(moduleTraitPvalue)
moduleTraitPvalue$ID = rownames(moduleTraitPvalue)
write.xlsx(moduleTraitPvalue, "moduleTraitPvalue.xlsx")
#putting aterisk in the modules
sig_asterisks <- matrix("",nrow(moduleTraitCor),ncol(moduleTraitCor))
ast_idx = which(moduleTraitPvalue<0.05, arr.ind=T)
for (i in 1:nrow(ast_idx)){sig_asterisks[ast_idx[i,1], ast_idx[i,2]] <- "*"}

#visualization of the cormat structure
cormat = as.matrix(moduleTraitCor)
library(circlize)
f1 = colorRamp2(seq(min(cormat), max(cormat), length = 7), c("blue", "purple", "cyan","yellow", "#FF7F00", "red", "#7F0000"), space = "RGB")
Heatmap(cormat, col = f1,rect_gp = gpar(col = "gray", lwd = 1), row_dend_reorder = TRUE,row_names_gp = gpar(fontsize = 16))

#human eigengene modules with correlation value
ratEG = moduleTraitCor
head(rownames(ratEG))
rownames(ratEG) = gsub("ME", "",rownames(ratEG))

#load significant modules that are overlapped with DEGs
Rat.DEG.overlapped.mod = read.xlsx("output/rat.MEGENA_mod_overlap_DEGs.xlsx")
Rat.DEG.overlapped.mod.0.05 = Rat.DEG.overlapped.mod[Rat.DEG.overlapped.mod$Pvalue<0.05, ]
#find eigenegene modules that are sig for rat deg overlapped modules
R.EG.DE.ME = Rat.DEG.overlapped.mod.0.05 %>%
  filter(Input %in% rownames(ratEG))
#next considering modules that are above 2 overlapped sizes or number of genes
R.EG.DE.ME = R.EG.DE.ME[R.EG.DE.ME$Overlap.Size >2,]

GO.rat = read.xlsx("output/rat.THC.MEGENA_mod_GO.xlsx")
#find the eigenegene modules that have GO ontologies and also significant for overlapped DEG modules
R.EG.DE.ME.GO = R.EG.DE.ME %>%
  filter(Input %in% c(GO.rat$Input))

ratEG = ratEG[R.EG.DE.ME.GO$Input,]
#for writing
ratEG$ID = rownames(ratEG)
write.xlsx(ratEG,"output/Rat.EG.ME.DE.GO.xlsx")
