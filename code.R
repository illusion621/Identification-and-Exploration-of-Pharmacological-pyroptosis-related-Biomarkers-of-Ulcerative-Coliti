library(limma)
library(ggThemeAssist)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(pheatmap)
library(RColorBrewer)
library(GOplot)
library(clusterProfiler)
library(ggplot2)
library(AnnotationDbi)
library(AnnotationHub)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(gplots)
library(data.table)
library(GenVisR)
library(descriptr)
library(vioplot)
library(survival)
library(corrplot)
library(preprocessCore)
source("Cibersort.R")
library(e1071)
library(psych)
library(ggpubr)
library(ggcorrplot)
library(WGCNA)
library(VennDiagram)
library(enrichplot)
library(pROC)
library(pacman)
library(ggsignif)
library(viridis)
library(Seurat)
library(clustree)
library(cowplot)
library(scHCL)
library(shinythemes)
library(celldex)
library(scRNAseq)
library(SingleR)
library(monocle)
library(AUCell)


wkdir="H:\UC\UC2"
getwd()
GPL13158_anno=data.table::fread("GPL13158-5065.txt",skip = 16) #读取注释，skip表示跳过开头的
GPL13158_anno=GPL13158_anno[,c(1,11)] 
GSE87466_matrix=read.table("GSE87466_series_matrix.txt",header = T,sep="\t",row.names = 1,comment.char = "!")
GSE87466_matrix$ID=rownames(GSE87466_matrix)
GSE87466=merge(x=GSE87466_matrix,y=GPL13158_anno,by="ID",all.x = T)
GSE87466$ID=NULL
rowMeans=apply(GSE87466,1,function(x)mean(as.numeric(x),na.rm=T))
GSE87466=GSE87466[order(rowMeans,decreasing = T),]
GSE87466=GSE87466[!duplicated(GSE87466[,dim(GSE87466)[2]]),]
GSE87466=GSE87466[!grepl("///",GSE87466$`Gene Symbol`),]
rownames(GSE87466)=GSE87466[,dim(GSE87466)[2]]
GSE87466=GSE87466[,-dim(GSE87466)[2]]
write.csv(GSE87466,"GSE87466_Symbol_expr.csv",quote=F)

heatmap(cor(GSE87466),scale ="none",Colv = NA,Rowv = NA)

par(cex = 0.7)
n.sample=ncol(GSE87466)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
boxplot(GSE87466, col = cols,main="expression value",las=2)

#DEG
library(limma)
group=read.csv("DEGfenlei.csv",sep = ",",row.names = 1,header = T,check.names = F)
design=model.matrix(~0+factor(group$type)) 
colnames(design)=levels(factor(group$type))
rownames(design)=colnames(GSE87466)
fit=lmFit(GSE87466,design) 
cont.matrix=makeContrasts(UC-Control,levels = design)
fit2=contrasts.fit(fit,cont.matrix)  
fit2=eBayes(fit2)
temp0output=topTable(fit2,coef=1,n=Inf,adjust="BH",sort.by = "B",resort.by = "M")
foldchange=1
pvalue=0.05
diff=temp0output
diff=diff[(diff$P.Value<pvalue&(diff$logFC>foldchange | diff$logFC<(-foldchange))),]

x=read.csv("limmaOut.csv",row.names = 1)
x$label=rownames(x)
logFCcut=1
pvaluecut=0.05
x[,7]=ifelse((x$P.Value<pvaluecut & x$logFC>logFCcut),"red",ifelse((x$P.Value<0.05 & x$logFC<(-logFCcut)),"blue","grey30"))
size=ifelse((x$P.Value<pvaluecut & abs(x$logFC)>logFCcut),4,2)
xmin=-3
xmax=3
ymax=45
ymin=0
valo=ggplot(data=x,aes(x=logFC,y=-log10(P.Value),label=label)) +
  geom_point(alpha=0.6,size=size,colour=x[,7]) +
  scale_color_manual(values=c("lightgrey","navy","red")) +
  labs(x=bquote(~log[2]~"(fold Change)"),y=bquote(~-log[10]~italic("p-value"))) +
  ylim(c(ymin,ymax)) +
  scale_x_continuous(
    breaks = c(-5,-3,-logFCcut,0,logFCcut,3,5),
    labels = c(-5,-3,-logFCcut,0,logFCcut,3,5),
    limits = c(-6,6)
  ) +
  geom_vline(xintercept = logFCcut,color="grey40",linetype="longdash",size=0.5) +
  geom_vline(xintercept = -logFCcut,color="grey40",linetype="longdash",size=0.5) +
  geom_hline(yintercept = -log10(pvaluecut),color="grey40",linetype="longdash",size=0.5) +
  guides(colour=guide_legend(override.aes = list(shape=16))) +
  theme_bw(base_size=12,base_family = "Times") +
  theme(legend.position = "right",
        panel.grid=element_blank(),
        legend.title=element_blank(),
        legend.text=element_text(face="bold",color="black",family = "Times_New_Roman",size=8),
        plot.title=element_text(hjust=0.8),
        axis.text.x=element_text(face="bold",color="black",size=15),
        axis.text.y=element_text(face="bold",color="black",size=15))

allp=GSE87466[rownames(diff),]
pall=pheatmap(allp,cluster_cols = F,cluster_rows=F,scale="row",color = colorRampPalette(c("navy","white","firebrick3"))(100),
              main = "DEGheatmap",annotation_col = group,show_colnames = F,show_rownames = F,
              annotation_row = annotation_row)

#GO KEGG
diff$SYMBOL=rownames(diff)
diff_id=bitr(diff$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
diff_id=inner_join(diff,diff_id,by="SYMBOL")
diff_BP=enrichGO(gene=diff_id$ENTREZID,
                 OrgDb = "org.Hs.eg.db",
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05)
diff_BP=simplify(diff_BP,cutoff=0.7,by="p.adjust",select_fun=min)

diff_kegg=enrichKEGG(gene=diff_id$ENTREZID,
                     keyType = "kegg",
                     organism = "hsa",
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",)

#Cibersort
cs_result=CIBERSORT("LM22.txt","GSE87466.txt",perm = 100,QN=T) 
write.csv(cs_result,"cs_result.csv")
cs_result=read.csv("cs_result.csv",header = T,row.names = 1)
cs_result=cs_result[(cs_result$P.value<0.05),] 
cs_result_cl=cs_result[,-c(23,24,25)] 
cs_result_t=t(cs_result_cl)

corrplot(corr=cor(cs_result_cl),
         method="ellipse",
         order="hclust",
         tl.col="black",
         addCoef.col = "blue",
         number.cex = 0.7,
         tl.srt = 45,
         col=colorRampPalette(c("green","white","red"))(100))

cs_result_xx=cs_result_xx%>%rownames_to_column("sample")
cs_result_xx$Group=NULL
cs_result_xx$Group=c(rep("Control",21),rep("UC",86))
cs_result_plot=gather(cs_result_xx,key = CIBERSORT,value=Proportion,-c(Group,sample))
ggboxplot(cs_result_plot,x="CIBERSORT",y="Proportion",
          fill="Group",palette = "Lancet") +
  stat_compare_means(aes(group=Group),
                     method = "t.test",
                     label="p.signif",
                     symnum.args = list(cutpoints=c(0,0.001,0.01,0.05,1),
                                        symbols=c("***","**","*"," "))) +
  theme(text=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust = 1))

sig_gene <- c("IL1B","IL1A","NLRP7","IL18R1","ZBP1","GZMB","TNF","IL6","TREM1","AIM2")
xi=as.data.frame(t(GSE87466))
xi=xi[-18,]
xi <- xi[,sig_gene]
yi <- cs_result_cl
xi=xi[-c(1:21),]
yi=yi[-c(1:20),]
di <- corr.test(xi,yi,use="complete",method = 'spearman')

ri <- di$r
pi <- di$p
ggcorrplot(t(di$r), show.legend = T, 
           p.mat = t(di$p.adj), digits = 1,  sig.level = 0.05,insig = 'blank',lab = T)


#WGCNA
WGCNA=GSE87466
WGCNA$CV=apply(WGCNA,1,function(x) sd(x)/mean(x))
WGCNA=WGCNA[(WGCNA$CV>=0.08),]
WGCNA$CV=NULL
write.csv(WGCNA,"WGCNA.csv")

WGCNA=as.data.frame(t(WGCNA))
sampleTree = hclust(dist(WGCNA), method = "average");
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
WGCNA=WGCNA[-47,]
WGCNA=WGCNA[-73,]
powers = c(c(1:10), seq(from = 11, to=20, by=1))
sft = pickSoftThreshold(WGCNA, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");

plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

cor <- WGCNA::cor
net = blockwiseModules(WGCNA, power =12 ,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.2,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,saveTOMFileBase = "femaleMouseTOM", 
                       verbose = 3,)
cor<-stats::cor 
table(net$colors)
sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
nGenes = ncol(WGCNA);
nSamples = nrow(WGCNA);
MEs0 = moduleEigengenes(WGCNA, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs,group_WGCNA, use = "p");
#group_WGCNA=read.csv("group_wgcna.csv",sep=",",header = T,row.names = 1)
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

par(mar = c(9, 10, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(group_WGCNA),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

geneModuleMembership = as.data.frame(cor(WGCNA, MEs, use = "p"));

turquoise_gene=names(WGCNA)[moduleColors=="turquoise"]
turquoise_gene=as.data.frame(turquoise_gene)

hub_UC<- abs(geneModuleMembership$MEturquoise)>0.8 & abs(geneTraitSignificance)>0.2
hub_UC=as.data.frame(hub_UC)
hub_UC$turquoise_gene=rownames(hub_UC)
hub_UC=merge(hub_UC,turquoise_gene,by="turquoise_gene",all.x=F)
rownames(hub_UC)=hub_UC$turquoise_gene
hub_UC=hub_UC[hub_UC$`group_WGCNA$UC`==TRUE,]

modNames = substring(names(MEs),3)
module = "turquoise"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance_nc[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

#VENN
Pyroptosis=read.csv("焦亡.csv",row.names = 1,header = T,sep=",")
turquoise=names(WGCNA)[moduleColors=="turquoise"]
List_ID=list("DEG"=rownames(diff),
             "turquoise"=turquoise,
             "Pyroptosis"=rownames(Pyroptosis)
)

venn.plot=venn.diagram(x=List_ID,filename="Venn.pdf",
                       col="black",lwd=2,fontface="bold",
                       fill = c("cornflowerblue", "turquoise", "red"))

inter <- get.venn.partitions(List_ID)

for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')

write.table(inter[-c(5, 6)], 'vennup_inter.txt', row.names = FALSE, sep = '\t', quote = FALSE)


List_hub=list("DEG"=rownames(diff),
              "turquoise_hub"=rownames(hub_UC),
              "Pyroptosis"=rownames(Pyroptosis)
)
venn.plot=venn.diagram(x=List_hub,filename="Venn_hub.png",
                       col="black",lwd=2,fontface="bold",
                       fill = c("cornflowerblue", "turquoise", "red"))

inter_hub <- get.venn.partitions(List_hub)

for (i in 1:nrow(inter_hub)) inter_hub[i,'values'] <- paste(inter_hub[[i,'..values..']], collapse = ', ')

write.table(inter_hub[-c(5, 6)], 'vennhub_inter.txt', row.names = FALSE, sep = '\t', quote = FALSE)



ROC cruve and Hub gene correlation coefficient completed by sangerbox

# Verification
GSE92415_matrix=read.table("GSE92415_series_matrix.txt",header = T,sep="\t",row.names = 1,comment.char = "!")
GSE92415_matrix=GSE92415_matrix[,-c(1:109)]
GPL13158_anno=data.table::fread("GPL13158-5065.txt",skip = 16)
GPL13158_anno=GPL13158_anno[,c(1,11)]
GSE92415_matrix$ID=rownames(GSE92415_matrix)
GSE92415=merge(x=GSE92415_matrix,y=GPL13158_anno,by="ID",all.x = T)
GSE92415$ID=NULL
rowMeans=apply(GSE92415,1,function(x)mean(as.numeric(x),na.rm=T))
GSE92415=GSE92415[order(rowMeans,decreasing = T),]
GSE92415=GSE92415[!duplicated(GSE92415[,dim(GSE92415)[2]]),]
GSE92415=GSE92415[!grepl("///",GSE92415$`Gene Symbol`),]
rownames(GSE92415)=GSE92415[,dim(GSE92415)[2]]
GSE92415=GSE92415[,-dim(GSE92415)[2]]
write.csv(GSE92415,"GSE92415_Symbol_expr.csv",quote=F)

GSE92415_yanzhenguc=GSE92415_uc[,c("IL1B","IL1A","TNF","TREM1","IL6","NLRP7","GZMB","ZBP1","AIM2","IL18R1")]
GSE92415_control=GSE92415_control[,c("IL1B","IL1A","TNF","TREM1","IL6","NLRP7","GZMB","ZBP1","AIM2","IL18R1")]
write.csv(GSE92415_yanzhenguc,"GSE92415_yanzhenguc.csv")

GSE92415_control=as.data.frame(t(GSE92415[,c(1:21)]))
GSE92415_uc=as.data.frame(t(GSE92415[,-c(1:21)]))
AIM2_control=GSE92415_control$AIM2
AIM2_uc=GSE92415_uc$AIM2
Group=c(rep("Control",21),rep("UC",53))
AIM2=c(AIM2_control,AIM2_uc)
data_AIM2 <- data.frame(Group,AIM2)
compaired <- list(c("Control","UC"))
ggplot(data_AIM2) +
  aes(x=Group,y=AIM2) +
  geom_boxplot(alpha=.25) +
  geom_jitter(alpha=.5,
              height=0,
              width=.25) +
  aes(color=names,fill=names) +
  theme_bw() +
  stat_boxplot(geom="errorbar",width=0.5) +
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test,colour="black")

IL1B_control=GSE92415_control$IL1B
IL1B_uc=GSE92415_uc$IL1B
IL1B=c(IL1B_control,IL1B_uc)
data_IL1B <- data.frame(Group,IL1B)
ggplot(data_IL1B) +
  aes(x=Group,y=IL1B) +
  geom_boxplot(alpha=.25) +
  geom_jitter(alpha=.5,
              height=0,
              width=.25) +
  aes(color=names,fill=names) +
  theme_bw() +
  stat_boxplot(geom="errorbar",width=0.5) +
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test,colour="black")

IL6_control=GSE92415_control$IL6
IL6_uc=GSE92415_uc$IL6
IL6=c(IL6_control,IL6_uc)
data_IL6 <- data.frame(Group,IL6)
ggplot(data_IL6) +
  aes(x=Group,y=IL6) +
  geom_boxplot(alpha=.25) +
  geom_jitter(alpha=.5,
              height=0,
              width=.25) +
  aes(color=names,fill=names) +
  theme_bw() +
  stat_boxplot(geom="errorbar",width=0.5) +
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test,colour="black")

NLRP7_control=GSE92415_control$NLRP7
NLRP7_uc=GSE92415_uc$NLRP7
NLRP7=c(NLRP7_control,NLRP7_uc)
data_NLRP7 <- data.frame(Group,NLRP7)
ggplot(data_NLRP7) +
  aes(x=Group,y=NLRP7) +
  geom_boxplot(alpha=.25) +
  geom_jitter(alpha=.5,
              height=0,
              width=.25) +
  aes(color=names,fill=names) +
  theme_bw() +
  stat_boxplot(geom="errorbar",width=0.5) +
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test,colour="black")

TNF_control=GSE92415_control$TNF
TNF_uc=GSE92415_uc$TNF
TNF=c(TNF_control,TNF_uc)
data_TNF <- data.frame(Group,TNF)
ggplot(data_TNF) +
  aes(x=Group,y=TNF) +
  geom_boxplot(alpha=.25) +
  geom_jitter(alpha=.5,
              height=0,
              width=.25) +
  aes(color=names,fill=names) +
  theme_bw() +
  stat_boxplot(geom="errorbar",width=0.5) +
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test,colour="black")


IL1A_control=GSE92415_control$IL1A
IL1A_uc=GSE92415_uc$IL1A
IL1A=c(IL1A_control,IL1A_uc)
data_IL1A <- data.frame(Group,IL1A)
ggplot(data_IL1A) +
  aes(x=Group,y=IL1A) +
  geom_boxplot(alpha=.25) +
  geom_jitter(alpha=.5,
              height=0,
              width=.25) +
  aes(color=names,fill=names) +
  theme_bw() +
  stat_boxplot(geom="errorbar",width=0.5) +
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test,colour="black")

IL18R1_control=GSE92415_control$IL18R1
IL18R1_uc=GSE92415_uc$IL18R1
IL18R1=c(IL18R1_control,IL18R1_uc)
data_IL18R1 <- data.frame(Group,IL18R1)
ggplot(data_IL18R1) +
  aes(x=Group,y=IL18R1) +
  geom_boxplot(alpha=.25) +
  geom_jitter(alpha=.5,
              height=0,
              width=.25) +
  aes(color=names,fill=names) +
  theme_bw() +
  stat_boxplot(geom="errorbar",width=0.5) +
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test,colour="black")

ZBP1_control=GSE92415_control$ZBP1
ZBP1_uc=GSE92415_uc$ZBP1
ZBP1=c(ZBP1_control,ZBP1_uc)
data_ZBP1 <- data.frame(Group,ZBP1)
ggplot(data_ZBP1) +
  aes(x=Group,y=ZBP1) +
  geom_boxplot(alpha=.25) +
  geom_jitter(alpha=.5,
              height=0,
              width=.25) +
  aes(color=names,fill=names) +
  theme_bw() +
  stat_boxplot(geom="errorbar",width=0.5) +
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test,colour="black")

GZMB_control=GSE92415_control$GZMB
GZMB_uc=GSE92415_uc$GZMB
GZMB=c(GZMB_control,GZMB_uc)
data_GZMB <- data.frame(Group,GZMB)
ggplot(data_GZMB) +
  aes(x=Group,y=GZMB) +
  geom_boxplot(alpha=.25) +
  geom_jitter(alpha=.5,
              height=0,
              width=.25) +
  aes(color=names,fill=names) +
  theme_bw() +
  stat_boxplot(geom="errorbar",width=0.5) +
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test,colour="black")

GSE107499_matrix=read.table("GSE107499_series_matrix.txt",header = T,sep="\t",row.names = 1,comment.char = "!")
GPL15207_anno=data.table::fread("GPL15207-17536.txt",skip = 36)
GPL15207_anno=GPL15207_anno[,c(1,17)]
GSE107499_matrix$ID=rownames(GSE107499_matrix)
GSE107499=merge(x=GSE107499_matrix,y=GPL15207_anno,by="ID",all.x = T)
GSE107499$ID=NULL
rowMeans=apply(GSE107499,1,function(x)mean(as.numeric(x),na.rm=T))
GSE107499=GSE107499[order(rowMeans,decreasing = T),]
GSE107499=GSE107499[!duplicated(GSE107499[,dim(GSE107499)[2]]),]
GSE107499=GSE107499[!grepl("///",GSE107499$`Gene Symbol`),]
rownames(GSE107499)=GSE107499[,dim(GSE107499)[2]]
GSE107499=GSE107499[,-dim(GSE107499)[2]]
write.csv(GSE107499,"GSE107499.csv")

GENE_NonLesional=as.numeric(c(GSE107499_Nonlesional$AIM2,GSE107499_Nonlesional$IL1B,GSE107499_Nonlesional$NLRP7,
                              GSE107499_Nonlesional$IL6,GSE107499_Nonlesional$TNF,GSE107499_Nonlesional$IL1A,
                              GSE107499_Nonlesional$IL18R1,GSE107499_Nonlesional$ZBP1,GSE107499_Nonlesional$GZMB,GSE107499_Nonlesional$TREM1))
z=c(rep("AIM2",22),rep("IL1B",22),rep("NLRP7",22),rep("IL6",22),rep("TNF",22),rep("IL1A",22),rep("IL18R1",22),
    rep("ZBP1",22),rep("GZMB",22),rep("TREM1",22))
NonLesional=as.data.frame(z)
NonLesional$value=GENE_NonLesional
write.csv(NonLesional,"NonLesional.csv")

GENE_Lesional=as.numeric(c(GSE107499_Lesional$AIM2,GSE107499_Lesional$IL1B,GSE107499_Lesional$NLRP7,
                           GSE107499_Lesional$IL6,GSE107499_Lesional$TNF,GSE107499_Lesional$IL1A,
                           GSE107499_Lesional$IL18R1,GSE107499_Lesional$ZBP1,GSE107499_Lesional$GZMB,GSE107499_Lesional$TREM1))
z1=c(rep("AIM2",97),rep("IL1B",97),rep("NLRP7",97),rep("IL6",97),rep("TNF",97),rep("IL1A",97),rep("IL18R1",97),
     rep("ZBP1",97),rep("GZMB",97),rep("TREM1",97))
Lesional=as.data.frame(z1)
Lesional$value=GENE_Lesional
write.csv(Lesional,"Lesional.csv")

#Part of the visualization work is done by sangerbox


GSE178753_matrix=read.table("GSE178753_series_matrix.txt",header = T,sep="\t",row.names = 1,comment.char = "!")
GPL30304_anno=data.table::fread("GPL30304_family.soft.gz",skip = 46)
GPL30304_anno=GPL30304_anno[,c(1,2)]
GSE178753_matrix$ID=rownames(GSE178753_matrix)
refGene=data.table::fread("refGene.txt")
refGene=refGene[,c(2,13)]
refGene$GB_ACC=refGene$V2
refGene$SYMBOL=refGene$V13
refGene$V2=NULL
refGene$V13=NULL
GPL30304_anno=merge(x=GPL30304_anno,y=refGene,by="GB_ACC",all.x = F)
GSE178753=merge(x=GSE178753_matrix,y=GPL30304_anno,by="ID",all.x = F)
GSE178753$ID=NULL
GSE178753$GB_ACC=NULL
rowMeans=apply(GSE178753,1,function(x)mean(as.numeric(x),na.rm=T))
GSE178753=GSE178753[order(rowMeans,decreasing = T),]
GSE178753=GSE92415[!duplicated(GSE178753[,dim(GSE178753)[2]]),]
GSE178753=GSE92415[!grepl("///",GSE178753$SYMBOL),]
rownames(GSE178753)=GSE178753[,dim(GSE178753)[2]]
GSE178753=GSE178753[,-dim(GSE178753)[2]]

#The following CSV data were obtained from the processed GSE178753 and the corresponding gene expression data were extracted
Rain_AIM2=read.csv("AIM2.csv")
Rain_AIM2AOV<-aov(Expression~Group,data=Rain_AIM2)
summary(Rain_AIM2AOV)
TukeyHSD(Rain_AIM2AOV)
Rain_AIM2 %>% ggplot(aes(Group,Expression,fill=Group)) +
  geom_half_violin(position = position_nudge(x = 0.2),side=2,alpha = 0.8) +
  geom_point(aes(y=Expression, color = Group), 
             position = position_jitter(width = 0.15),
             size = 1,alpha = 0.8)+labs(x=NULL) +
  geom_boxplot(width = 0.1, outlier.shape = NA,alpha = 0.8) +
  coord_flip() +
  theme_bw() +
  xlab("Experimental group")+ylab("AIM2") +
  geom_signif(annotations = c("***","***"), y_position = c(8,8),xmin = c(2,2), xmax = c(1,3))


Rain_IL1B=read.csv("IL1B.csv")
Rain_IL1BAOV<-aov(Expression~Group,data=Rain_IL1B)
summary(Rain_IL1BAOV)
TukeyHSD(Rain_IL1BAOV)
Rain_IL1B %>% ggplot(aes(Group,Expression,fill=Group)) +
  geom_half_violin(position = position_nudge(x = 0.2),side=2,alpha = 0.8) +
  geom_point(aes(y=Expression, color = Group), 
             position = position_jitter(width = 0.15),
             size = 1,alpha = 0.8)+labs(x=NULL) +
  geom_boxplot(width = 0.1, outlier.shape = NA,alpha = 0.8) +
  coord_flip() +
  theme_bw() +
  xlab("Experimental group")+ylab("IL1B") +
  geom_signif(annotations = c("***","***","*"), y_position = c(11.5,11.5,12.5),xmin = c(2,2,1), xmax = c(1,3,3))


Rain_TNF=read.csv("TNF.csv")
Rain_TNFAOV<-aov(Expression~Group,data=Rain_TNF)
summary(Rain_TNFAOV)
TukeyHSD(Rain_TNFAOV)
Rain_TNF %>% ggplot(aes(Group,Expression,fill=Group)) +
  geom_half_violin(position = position_nudge(x = 0.2),side=2,alpha = 0.8) +
  geom_point(aes(y=Expression, color = Group), 
             position = position_jitter(width = 0.15),
             size = 1,alpha = 0.8)+labs(x=NULL) +
  geom_boxplot(width = 0.1, outlier.shape = NA,alpha = 0.8) +
  coord_flip() +
  theme_bw() +
  xlab("Experimental group")+ylab("TNF") +
  geom_signif(annotations = c("***","***"), y_position = c(8,8),xmin = c(2,2), xmax = c(1,3))



Rain_NLRP7=read.csv("NLRP7.csv")
Rain_NLRP7AOV<-aov(Expression~Group,data=Rain_NLRP7)
summary(Rain_NLRP7AOV)
TukeyHSD(Rain_NLRP7AOV)
Rain_NLRP7 %>% ggplot(aes(Group,Expression,fill=Group)) +
  geom_half_violin(position = position_nudge(x = 0.2),side=2,alpha = 0.8) +
  geom_point(aes(y=Expression, color = Group), 
             position = position_jitter(width = 0.15),
             size = 1,alpha = 0.8)+labs(x=NULL) +
  geom_boxplot(width = 0.1, outlier.shape = NA,alpha = 0.8) +
  coord_flip() +
  theme_bw() +
  xlab("Experimental group")+ylab("NLRP7") +
  geom_signif(annotations = c("***","***"), y_position = c(6,6),xmin = c(2,2), xmax = c(1,3))

Rain_IL6=read.csv("IL6.csv")
Rain_IL6AOV<-aov(Expression~Group,data=Rain_IL6)
summary(Rain_IL6AOV)
TukeyHSD(Rain_IL6AOV)
Rain_IL6 %>% ggplot(aes(Group,Expression,fill=Group)) +
  geom_half_violin(position = position_nudge(x = 0.2),side=2,alpha = 0.8) +
  geom_point(aes(y=Expression, color = Group), 
             position = position_jitter(width = 0.15),
             size = 1,alpha = 0.8)+labs(x=NULL) +
  geom_boxplot(width = 0.1, outlier.shape = NA,alpha = 0.8) +
  coord_flip() +
  theme_bw() +
  xlab("Experimental group")+ylab("IL6") +
  geom_signif(annotations = c("***","***"), y_position = c(11.5,11.5),xmin = c(2,2), xmax = c(1,3))

Rain_IL1A=read.csv("IL1A.csv")
Rain_IL1AAOV<-aov(Expression~Group,data=Rain_IL1A)
summary(Rain_IL1AAOV)
TukeyHSD(Rain_IL1AAOV)
Rain_IL1A %>% ggplot(aes(Group,Expression,fill=Group)) +
  geom_half_violin(position = position_nudge(x = 0.2),side=2,alpha = 0.8) +
  geom_point(aes(y=Expression, color = Group), 
             position = position_jitter(width = 0.15),
             size = 1,alpha = 0.8)+labs(x=NULL) +
  geom_boxplot(width = 0.1, outlier.shape = NA,alpha = 0.8) +
  coord_flip() +
  theme_bw() +
  xlab("Experimental group")+ylab("NLRP7") +
  geom_signif(annotations = c("***","***"), y_position = c(10,10),xmin = c(2,2), xmax = c(1,3))

Rain_IL18R1=read.csv("IL18R1.csv")
Rain_IL18R1AOV<-aov(Expression~Group,data=Rain_IL18R1)
summary(Rain_IL18R1AOV)
TukeyHSD(Rain_IL18R1AOV)
Rain_IL18R1 %>% ggplot(aes(Group,Expression,fill=Group)) +
  geom_half_violin(position = position_nudge(x = 0.2),side=2,alpha = 0.8) +
  geom_point(aes(y=Expression, color = Group), 
             position = position_jitter(width = 0.15),
             size = 1,alpha = 0.8)+labs(x=NULL) +
  geom_boxplot(width = 0.1, outlier.shape = NA,alpha = 0.8) +
  coord_flip() +
  theme_bw() +
  xlab("Experimental group")+ylab("IL18R1") +
  geom_signif(annotations = c("***","***"), y_position = c(7.5,7.5),xmin = c(2,2), xmax = c(1,3))

Rain_ZBP1=read.csv("ZBP1.csv")
Rain_ZBP1AOV<-aov(Expression~Group,data=Rain_ZBP1)
summary(Rain_ZBP1AOV)
TukeyHSD(Rain_ZBP1AOV)
Rain_ZBP1 %>% ggplot(aes(Group,Expression,fill=Group)) +
  geom_half_violin(position = position_nudge(x = 0.2),side=2,alpha = 0.8) +
  geom_point(aes(y=Expression, color = Group), 
             position = position_jitter(width = 0.15),
             size = 1,alpha = 0.8)+labs(x=NULL) +
  geom_boxplot(width = 0.1, outlier.shape = NA,alpha = 0.8) +
  coord_flip() +
  theme_bw() +
  xlab("Experimental group")+ylab("ZBP1") +
  geom_signif(annotations = c("***","***"), y_position = c(9,9),xmin = c(2,2), xmax = c(1,3))

Rain_GZMB=read.csv("GZMB.csv")
Rain_GZMBAOV<-aov(Expression~Group,data=Rain_GZMB)
summary(Rain_GZMBAOV)
TukeyHSD(Rain_GZMBAOV)
Rain_GZMB %>% ggplot(aes(Group,Expression,fill=Group)) +
  geom_half_violin(position = position_nudge(x = 0.2),side=2,alpha = 0.8) +
  geom_point(aes(y=Expression, color = Group), 
             position = position_jitter(width = 0.15),
             size = 1,alpha = 0.8)+labs(x=NULL) +
  geom_boxplot(width = 0.1, outlier.shape = NA,alpha = 0.8) +
  coord_flip() +
  theme_bw() +
  xlab("Experimental group")+ylab("NLRP7") +
  geom_signif(annotations = c("***","***"), y_position = c(8.5,8.5),xmin = c(2,2), xmax = c(1,3))

Rain_TREM1=read.csv("TREM1.csv")
Rain_TREM1AOV<-aov(Expression~Group,data=Rain_TREM1)
summary(Rain_TREM1AOV)
TukeyHSD(Rain_TREM1AOV)
Rain_TREM1 %>% ggplot(aes(Group,Expression,fill=Group)) +
  geom_half_violin(position = position_nudge(x = 0.2),side=2,alpha = 0.8) +
  geom_point(aes(y=Expression, color = Group), 
             position = position_jitter(width = 0.15),
             size = 1,alpha = 0.8)+labs(x=NULL) +
  geom_boxplot(width = 0.1, outlier.shape = NA,alpha = 0.8) +
  coord_flip() +
  theme_bw() +
  xlab("Experimental group")+ylab("TREM1") +
  geom_signif(annotations = c("***","***"), y_position = c(10,10),xmin = c(2,2), xmax = c(1,3))

Rain_IL18=read.csv("IL18.csv")
Rain_IL18AOV<-aov(Expression~Group,data=Rain_IL18)
summary(Rain_IL18AOV)
TukeyHSD(Rain_IL18AOV)
Rain_IL18 %>% ggplot(aes(Group,Expression,fill=Group)) +
  geom_half_violin(position = position_nudge(x = 0.2),side=2,alpha = 0.8) +
  geom_point(aes(y=Expression, color = Group), 
             position = position_jitter(width = 0.15),
             size = 1,alpha = 0.8)+labs(x=NULL) +
  geom_boxplot(width = 0.1, outlier.shape = NA,alpha = 0.8) +
  coord_flip() +
  theme_bw() +
  xlab("Experimental group")+ylab("IL18") +
  geom_signif(annotations = c("***","***"), y_position = c(10,10),xmin = c(2,2), xmax = c(1,3))


Rain_IL33=read.csv("IL33.csv")
Rain_IL33AOV<-aov(Expression~Group,data=Rain_IL33)
summary(Rain_IL33AOV)
TukeyHSD(Rain_IL33AOV)
Rain_IL33 %>% ggplot(aes(Group,Expression,fill=Group)) +
  geom_half_violin(position = position_nudge(x = 0.2),side=2,alpha = 0.8) +
  geom_point(aes(y=Expression, color = Group), 
             position = position_jitter(width = 0.15),
             size = 1,alpha = 0.8)+labs(x=NULL) +
  geom_boxplot(width = 0.1, outlier.shape = NA,alpha = 0.8) +
  coord_flip() +
  theme_bw() +
  xlab("Experimental group")+ylab("IL33") +
  geom_signif(annotations = c("***","***","*"), y_position = c(7,7,8),xmin = c(2,2,1), xmax = c(1,3,3))

Rain_IL17A=read.csv("IL17A.csv")
Rain_IL17AAOV<-aov(Expression~Group,data=Rain_IL17A)
summary(Rain_IL17AAOV)
TukeyHSD(Rain_IL17AAOV)
Rain_IL17A %>% ggplot(aes(Group,Expression,fill=Group)) +
  geom_half_violin(position = position_nudge(x = 0.2),side=2,alpha = 0.8) +
  geom_point(aes(y=Expression, color = Group), 
             position = position_jitter(width = 0.15),
             size = 1,alpha = 0.8)+labs(x=NULL) +
  geom_boxplot(width = 0.1, outlier.shape = NA,alpha = 0.8) +
  coord_flip() +
  theme_bw() +
  xlab("Experimental group")+ylab("IL17A") +
  geom_signif(annotations = c("***","***"), y_position = c(8,8),xmin = c(2,2), xmax = c(1,3))

Rain_IL8=read.csv("CXCL8.csv")
Rain_IL8AAOV<-aov(Expression~Group,data=Rain_IL8)
summary(Rain_IL8AOV)
TukeyHSD(Rain_IL8AAOV)
Rain_IL8 %>% ggplot(aes(Group,Expression,fill=Group)) +
  geom_half_violin(position = position_nudge(x = 0.2),side=2,alpha = 0.8) +
  geom_point(aes(y=Expression, color = Group), 
             position = position_jitter(width = 0.15),
             size = 1,alpha = 0.8)+labs(x=NULL) +
  geom_boxplot(width = 0.1, outlier.shape = NA,alpha = 0.8) +
  coord_flip() +
  theme_bw() +
  xlab("Experimental group")+ylab("IL8") +
  geom_signif(annotations = c("***","***"), y_position = c(11,11),xmin = c(2,2), xmax = c(1,3))


#5-ASA
SE46451_matrix=read.table("GSE46451_series_matrix.txt",header = T,sep="\t",row.names = 1,comment.char = "!")
GSE46451_matrix=log2(GSE46451_matrix)
GPL10558_anno=data.table::fread("GPL10558-50081.txt",skip = 30) #读取注释，skip表示跳过开头的
GSE46451_matrix$ID=rownames(GSE46451_matrix)
GPL10558_anno=GPL10558_anno[,c(1,13)]
GSE46451=merge(x=GSE46451_matrix,y=GPL10558_anno,by="ID",all.x = F)
GSE46451$ID=NULL
rowMeans=apply(GSE46451,1,function(x)mean(as.numeric(x),na.rm=T))
GSE46451=GSE46451[order(rowMeans,decreasing = T),]
GSE46451=GSE46451[!duplicated(GSE46451[,dim(GSE46451)[2]]),]
GSE46451=na.omit(GSE46451)
rownames(GSE46451)=GSE46451[,dim(GSE46451)[2]]
GSE46451=GSE46451[,-dim(GSE46451)[2]]
write.csv(GSE46451,"GSE46451.csv")


AIM2_5ASA=read.csv("AIM2_5ASA.csv")
AIM2_5ASA=stack(AIM2_5ASA)
names(AIM2_5ASA)=c("Expression","Group")
attach(AIM2_5ASA)
detach(AIM2_5ASA)
tapply(Expression,Group,shapiro.test)
bartlett.test(Expression~Group)
AIM2_5ASAAOV<-aov(Expression~Group,data=AIM2_5ASA)
summary(AIM2_5ASAAOV)
TukeyHSD(AIM2_5ASAAOV)
library(PMCMRplus)
AIM2_5ASACOMPARE=bwsAllPairsTest(Expression~Group,data=AIM2_5ASA)
summary(AIM2_5ASACOMPARE)
AIM2_5ASAplot=AIM2_5ASA%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(AIM2_5ASAplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("AIM2")+ylim(0,10) +
  #geom_signif(annotations = c("***","***"), y_position = c(7.5,8),xmin = c(1, 1), xmax = c(2, 3)) +
  #geom_signif(annotations = c("**"), y_position = c(8.5),xmin = c(1), xmax = c(4)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))




IL1B_5ASA=read.csv("IL1B_5ASA.csv")
IL1B_5ASA=stack(IL1B_5ASA)
names(IL1B_5ASA)=c("Expression","Group")
attach(IL1B_5ASA)
detach(IL1B_5ASA)
tapply(Expression,Group,shapiro.test)
bartlett.test(Expression~Group)
IL1B_5ASAAOV<-aov(Expression~Group,data=IL1B_5ASA)
summary(IL1B_5ASAAOV)
TukeyHSD(IL1B_5ASAAOV)
IL1B_5ASAplot=IL1B_5ASA%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(IL1B_5ASAplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("IL1B")+ylim(0,15) +
  geom_signif(annotations = c("**","***"), y_position = c(12,14),xmin = c(1, 1), xmax = c(2, 4)) +
  geom_signif(annotations = c("**"), y_position = c(13),xmin = c(1), xmax = c(3)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))

TNF_5ASA=read.csv("TNF_5ASA.csv")
TNF_5ASA=stack(TNF_5ASA)
names(TNF_5ASA)=c("Expression","Group")
TNF_5ASAAOV<-aov(Expression~Group,data=TNF_5ASA)
summary(TNF_5ASAAOV)
TukeyHSD(TNF_5ASAAOV)
TNF_5ASAplot=TNF_5ASA%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(TNF_5ASAplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("TNF")+ylim(0,10) +
  #geom_signif(annotations = c("**","***"), y_position = c(12,14),xmin = c(1, 1), xmax = c(2, 4)) +
  #geom_signif(annotations = c("**"), y_position = c(13),xmin = c(1), xmax = c(3)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))


NLRP7_5ASA=read.csv("NLRP7_5ASA.csv")
NLRP7_5ASA=stack(NLRP7_5ASA)
names(NLRP7_5ASA)=c("Expression","Group")
NLRP7_5ASAAOV<-aov(Expression~Group,data=NLRP7_5ASA)
summary(NLRP7_5ASAAOV)
TukeyHSD(NLRP7_5ASAAOV)
NLRP7_5ASAplot=NLRP7_5ASA%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(NLRP7_5ASAplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("NLRP7")+ylim(0,10) +
  #geom_signif(annotations = c("**","***"), y_position = c(12,14),xmin = c(1, 1), xmax = c(2, 4)) +
  #geom_signif(annotations = c("**"), y_position = c(13),xmin = c(1), xmax = c(3)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))

IL6_5ASA=read.csv("IL6_5ASA.csv")
IL6_5ASA=stack(IL6_5ASA)
names(IL6_5ASA)=c("Expression","Group")
IL6_5ASAAOV<-aov(Expression~Group,data=IL6_5ASA)
summary(IL6_5ASAAOV)
TukeyHSD(IL6_5ASAAOV)
IL6_5ASAplot=IL6_5ASA%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(IL6_5ASAplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("IL6")+ylim(0,15) +
  geom_signif(annotations = c("***","***"), y_position = c(13,15),xmin = c(1, 1), xmax = c(2, 4)) +
  geom_signif(annotations = c("***"), y_position = c(14),xmin = c(1), xmax = c(3)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))

IL18R1_5ASA=read.csv("IL18R1_5ASA.csv")
IL18R1_5ASA=stack(IL18R1_5ASA)
names(IL18R1_5ASA)=c("Expression","Group")
IL18R1_5ASAplot=IL18R1_5ASA%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(IL18R1_5ASAplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("IL18R1")+ylim(0,10) +
  #geom_signif(annotations = c("***","***"), y_position = c(13,15),xmin = c(1, 1), xmax = c(2, 4)) +
  #geom_signif(annotations = c("***"), y_position = c(14),xmin = c(1), xmax = c(3)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))

IL1A_5ASA=read.csv("IL1A_5ASA.csv")
IL1A_5ASA=stack(IL1A_5ASA)
names(IL1A_5ASA)=c("Expression","Group")
IL1A_5ASAAOV<-aov(Expression~Group,data=IL1A_5ASA)
summary(IL1A_5ASAAOV)
TukeyHSD(IL1A_5ASAAOV)
IL1A_5ASAplot=IL1A_5ASA%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(IL1A_5ASAplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("IL1A")+ylim(0,12) +
  geom_signif(annotations = c("***","***"), y_position = c(10,11),xmin = c(1, 1), xmax = c(2, 4)) +
  geom_signif(annotations = c("***"), y_position = c(10.5),xmin = c(1), xmax = c(3)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))

ZBP1_5ASA=read.csv("ZBP1_5ASA.csv")
ZBP1_5ASA=stack(ZBP1_5ASA)
names(ZBP1_5ASA)=c("Expression","Group")
ZBP1_5ASAAOV<-aov(Expression~Group,data=ZBP1_5ASA)
summary(ZBP1_5ASAAOV)
TukeyHSD(ZBP1_5ASAAOV)
ZBP1_5ASAplot=ZBP1_5ASA%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(ZBP1_5ASAplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("ZBP1")+ylim(0,10) +
  #geom_signif(annotations = c("***","***"), y_position = c(10,11),xmin = c(1, 1), xmax = c(2, 4)) +
  #geom_signif(annotations = c("***"), y_position = c(10.5),xmin = c(1), xmax = c(3)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))

GZMB_5ASA=read.csv("GZMB_5ASA.csv")
GZMB_5ASA=stack(GZMB_5ASA)
names(GZMB_5ASA)=c("Expression","Group")
GZMB_5ASAAOV<-aov(Expression~Group,data=GZMB_5ASA)
summary(GZMB_5ASAAOV)
TukeyHSD(GZMB_5ASAAOV)
GZMB_5ASAplot=GZMB_5ASA%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(GZMB_5ASAplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("GZMB")+ylim(0,10) +
  #geom_signif(annotations = c("***","*"), y_position = c(10,11),xmin = c(1, 1), xmax = c(2, 4)) +
  geom_signif(annotations = c("*"), y_position = c(10),xmin = c(1), xmax = c(4)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))


TREM1_5ASA=read.csv("TREM1_5ASA.csv")
TREM1_5ASA=stack(TREM1_5ASA)
names(TREM1_5ASA)=c("Expression","Group")
TREM1_5ASAAOV<-aov(Expression~Group,data=TREM1_5ASA)
summary(TREM1_5ASAAOV)
TukeyHSD(TREM1_5ASAAOV)
TREM1_5ASAplot=TREM1_5ASA%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(TREM1_5ASAplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("TREM1")+ylim(0,10) +
  #geom_signif(annotations = c("***","*"), y_position = c(10,11),xmin = c(1, 1), xmax = c(2, 4)) +
  #geom_signif(annotations = c("*"), y_position = c(10),xmin = c(1), xmax = c(4)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))


#IFX VDZ
GSE73661_matrix=read.table("GSE73661_series_matrix.txt",header = T,sep="\t",row.names = 1,comment.char = "!")
GPL6244_anno=data.table::fread("GPL6244-17930.txt",skip = 12)
GSE73661_matrix$ID=rownames(GSE73661_matrix)
GSE73661=merge(x=GSE73661_matrix,y=GPL6244_anno,by="ID",all.x = F)
GSE73661$ID=NULL
rowMeans=apply(GSE73661,1,function(x)mean(as.numeric(x),na.rm=T))
GSE73661=GSE73661[order(rowMeans,decreasing = T),]
GSE73661=GSE73661[!duplicated(GSE73661[,dim(GSE73661)[2]]),]
GSE73661=na.omit(GSE73661)
rownames(GSE73661)=GSE73661[,dim(GSE73661)[2]]
GSE73661=GSE73661[,-dim(GSE73661)[2]]
GSE73661_Control=c(GSE73661$GSM1900159,GSE73661$1900225,
                   GSE73661$1900222,GSE73661$GSM1900150,
                   GSE73661$1900152,GSE73661$1900153,GSE73661$1900221,
                   GSE73661$1900224,GSE73661$GSM1900151,GSE73661$GSM1900223,
                   GSE73661$GSM1900226,GSE73661$GSM1900211)
GSE73661_Control=GSE73661[,c("GSM1900159","GSM1900225","GSM1900222","GSM1900150","GSM1900152",
                             "GSM1900153","GSM1900221","GSM1900224","GSM1900151","GSM1900223",
                             "GSM1900226","GSM1900211")]

GSE73661_IFX_qian=GSE73661[,c("GSM1900175","GSM1900208","GSM1900176","GSM1900202","GSM1900215","GSM1900195",
                              "GSM1900172","GSM1900185","GSM1900204","GSM1900148","GSM1900155","GSM1900210",
                              "GSM1900214","GSM1900154","GSM1900181","GSM1900158","GSM1900186","GSM1900180",
                              "GSM1900192","GSM1900184","GSM1900213","GSM1900217","GSM1900206")]
GSE73661_VDZ_qian=GSE73661[,c("GSM1900288","GSM1900227","GSM1900291","GSM1900299","GSM1900274","GSM1900230",
                              "GSM1900233","GSM1900317","GSM1900325","GSM1900236","GSM1900263","GSM1900318",
                              "GSM1900265","GSM1900323","GSM1900313","GSM1900285","GSM1900310","GSM1900307",
                              "GSM1900271","GSM1900314","GSM1900277","GSM1900280","GSM1900268","GSM1900319",
                              "GSM1900241","GSM1900297","GSM1900321","GSM1900283","GSM1900244","GSM1900305",
                              "GSM1900248","GSM1900320","GSM1900251","GSM1900324","GSM1900322","GSM1900260",
                              "GSM1900257","GSM1900315","GSM1900294","GSM1900302","GSM1900316")]
GSE73661_IFX_hou=GSE73661[,c("GSM1900177","GSM1900205","GSM1900212","GSM1900160","GSM1900190",
                             "GSM1900182","GSM1900193","GSM1900216")]
GSE73661_IFX_NR=GSE73661[,c("GSM1900178","GSM1900209","GSM1900179","GSM1900203","GSM1900218","GSM1900199",
                            "GSM1900188","GSM1900149","GSM1900157","GSM1900219","GSM1900156","GSM1900183",
                            "GSM1900187","GSM1900220","GSM1900207")]
write.csv(GSE73661_IFX_NR,"GSE73661_IFX_NR.csv")
write.csv(GSE73661_IFX_hou,"GSE73661_IFX_hou.csv")
write.csv(GSE73661_IFX_qian,"GSE73661_IFX_qian.csv")
write.csv(GSE73661_Control,"GSE73661_Control.csv")

GSE73661_VDZ_NR=GSE73661[,c("GSM1900289","GSM1900228","GSM1900292","GSM1900275","GSM1900231","GSM1900234",
                            "GSM1900286","GSM1900311","GSM1900308","GSM1900272","GSM1900278","GSM1900269",
                            "GSM1900242","GSM1900298","GSM1900284","GSM1900245","GSM1900249","GSM1900252",
                            "GSM1900258","GSM1900295","GSM1900303")]
GSE73661_VDZ_hou6=GSE73661[,c("GSM1900237","GSM1900266","GSM1900281","GSM1900164","GSM1900170","GSM1900163")]
GSE73661_VDZ_hou52=GSE73661[,c("GSM1900229","GSM1900198","GSM1900238","GSM1900267","GSM1900194","GSM1900273",
                               "GSM1900279","GSM1900282","GSM1900270","GSM1900201","GSM1900197")]
write.csv(GSE73661_VDZ_NR,"GSE73661_VCD_NR.csv")
write.csv(GSE73661_VDZ_hou6,"GSE73661_VCD_hou6.csv")
write.csv(GSE73661_VDZ_hou52,"GSE73661_VCD_hou52.csv")
write.csv(GSE73661_VDZ_qian,"GSE73661_VDZ_qian.csv")

#The following CSV data were obtained from the processed dataset and the corresponding gene expression data were extracted

AIM2_IFX=read.csv("AIM2_IFX.csv")
AIM2_IFX=stack(AIM2_IFX)
names(AIM2_IFX)=c("Expression","Group")
AIM2_IFXAOV<-aov(Expression~Group,data=AIM2_IFX)
summary(AIM2_IFXAOV)
TukeyHSD(AIM2_IFXAOV)
AIM2_IFXplot=AIM2_IFX%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(AIM2_IFXplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("AIM2")+ylim(0,10) +
  geom_signif(annotations = c("***","***"), y_position = c(8,9),xmin = c(1, 1), xmax = c(2, 4)) +
  geom_signif(annotations = c("**","*"), y_position = c(8,8),xmin = c(2,3), xmax = c(3,4)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))

IL6_IFX=read.csv("IL6_IFX.csv")
IL6_IFX=stack(IL6_IFX)
names(IL6_IFX)=c("Expression","Group")
IL6_IFXAOV<-aov(Expression~Group,data=IL6_IFX)
summary(IL6_IFXAOV)
TukeyHSD(IL6_IFXAOV)
IL6_IFXplot=IL6_IFX%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(IL6_IFXplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("IL6")+ylim(0,12) +
  geom_signif(annotations = c("***","**"), y_position = c(10,11),xmin = c(1, 1), xmax = c(2, 4)) +
  geom_signif(annotations = c("***","*"), y_position = c(10,10),xmin = c(2,3), xmax = c(3,4)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))

IL1B_IFX=read.csv("IL1B_IFX.csv")
IL1B_IFX=stack(IL1B_IFX)
names(IL1B_IFX)=c("Expression","Group")
IL1B_IFXAOV<-aov(Expression~Group,data=IL1B_IFX)
summary(IL1B_IFXAOV)
TukeyHSD(IL1B_IFXAOV)
IL1B_IFXplot=IL1B_IFX%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(IL1B_IFXplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("IL1B")+ylim(0,15) +
  geom_signif(annotations = c("***","***"), y_position = c(12,13),xmin = c(1, 1), xmax = c(2, 4)) +
  geom_signif(annotations = c("***","***","***"), y_position = c(12,12,14),xmin = c(2,3,1), xmax = c(3,4,3)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))


TNF_IFX=read.csv("TNF_IFX.csv")
TNF_IFX=stack(TNF_IFX)
names(TNF_IFX)=c("Expression","Group")
TNF_IFXAOV<-aov(Expression~Group,data=TNF_IFX)
summary(TNF_IFXAOV)
TukeyHSD(TNF_IFXAOV)
TNF_IFXplot=TNF_IFX%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(TNF_IFXplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("TNF")+ylim(0,10) +
  geom_signif(annotations = c("***","***"), y_position = c(8,9),xmin = c(1, 1), xmax = c(2, 4)) +
  geom_signif(annotations = c("***","**"), y_position = c(8,8),xmin = c(2,3), xmax = c(3,4)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))



TNF_IFX=read.csv("TNF_IFX.csv")
TNF_IFX=stack(TNF_IFX)
names(TNF_IFX)=c("Expression","Group")
TNF_IFXAOV<-aov(Expression~Group,data=TNF_IFX)
summary(TNF_IFXAOV)
TukeyHSD(TNF_IFXAOV)
TNF_IFXplot=TNF_IFX%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(TNF_IFXplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("TNF")+ylim(0,10) +
  geom_signif(annotations = c("***","***"), y_position = c(8,9),xmin = c(1, 1), xmax = c(2, 4)) +
  geom_signif(annotations = c("***","**"), y_position = c(8,8),xmin = c(2,3), xmax = c(3,4)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))



NLRP7_IFX=read.csv("NLRP7_IFX.csv")
NLRP7_IFX=stack(NLRP7_IFX)
names(NLRP7_IFX)=c("Expression","Group")
NLRP7_IFXAOV<-aov(Expression~Group,data=NLRP7_IFX)
summary(NLRP7_IFXAOV)
TukeyHSD(NLRP7_IFXAOV)
NLRP7_IFXplot=NLRP7_IFX%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(NLRP7_IFXplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("NLRP7")+ylim(0,8) +
  geom_signif(annotations = c("**","***"), y_position = c(6.5,7),xmin = c(1, 2), xmax = c(2, 3)) +
  geom_signif(annotations = c("**"), y_position = c(6.5),xmin = c(3), xmax = c(4)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))

IL1A_IFX=read.csv("IL1A_IFX.csv")
IL1A_IFX=stack(IL1A_IFX)
names(IL1A_IFX)=c("Expression","Group")
IL1A_IFXAOV<-aov(Expression~Group,data=IL1A_IFX)
summary(IL1A_IFXAOV)
TukeyHSD(IL1A_IFXAOV)
IL1A_IFXplot=IL1A_IFX%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(IL1A_IFXplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("IL1A")+ylim(0,12) +
  geom_signif(annotations = c("***","***"), y_position = c(10,11),xmin = c(1, 1), xmax = c(2, 4)) +
  geom_signif(annotations = c("***","**"), y_position = c(10,10),xmin = c(2,3), xmax = c(3,4)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))


IL18R1_IFX=read.csv("IL18R1_IFX.csv")
IL18R1_IFX=stack(IL18R1_IFX)
names(IL18R1_IFX)=c("Expression","Group")
IL18R1_IFXAOV<-aov(Expression~Group,data=IL18R1_IFX)
summary(IL18R1_IFXAOV)
TukeyHSD(IL18R1_IFXAOV)
IL18R1_IFXplot=IL18R1_IFX%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(IL18R1_IFXplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("IL18R1")+ylim(0,10) +
  geom_signif(annotations = c("***","***"), y_position = c(8,9),xmin = c(1, 1), xmax = c(2, 4)) +
  geom_signif(annotations = c("***","*"), y_position = c(8,8),xmin = c(2,3), xmax = c(3,4)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))

ZBP1_IFX=read.csv("ZBP1_IFX.csv")
ZBP1_IFX=stack(ZBP1_IFX)
names(ZBP1_IFX)=c("Expression","Group")
ZBP1_IFXAOV<-aov(Expression~Group,data=ZBP1_IFX)
summary(ZBP1_IFXAOV)
TukeyHSD(ZBP1_IFXAOV)
ZBP1_IFXplot=ZBP1_IFX%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(ZBP1_IFXplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("ZBP1")+ylim(0,10) +
  geom_signif(annotations = c("***","***"), y_position = c(9.2,10),xmin = c(1, 1), xmax = c(2, 4)) +
  geom_signif(annotations = c("***","***"), y_position = c(9.2,9.2),xmin = c(2,3), xmax = c(3,4)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))

GZMB_IFX=read.csv("GZMB_IFX.csv")
GZMB_IFX=stack(GZMB_IFX)
names(GZMB_IFX)=c("Expression","Group")
GZMB_IFXAOV<-aov(Expression~Group,data=GZMB_IFX)
summary(GZMB_IFXAOV)
TukeyHSD(GZMB_IFXAOV)
GZMB_IFXplot=GZMB_IFX%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(GZMB_IFXplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("GZMB")+ylim(0,10) +
  geom_signif(annotations = c("***","***"), y_position = c(8,9),xmin = c(1, 1), xmax = c(2, 4)) +
  geom_signif(annotations = c("***","**"), y_position = c(8,8),xmin = c(2,3), xmax = c(3,4)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))


TREM1_IFX=read.csv("TREM1_IFX.csv")
TREM1_IFX=stack(TREM1_IFX)
names(TREM1_IFX)=c("Expression","Group")
TREM1_IFXAOV<-aov(Expression~Group,data=TREM1_IFX)
summary(TREM1_IFXAOV)
TukeyHSD(TREM1_IFXAOV)
TREM1_IFXplot=TREM1_IFX%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(TREM1_IFXplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("TREM1")+ylim(0,10) +
  geom_signif(annotations = c("***","***"), y_position = c(9,10),xmin = c(1, 1), xmax = c(2, 4)) +
  geom_signif(annotations = c("***","**"), y_position = c(9,9),xmin = c(2,3), xmax = c(3,4)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))

AIM2_VDZ=read.csv("AIM2_VDZ.csv")
AIM2_VDZ=stack(AIM2_VDZ)
names(AIM2_VDZ)=c("Expression","Group")
AIM2_VDZAOV<-aov(Expression~Group,data=AIM2_VDZ)
summary(AIM2_VDZAOV)
TukeyHSD(AIM2_VDZAOV)
AIM2_VDZplot=AIM2_VDZ%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(AIM2_VDZplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("AIM2")+ylim(0,10) +
  geom_signif(annotations = c("***","***"), y_position = c(8,9),xmin = c(1, 1), xmax = c(2, 4)) +
  geom_signif(annotations = c("***","***"), y_position = c(8,8),xmin = c(2,3), xmax = c(3,4)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))

IL6_VDZ=read.csv("IL6_VDZ.csv")
IL6_VDZ=stack(IL6_VDZ)
names(IL6_VDZ)=c("Expression","Group")
IL6_VDZAOV<-aov(Expression~Group,data=IL6_VDZ)
summary(IL6_VDZAOV)
TukeyHSD(IL6_VDZAOV)
IL6_VDZplot=IL6_VDZ%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(IL6_VDZplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("IL6")+ylim(0,10) +
  geom_signif(annotations = c("**","*"), y_position = c(8.7,9.3),xmin = c(1, 1), xmax = c(2, 4)) +
  geom_signif(annotations = c("***","*"), y_position = c(8.7,8.7),xmin = c(2,3), xmax = c(3,4)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))

IL1B_VDZ=read.csv("IL1B_VDZ.csv")
IL1B_VDZ=stack(IL1B_VDZ)
names(IL1B_VDZ)=c("Expression","Group")
IL1B_VDZAOV<-aov(Expression~Group,data=IL1B_VDZ)
summary(IL1B_VDZAOV)
TukeyHSD(IL1B_VDZAOV)
IL1B_VDZplot=IL1B_VDZ%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(IL1B_IFXplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("IL1B")+ylim(0,15) +
  geom_signif(annotations = c("***","***"), y_position = c(12,13),xmin = c(1, 1), xmax = c(2, 4)) +
  geom_signif(annotations = c("***","***"), y_position = c(12,12),xmin = c(2,3), xmax = c(3,4)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))

TNF_VDZ=read.csv("TNF_VDZ.csv")
TNF_VDZ=stack(TNF_VDZ)
names(TNF_VDZ)=c("Expression","Group")
TNF_VDZAOV<-aov(Expression~Group,data=TNF_VDZ)
summary(TNF_VDZAOV)
TukeyHSD(TNF_VDZAOV)
TNF_VDZplot=TNF_VDZ%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(TNF_VDZplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("TNF")+ylim(0,10) +
  geom_signif(annotations = c("**","*"), y_position = c(8,9),xmin = c(1, 1), xmax = c(2, 4)) +
  geom_signif(annotations = c("**"), y_position = c(8),xmin = c(2), xmax = c(3)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))

NLRP7_VDZ=read.csv("NLRP7_VDZ.csv")
NLRP7_VDZ=stack(NLRP7_VDZ)
names(NLRP7_VDZ)=c("Expression","Group")
NLRP7_VDZAOV<-aov(Expression~Group,data=NLRP7_VDZ)
summary(NLRP7_VDZAOV)
TukeyHSD(NLRP7_VDZAOV)
NLRP7_VDZplot=NLRP7_VDZ%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(NLRP7_VDZplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("NLRP7")+ylim(0,8) +
  geom_signif(annotations = c("**","**"), y_position = c(6.5,7),xmin = c(1, 1), xmax = c(2, 4)) +
  geom_signif(annotations = c("***","**"), y_position = c(6.5,6.5),xmin = c(2,3), xmax = c(3,4)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))


IL1A_VDZ=read.csv("IL1A_VDZ.csv")
IL1A_VDZ=stack(IL1A_VDZ)
names(IL1A_VDZ)=c("Expression","Group")
IL1A_VDZAOV<-aov(Expression~Group,data=IL1A_VDZ)
summary(IL1A_VDZAOV)
TukeyHSD(IL1A_VDZAOV)
IL1A_VDZplot=IL1A_VDZ%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(IL1A_VDZplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("IL1A")+ylim(0,10) +
  geom_signif(annotations = c("***","***"), y_position = c(8.5,9.5),xmin = c(1, 1), xmax = c(2, 4)) +
  geom_signif(annotations = c("***","***"), y_position = c(8.5,8.5),xmin = c(2,3), xmax = c(3,4)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))

IL18R1_VDZ=read.csv("IL18R1_VDZ.csv")
IL18R1_VDZ=stack(IL18R1_VDZ)
names(IL18R1_VDZ)=c("Expression","Group")
IL18R1_VDZAOV<-aov(Expression~Group,data=IL18R1_VDZ)
summary(IL18R1_VDZAOV)
TukeyHSD(IL18R1_VDZAOV)
IL18R1_VDZplot=IL18R1_VDZ%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(IL18R1_VDZplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("IL18R1")+ylim(0,8) +
  geom_signif(annotations = c("***","***"), y_position = c(7,7.5),xmin = c(1, 1), xmax = c(2, 4)) +
  geom_signif(annotations = c("***","**"), y_position = c(7,7),xmin = c(2,3), xmax = c(3,4)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))

ZBP1_VDZ=read.csv("ZBP1_VDZ.csv")
ZBP1_VDZ=stack(ZBP1_VDZ)
names(ZBP1_VDZ)=c("Expression","Group")
ZBP1_VDZAOV<-aov(Expression~Group,data=ZBP1_VDZ)
summary(ZBP1_VDZAOV)
TukeyHSD(ZBP1_VDZAOV)
ZBP1_VDZplot=ZBP1_VDZ%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(ZBP1_VDZplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("ZBP1")+ylim(0,10) +
  geom_signif(annotations = c("***","***"), y_position = c(9.2,10),xmin = c(1, 1), xmax = c(2, 4)) +
  geom_signif(annotations = c("***","***"), y_position = c(9.2,9.2),xmin = c(2,3), xmax = c(3,4)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))


GZMB_VDZ=read.csv("GZMB_VDZ.csv")
GZMB_VDZ=stack(GZMB_VDZ)
names(GZMB_VDZ)=c("Expression","Group")
GZMB_VDZAOV<-aov(Expression~Group,data=GZMB_VDZ)
summary(GZMB_VDZAOV)
TukeyHSD(GZMB_VDZAOV)
GZMB_VDZplot=GZMB_VDZ%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(GZMB_VDZplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("GZMB")+ylim(0,10) +
  geom_signif(annotations = c("***","***"), y_position = c(8,9),xmin = c(1, 1), xmax = c(2, 4)) +
  geom_signif(annotations = c("***","***"), y_position = c(8,8),xmin = c(2,3), xmax = c(3,4)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))

TREM1_VDZ=read.csv("TREM1_VDZ.csv")
TREM1_VDZ=stack(TREM1_VDZ)
names(TREM1_VDZ)=c("Expression","Group")
TREM1_VDZAOV<-aov(Expression~Group,data=TREM1_VDZ)
summary(TREM1_VDZAOV)
TukeyHSD(TREM1_VDZAOV)
TREM1_VDZplot=TREM1_VDZ%>%group_by(Group)%>%summarise(n=n(),mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE))
ggplot(TREM1_VDZplot,aes(Group,mean,fill=Group))+geom_col() +
  geom_errorbar(aes(Group,ymin=mean-sd,ymax=mean+sd,color=Group),width=.5,size=1) +
  xlab("Experimental group")+ylab("TREM1")+ylim(0,10) +
  geom_signif(annotations = c("***","***"), y_position = c(9,9.5),xmin = c(1, 1), xmax = c(2, 4)) +
  geom_signif(annotations = c("***","***"), y_position = c(9,9),xmin = c(2,3), xmax = c(3,4)) +
  theme_bw() +
  scale_color_manual(values=c("black", "black", "black","black")) +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust = 0.5))

#singel-RNA data analysis
a1=read.table("GSM4949927_SingleCellCounts_UCI_1.txt",header=T,row.names = 1)
a2=read.table("GSM4949928_SingleCellCounts_UCI_2.txt",header=T,row.names = 1)
a3=read.table("GSM4949929_SingleCellCounts_UCI_3.txt",header=T,row.names = 1)
a4=read.table("GSM4949930_SingleCellCounts_UCI_4.txt",header=T,row.names = 1)
a5=read.table("GSM4949931_SingleCellCounts_UCI_5.txt",header=T,row.names = 1)
a6=read.table("GSM4949932_SingleCellCounts_UCI_6.txt",header=T,row.names = 1)
a7=read.table("GSM4949933_SingleCellCounts_UCI_7.txt",header=T,row.names = 1)
a8=read.table("GSM4949934_SingleCellCounts_UCI_8.txt",header=T,row.names = 1)
a9=read.table("GSM4949935_SingleCellCounts_UCI_9.txt",header=T,row.names = 1)
a10=read.table("GSM4949936_SingleCellCounts_UCI_10.txt",header=T,row.names = 1)
a11=read.table("GSM4949937_SingleCellCounts_UCI_11.txt",header=T,row.names = 1)
a2=merge(x=a1,y=a2,by="gene")
rownames(a2)=a2$gene
a1$gene=rownames(a1)
a2$gene=rownames(a2)
a3$gene=rownames(a3)
a4$gene=rownames(a4)
a5$gene=rownames(a5)
a6$gene=rownames(a6)
a7$gene=rownames(a7)
a8$gene=rownames(a9)
a9$gene=rownames(a9)
a10$gene=rownames(a10)
a11$gene=rownames(a11)
a2=merge(x=a3,y=a2,by="gene")
a2=merge(x=a4,y=a2,by="gene")
a2=merge(x=a5,y=a2,by="gene")
a2=merge(x=a6,y=a2,by="gene")
a2=merge(x=a7,y=a2,by="gene")
a2=merge(x=a8,y=a2,by="gene")
a2=merge(x=a9,y=a2,by="gene")
a2=merge(x=a10,y=a2,by="gene")
a2=merge(x=a11,y=a2,by="gene")

a2$gene=NULL
write.csv(a2,"a2.csv",quote = F)
a=read.csv("a2.csv",header=T,row.names = 1)
rownames(a)=a$gene
a$gene=NULL
pbmc=CreateSeuratObject(counts = a,min.cells = 100,min.features = 50)
pbmc[["precent.mt"]]=PercentageFeatureSet(pbmc,pattern="^MT-")
pbmc=subset(pbmc,subset=nFeature_RNA>300 & nFeature_RNA<6000 & precent.mt < 20)
pbmc=NormalizeData(pbmc,normalization.method = "LogNormalize",scale.factor = 10000)
pbmc=FindVariableFeatures(pbmc,selection.method = "vst",nfeatures = 3000)

all.genes=rownames(pbmc)
pbmc=ScaleData(pbmc,features = all.genes)
pbmc=RunPCA(pbmc,features = VariableFeatures(object=pbmc))
print(pbmc[["pca"]], dims = 1:30, nfeatures = 5)
pbmc=FindNeighbors(pbmc,dims=1:25)
pbmc=FindClusters(pbmc,resolution = 2)
pbmc=RunTSNE(object=pbmc,dims = 1:25)
DimPlot(pbmc,reduction = "tsne")

CM.mt<-cbind(as.data.frame(pbmc@meta.data$seurat_clusters),as.matrix(t(GetAssayData(pbmc,slot = "data"))))
#CM.mt=CM.mt[CM.mt$`pbmc@meta.data$seurat_clusters`%in%c(7),]
CM.mt=CM.mt[CM.mt$`pbmc@meta.data$seurat_clusters`%in%c(17,20,25,26,28),]
CM.mt=as.data.frame(CM.mt)
colnames(CM.mt)[1]<-"cluster"
CM.mt<-aggregate(.~cluster,CM.mt,mean)
rownames(CM.mt)<-CM.mt$cluster
CM.mt<-CM.mt[,-1]
CM.mt<-t(CM.mt)
hcl_result <- scHCL(scdata = CM.mt, numbers_plot = 6)
head(hcl_result$scHCL_probility,15)
write.csv(hcl_result$scHCL_probility,"7131516.csv")

hpca.se=celldex::BlueprintEncodeData()
table(hpca.se@colData@listData$label.fine)
cluster_pbmc=pbmc@meta.data$seurat_clusters
table(pbmc@meta.data$seurat_clusters)
datas=as.SingleCellExperiment(pbmc)
pred.hesc=SingleR(test=datas,ref=hpca.se,
                  labels=hpca.se$label.main,clusters=cluster_pbmc,
                  assay.type.test="logcounts",assay.type.ref="logcounts")
table(pred.hesc$labels)
celltype=data.frame(cluster=rownames(pred.hesc),celltype=pred.hesc$labels)
celltype$celltype[c(18,21,26,27,29)]="Macrophages"

for(i in 1:18375){
  index=pbmc@meta.data$seurat_clusters[i]
  pbmc@meta.data$celltype[i]=celltype[index,2]
}
DimPlot(pbmc,group.by = "celltype",reduction = "tsne",label = T,pt.size = 1)
pbmc.markers=FindAllMarkers(pbmc,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
write.table(pbmc.markers,file="pbmc.marker.txt",sep="\t",row.names = F,quote = F)
library(tidyverse)
sig.markers=pbmc.markers%>%select(gene,everything())%>%
  subset(p_val<0.05&abs(pbmc.markers$avg_log2FC)>1)
for(i in 1:nrow(sig.markers)){
  clusterID=as.numeric(sig.markers[i,]$cluster)
  celltypes=celltype[clusterID,2]
  sig.markers$celltype[i]=celltypes
}
sig.markerssave=sig.markers[,-8]
write.csv(sig.markerssave,"sig.makrer")
Pyroptosis=read.csv("焦亡.csv")
diffgene=sig.markers[sig.markers$gene%in%c("IL1B","IL1A","NLRP7","IL18R1","ZBP1","GZMB","TNF","IL6","TREM1","AIM2"),]

Pyroptosisgene=sig.markers[sig.markers$gene%in%Pyroptosis$Genes,]

mar=sig.markers[sig.markers$cluster%in%c(17,20,25,26,28),]
#Enrichment analysis performed by sangerbox

VlnPlot(object=pbmc,features = c("IL1B","IL1A","NLRP7","IL18R1","ZBP1"),
        group.by = "celltype")
FeaturePlot(object = pbmc,features = c("IL1B","IL1A","NLRP7","IL18R1","ZBP1"),
            cols = cividis(10))

FeaturePlot(object = pbmc,features = c("IL1B","IL1A"),blend = T)
VlnPlot(object=pbmc,features = c("GZMB","TNF","IL6","TREM1","AIM2"),
        group.by = "celltype")
FeaturePlot(object = pbmc,features = c("GZMB","TNF","IL6","TREM1","AIM2"),cols = cividis(10))
DotPlot(object=pbmc,features = c("IL1B","IL1A","NLRP7","IL18R1","ZBP1","GZMB","TNF","IL6","TREM1","AIM2"),
        group.by = "celltype") +
  coord_flip() +
  scale_color_viridis(option="H") +
  RotatedAxis()

DotPlot(object=pbmc,features = Pyroptosis$Genes,
        group.by = "celltype") +
  coord_flip() +
  scale_color_viridis(option="H") +
  RotatedAxis()

VlnPlot(pbmc,features =c("IL1B","IL1A","NLRP7","IL18R1","ZBP1","GZMB","TNF","IL6","TREM1","AIM2"),
        stack = T,pt.size = 0,direction="horizontal",x.lab="",y.lab="") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())
col=c("#E5D2DD","#53A85F","#F1BB72","#F3B1A0","#D6E7A3","#476D87","#E59CC4","#23452F","#585658","#9FA3A8","#CCC9E6")
markers=c("IL1B","IL1A","NLRP7","IL18R1","ZBP1","GZMB","TNF","IL6","TREM1","AIM2")
VlnPlot(pbmc,features =markers,
        stacked = T,pt.size = 0,direction="horizontal",x.lab="",y.lab="",group.by ="celltype") +
  theme(axis.text.x = element_blank())

pamc@meta.data$seurat_clusters
pamc=subset(pbmc,idents=c(17,20,25,26,28))
VlnPlot(pamc,features =markers,
        stacked = T,pt.size = 0,x.lab="",y.lab="",group.by ="seurat_clusters") +
  theme(axis.text.y = element_blank())
VlnPlot(pamc,features =markers,
        stack=T,pt.size = 0.1,direction="horizontal",x.lab="",y.lab="") +
  theme(axis.text.x = element_blank())


DotPlot(object=pamc,features = Pyroptosis$Genes) +
  coord_flip() +
  scale_color_viridis(option="F") +
  RotatedAxis()

DotPlot(object=pamc,features = c("IL1B","IL1A","NLRP7","IL18R1","ZBP1","GZMB","TNF","IL6","TREM1","AIM2")) +
  coord_flip() +
  scale_color_viridis(option="F") +
  RotatedAxis()


#AUCELL The first column of 焦亡 is the pyorptosis  pathway, and the second column is the related gene ID
cells_rankings <- AUCell_buildRankings(pbmc@assays$RNA@data)
h <- read.csv("焦亡.csv")
geneSets<-lapply(unique(h$term),function(x){h$Genes[h$term==x]})
names(geneSets) <- unique(h$term)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings,nCores = 1, aucMaxRank=nrow(cells_rankings)*0.05)
length(rownames(cells_AUC@assays@data$AUC))
geneSet <- "Pyroptosis_pathway"  
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])  #提取这个通路在每一个细胞的得分
pbmc$AUC <- aucs
df<- data.frame(pbmc@meta.data, pbmc@reductions$tsne@cell.embeddings)

class_avg <- df %>%
  group_by(celltype) %>%        #这里可以改成cluster  seurat_clusters/或者其他的annotation
  summarise(
    tSNE_1 = median(tSNE_1),
    tSNE_2 = median(tSNE_2)
  )
ggplot(df, aes(tSNE_1, tSNE_2))  +
  geom_point(aes(colour  = AUC)) + viridis::scale_color_viridis(option="H") +
  ggrepel::geom_label_repel(aes(label = celltype),
                            data = class_avg,
                            size = 6,
                            label.size = 0,
                            segment.color = NA
  )+   theme(legend.position = "none") + theme_bw()

cells_assignment <- AUCell_exploreThresholds(cells_AUC,plotHist=TRUE, assign=TRUE,thrP = 0.05)

#pseudotime analysis
Macrophages=subset(pbmc,idents=c(17,20,25,26,28))
data=as(as.matrix(Macrophages@assays$RNA@counts),"sparseMatrix")
pd=new("AnnotatedDataFrame",data=Macrophages@meta.data)
fdata=data.frame(gene_short_name=row.names(data),row.names = row.names(data))
fd=new("AnnotatedDataFrame",data=fdata)
HSMM=newCellDataSet(data,phenoData = pd,featureData = fd,
                    expressionFamily = negbinomial.size())

HSMM=estimateSizeFactors(HSMM)
HSMM=estimateDispersions(HSMM)
HSMM=detectGenes(HSMM,min_expr = 0.1)
disp_table=dispersionTable(HSMM)
unsup_clustering_gene=subset(disp_table,mean_expression>=0.1)
HSMM_myo=setOrderingFilter(HSMM,unsup_clustering_gene$gene_id)
plot_ordering_genes(HSMM)

HSMM_myo=reduceDimension(HSMM_myo,max_components = 2,reduction_method = "DDRTree",verbose = T)
HSMM_myo=orderCells(HSMM_myo)

plot_cell_trajectory(HSMM_myo,color_by = "Pseudotime",cell_size = 0.75)
plot_cell_trajectory(HSMM_myo,color_by = "seurat_clusters",cell_size = 0.75) +facet_wrap(~seurat_clusters,nrow = 1)
plot_cell_trajectory(HSMM_myo,color_by = "State",cell_size = 0.75) +facet_wrap(~State,nrow = 1)
plot_cell_trajectory(HSMM_myo,color_by = "State")
plot_cell_trajectory(HSMM_myo,color_by = "seurat_clusters",cell_size = 0.75)
plot_cell_trajectory(HSMM_myo,color_by = "celltype",cell_size = 0.75)

cg1=as.character(c("IL1B","IL1A","NLRP7","IL18R1","ZBP1"))
cg2=as.character(c("GZMB","TNF","IL6","TREM1","AIM2"))

plot_genes_in_pseudotime(HSMM_myo[cg2,],color_by = "seurat_clusters")
plot_genes_in_pseudotime(HSMM_myo[cg1,],color_by = "celltype")
plot_genes_in_pseudotime(HSMM_myo[cg2,],color_by = "celltype")

Pyroptosisnsx=Pyroptosis[Pyroptosis$Genes%in%pbmc.markers$gene,]
gc=as.character(Pyroptosisnsx$Genes)
diff_test_res=differentialGeneTest(HSMM_myo[gc,],fullModelFormulaStr = "~sm.ns(Pseudotime)")
plot_pseudotime_heatmap(HSMM_myo[rownames(diff_test_res),], cores = 1,show_rownames = T)

beam_res2=BEAM(HSMM_myo[gc,],branch_point = 2,cores = 4)
plot_genes_branched_heatmap(HSMM_myo[rownames(beam_res2),],branch_point = 2,cores = 4,show_rownames = T)

beam_res3=BEAM(HSMM_myo[gc,],branch_point = 3,cores = 4)
plot_genes_branched_heatmap(HSMM_myo[rownames(beam_res3),],branch_point = 3,cores = 4,show_rownames = T)

plot_genes_branched_pseudotime(HSMM_myo[cg1,],color_by = "celltype",branch_point = 1)
plot_genes_branched_pseudotime(HSMM_myo[cg2,],color_by = "celltype",branch_point = 2)


Pyroptosispbmc=Pyroptosis[Pyroptosis$Genes%in%rownames(pbmc),]
gc2=as.character(Pyroptosispbmc$Genes)
diff_test_resgc2=differentialGeneTest(HSMM_myo[gc2,],fullModelFormulaStr = "~sm.ns(Pseudotime)")
plot_pseudotime_heatmap(HSMM_myo[rownames(diff_test_resgc2),], cores = 1,show_rownames = T)

beam_resgc1=BEAM(HSMM_myo[gc2,],branch_point = 1,cores = 4)
plot_genes_branched_heatmap(HSMM_myo[rownames(beam_resgc1),],branch_point = 1,cores = 4,show_rownames = T)

beam_resgc2=BEAM(HSMM_myo[c(cg1,cg2),],branch_point = 1)
plot_genes_branched_pseudotime(HSMM_myo[rownames(beam_resgc2),],color_by = "seurat_clusters",branch_point = 1)


