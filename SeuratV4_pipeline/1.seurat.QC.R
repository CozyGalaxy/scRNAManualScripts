## Rcode description: scRNA-seq QC by Seurat
## Author: Zhu Yuanzhen
##
## Modified by: Zhang Pei (zhangpei@genomics.cn), adding cell cycle regression and marker plot.
##
## Date: 2022.5.17

parser = argparse::ArgumentParser(description="scRNA-seq QC by Seurat")
parser$add_argument('-I','--input', help = 'input raw matrix or 10X-like dir')
parser$add_argument('-G','--mregex', help = 'regex for mitochondria')
parser$add_argument('-X','--mini', help = 'minimum nFeature_RNA cut off')
parser$add_argument('-Y','--maxi', help = 'maximum nFeature_RNA cut off')
parser$add_argument('-Z','--mt.percent', help = 'default 5, should be 1 for snRNA-seq')
parser$add_argument('-T','--transform', help = 'transform method, log (default) or sct')
parser$add_argument('-F','--nf', help = 'HVG number set')
parser$add_argument('-D','--dim', help = 'PCA dim')
parser$add_argument('-U','--dusage', help = 'PCA dim usage')
parser$add_argument('-P','--percentage', help = 'doublets percentage')
parser$add_argument('-R','--res', help = 'Map resolution usage')
parser$add_argument('-O','--out', help = 'out directory')
parser$add_argument("-M",'--marker', help= 'marker list')
parser$add_argument("-C",'--cc_gene', help= 'input list for cell cycle genes, 1 set for seurat default cc genes')
parser$add_argument('-S','--sample', help = 'sample name')
args = parser$parse_args()

### Load library
library(Seurat)
library(DoubletFinder)
library(ggplot2)
library(dplyr)
library(patchwork)

#EC.data <- readRDS(args$input)
EC.data <- Read10X(data.dir = args$input, gene.column = 1)
args$mregex <- if(!is.null(args$mregex)) args$mregex else "MT-"
args$mregex <- paste("^",args$mregex, sep="")
args$mt.percent <- as.numeric(if(!is.null(args$mt.percent)) args$mt.percent else 5)
nf.usage <- as.numeric(if(!is.null(args$nf)) args$nf else 3000)
dim.all <- as.numeric(if(!is.null(args$dim)) args$dim else 50)
nf_mini <- as.numeric(if(!is.null(args$mini)) args$mini else 400)
nf_maxi <- as.numeric(if(!is.null(args$maxi)) args$maxi else 10000)
dim.usage <- as.numeric(if(!is.null(args$dusage)) args$dusage else round(0.5*dim.all))
doublets.percentage <- if(!is.null(args$percentage)) args$percentage else 0.05
doublets.percentage <- as.numeric(doublets.percentage)
res.usage <- as.numeric(if(!is.null(args$res)) args$res else 0.6)
args$transform <- if(!is.null(args$transform)) args$transform else "log"

### Creat Seurat object, basic filtering and statistics
EC <- CreateSeuratObject(EC.data, project = args$sample, min.cells = 3, min.features = 200)
EC[["percent.mt"]] <- PercentageFeatureSet(EC, pattern = args$mregex)

pdf(paste0(args$out, "/count_mt.cor.pdf"))
FeatureScatter(EC, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()
pdf(paste0(args$out, "/count_gene.cor.pdf"))
FeatureScatter(EC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

p1 <- VlnPlot(EC, features = "nFeature_RNA",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#1B9E77") + NoLegend() + xlab("Gene") + labs(title="")+ theme(axis.text.x = element_blank())
p2 <- VlnPlot(EC, features = "nCount_RNA",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#D95F02") + NoLegend() + xlab("Transcript") + labs(title="")+ theme(axis.text.x = element_blank()) #+ ylim(0,15000)
p3 <- VlnPlot(EC, features = "percent.mt",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#7570B3") + NoLegend() + xlab("percent.mt") + labs(title="")+ theme(axis.text.x = element_blank())
p <-p1|p2|p3
pdf(paste0(args$out,"/QC_before.dis.pdf"))
print(p)
dev.off()

EC <- subset(EC, subset = nFeature_RNA > nf_mini & nFeature_RNA < nf_maxi & percent.mt < args$mt.percent)
p1 <- VlnPlot(EC, features = "nFeature_RNA",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#1B9E77") + NoLegend() + xlab("Gene") + labs(title="")+ theme(axis.text.x = element_blank())
p2 <- VlnPlot(EC, features = "nCount_RNA",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#D95F02") + NoLegend() + xlab("Transcript") + labs(title="")+ theme(axis.text.x = element_blank()) #+ ylim(0,15000)
p3 <- VlnPlot(EC, features = "percent.mt",pt.size = 0) + geom_boxplot(width=.3,col="black",fill="white",outlier.colour=NA) + scale_fill_manual(values = "#7570B3") + NoLegend() + xlab("percent.mt") + labs(title="")+ theme(axis.text.x = element_blank())
p <-p1|p2|p3
pdf(paste0(args$out,"/QC_after.dis.pdf"))
print(p)
dev.off()

if(args$transform == "sct"){
	library(sctransform)
	EC <- PercentageFeatureSet(EC, pattern = args$mregex, col.name = "percent.mt")
	EC <- SCTransform(EC, vars.to.regress = "percent.mt", verbose = FALSE)
}else{
	EC <- NormalizeData(EC, normalization.method = "LogNormalize", scale.factor = 10000)
	EC <- FindVariableFeatures(EC, selection.method = "vst", nfeatures = nf.usage)
	EC <- ScaleData(EC)
}

pdf(paste0(args$out, "/VariableFeaturePlot.pdf"))
VariableFeaturePlot(EC)
dev.off()

### PCA, statistics
EC <- RunPCA(EC, npcs = dim.all)
pdf(paste0(args$out,"/PCA_DimHeatmap.pdf"), height = 20, width = 20)
DimHeatmap(EC, dims = 1:dim.usage, cells = 100, balanced = TRUE)
dev.off()
if(args$transform != "sct"){
	EC <- JackStraw(EC, dims=dim.usage, num.replicate = 100)
	EC <- ScoreJackStraw(EC, dims = 1:dim.usage)
	pdf(paste0(args$out, "/PCA_JackStrawPlot.pdf"))
	JackStrawPlot(EC, dims = 1:20)
	dev.off()
}
pdf(paste0(args$out, "/PCA_ElbowPlot.pdf"))
ElbowPlot(EC)
dev.off()


### Find neighbors, UMAP and T-SNE
EC <- FindNeighbors(EC, dims = 1:dim.usage)
EC <- FindClusters(EC, resolution = res.usage)
EC <- RunUMAP(EC, dims = 1:dim.usage)
EC <- RunTSNE(EC, dims = 1:dim.usage)

p1 <- DimPlot(EC, reduction = "umap")
p2 <- DimPlot(EC, reduction = "tsne")
p <- p1|p2
ggsave(filename = paste0(args$out, "/Umap_TSNE.orig.pdf"), plot = p, device = "pdf", width = 15, height = 7)

## Define Find_doublet function
# the doublet may be two "intermediate" cell state, and find_doublet may mistakenly remove them !
Find_doublet <- function(data){
	if(args$transform == "sct"){
		sweep.res.list <- paramSweep_v3(data, PCs = 1:dim.usage, sct = TRUE)
	}else{
		sweep.res.list <- paramSweep_v3(data, PCs = 1:dim.usage, sct = FALSE)
	}
	sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
	bcmvn <- find.pK(sweep.stats)
	nExp_poi <- round(as.numeric(doublets.percentage)*ncol(data))
	p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
	if(args$transform == "sct"){
		data <- doubletFinder_v3(data, PCs = 1:dim.usage, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
	}else{
		data <- doubletFinder_v3(data, PCs = 1:dim.usage, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
	}
	colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
	#data<-subset(data,subset=doublet_info=="Singlet")
	data
}


#Find doublets
EC <- Find_doublet(EC)
write.table(EC@meta.data,paste0(args$out,"/",args$sample,"_doublets_info.txt"),sep="\t", quote=FALSE)
EC <- subset(EC,subset=doublet_info=="Singlet")
EC@meta.data$split = args$sample

if(args$transform == "sct"){
	EC <- PercentageFeatureSet(EC, pattern = args$mregex, col.name = "percent.mt")
	EC <- SCTransform(EC, vars.to.regress = "percent.mt", verbose = FALSE)
}else{
	EC <- FindVariableFeatures(EC, selection.method = "vst", nfeatures = nf.usage)
	EC <- ScaleData(EC)
}

if(!is.null(args$cc_gene)){
	if(args$cc_gene == 1){
		s.genes <- cc.genes$s.genes
		g2m.genes <- cc.genes$g2m.genes
	}else{
		manual_cc_genes <- read.table(args$cc_gene,header=F)
		s.genes <- manual_cc_genes[manual_cc_genes[2] == "s.genes",1]
		g2m.genes <- manual_cc_genes[manual_cc_genes[2] == "g2m.genes",1]
	}
	EC <- CellCycleScoring(EC, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
	if(args$transform == "sct"){
		library(glmGamPoi)
		EC <- SCTransform(EC, method = "glmGamPoi", vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
	}else{
		EC <- ScaleData(EC, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(EC))
	}
}

EC <- RunPCA(EC, features = VariableFeatures(EC), npcs = dim.all)

#pct <- EC[["pca"]]@stdev / sum(EC[["pca"]]@stdev) * 100
#cumu <- cumsum(pct)
#co1 <- which(cumu > 80 & pct < 5)[1]
#co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
#pcs <- min(co1, co2)


EC <- FindNeighbors(EC, dims = 1:dim.usage)
EC <- FindClusters(EC, resolution = res.usage)
EC <- RunUMAP(EC, dims = 1:dim.usage)
EC <- RunTSNE(EC, dims = 1:dim.usage)

p1 <- DimPlot(EC, reduction = "umap", label=T)
p2 <- DimPlot(EC, reduction = "tsne", label=T)
p <- p1|p2
ggsave(filename = paste0(args$out, "/Umap_TSNE_new.pdf"), plot = p, device = "pdf", width = 15, height = 7)

nfeature_quantile <- quantile(EC$nFeature_RNA, probs = c(.2, .8))
p1 <- FeaturePlot(EC, "nFeature_RNA", cols=c("blue","yellow","red"), min.cutoff=nfeature_quantile[[1]], max.cutoff=nfeature_quantile[[2]], reduction = "umap")
p2 <- FeaturePlot(EC, "nFeature_RNA", cols=c("blue","yellow","red"), min.cutoff=nfeature_quantile[[1]], max.cutoff=nfeature_quantile[[2]], reduction = "tsne")
p <- p1|p2
ggsave(filename = paste0(args$out, "/Umap_TSNE_new_nfeature.pdf"), plot = p, device = "pdf", width = 15, height = 7)


plot_theme<-theme(panel.background = element_blank(),axis.line = element_line(size=0.1),axis.ticks = element_line(size=0.1),axis.text = element_text(size=7),axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),axis.ticks.length = unit(1, "pt"),plot.title = element_text(size = 7, face = "bold",margin=margin(0,0,4,0)),plot.margin = unit(c(0,0,0,0),"cm"),legend.key.size = unit(0.1, 'cm'),legend.key.height = unit(0.1, 'cm'),legend.key.width = unit(0.1, 'cm'),legend.title = element_text(size=6),legend.text = element_text(size=6))

if(!is.null(args$marker)){
	sample_marker <- read.table(args$marker,header=F)
	colnames(sample_marker) = c("marker", "anno")
	dot_p <- DotPlot(EC, features=sample_marker$marker,cols = c("grey","red"),col.min= -0.5, col.max =2,dot.scale=3)+labs(x="",y="")+plot_theme
	ggsave(filename = paste0(args$out, "/Marker.DotPlot.pdf"), plot = dot_p, device = "pdf", width = round(nrow(sample_marker)/8+2), height = 4)
}

write.table(EC@meta.data,paste0(args$out,"/",args$sample,".cell_info.txt"),sep="\t", quote=FALSE)
markers <- FindAllMarkers(EC, only.pos=TRUE)
write.table(markers,paste0(args$out,"/",args$sample,".marker_info.txt"),sep="\t", quote=FALSE)

saveRDS(EC,paste(args$out,"/",args$sample,"_QC.RDS",sep=""))

