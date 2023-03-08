## Rcode description: scRNA-seq integration by Seurat
##
## Author: Zhang Pei (zhangpei@genomics.cn)
##
## Date: 2022.5.17

parser = argparse::ArgumentParser(description="scRNA-seq integration by Seurat")
parser$add_argument('-I','--input', help = 'seurat object name and directory, without header')
parser$add_argument('-G','--mregex', help = 'regex for mitochondria')
parser$add_argument('-T','--transform', help = 'transform method, LogNormalize (default) or SCT')
parser$add_argument('-F','--nf',help = 'number of variable features')
parser$add_argument('-D','--dim',help = 'PCA dim')
parser$add_argument('-U','--dusage',help = 'PCA dim usage')
parser$add_argument('-W','--workflow', help = 'workflow for integration, including cca (default), rpca and rlsi')
parser$add_argument('-A','--kanchor', help = 'number of k.anchor for integration, 5 as default for cca and 15 as default for rpca')
parser$add_argument('-B','--ref_based', help = 'reference based, NULL as default, should be ranks in input file and separated by comma, such as 1,2,3')
parser$add_argument('-R','--res',help = 'Map resolution usage, 1 as default')
parser$add_argument('-O','--out',help = 'out directory, ./ as default')
parser$add_argument("-M", '--marker',help= 'chosen marker list, NULL as default, file containing chosen marker genes for DotPlot, without header')
parser$add_argument("-C",'--cc_gene',help= 'input list for cell cycle genes, NULL as default, 1 for cc genes provided by seurat, or should be file including gene name and type (s.genes and g2m.genes) and without header')
parser$add_argument('-S','--sample',help = 'Sample name, also be used for slot name. There must be no replicate names for the total project')
parser$add_argument('-P','--processor',help = 'processor number, should not be larger than 5' )
parser$add_argument('-E','--step',help = 'step for integration, 1 for replicates of every sample, 2 for samples of the dataset, default 1')
parser$add_argument('--algorithm',help='clustering algorithm, leiden (default) or louvain')

args = parser$parse_args()
args$mregex <- if(!is.null(args$mregex)) args$mregex else "Mitochondria"
args$mregex <- paste("^",args$mregex, sep="")
nf.usage <- as.numeric(if(!is.null(args$nf)) args$nf else 4000)
dim.all <- as.numeric(if(!is.null(args$dim)) args$dim else 50)
dim.usage <- as.numeric(if(!is.null(args$dusage)) args$dusage else round(0.75*dim.all))
res.usage <- as.numeric(if(!is.null(args$res)) args$res else 1)
args$out <- if(!is.null(args$out)) args$out else "./"
args$step <- as.numeric(if(!is.null(args$step)) args$step else 1)
args$workflow <- if(!is.null(args$workflow)) args$workflow else "cca"
args$kanchor <- as.numeric(if(!is.null(args$kanchor)) args$kanchor else 5)
args$transform <- if(!is.null(args$transform)) args$transform else "LogNormalize"
if(args$workflow == "rpca"){args$kanchor <- 20}

args$algorithm <- if(!is.null(args$algorithm)) args$algorithm else "leiden"


library(Seurat)
library(DoubletFinder)
library(ggplot2)
library(dplyr)
library(patchwork)
library(plyr)

if(args$transform == "SCT"){
	library(sctransform)
	library(glmGamPoi)
}

if(args$algorithm == "leiden"){
    library(reticulate)
    use_python("/ldfssz1/ST_DIVERSITY/PUB/USER/zhangpei/bin/python3.7.13/bin/python3", required = T)
    py_config()
    library(igraph)
    library(leiden)
}


args$processor <- as.numeric(if(!is.null(args$processor)) args$processor else 1)
if(args$processor > 1){
	library(future)
	plan("multicore", workers = args$processor)
	options(future.globals.maxSize = args$processor * 1500 * 1000 * 1024^2)
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
}


data.list <- read.table(args$input, header=F)
colnames(data.list) = c("sample", "dir")
all_rds <- tapply(data.list$dir,data.list$sample, readRDS)

if(args$step == 1){
	if(args$transform == "SCT"){
		all_rds <- lapply(X = all_rds, FUN = function(x){
			x <- PercentageFeatureSet(x, pattern = args$mregex, col.name = "percent.mt")
			x <- SCTransform(x, method = "glmGamPoi", vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA'), verbose = FALSE)
			if(!is.null(args$cc_gene)){
				x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, assay = 'SCT', set.ident = TRUE)
				x <- SCTransform(x, method = "glmGamPoi", vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'), verbose = FALSE)
			}
			x <- RunPCA(x, npcs = dim.all)
		})
	}else{
		all_rds <- lapply(X = all_rds, FUN = function(x) {
			x <- PercentageFeatureSet(x, pattern = args$mregex, col.name = "percent.mt")
			x <- NormalizeData(x)
			x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nf.usage)
			x <- ScaleData(x, vars.to.regress = c("percent.mt", "nCount_RNA", "nFeature_RNA"))
			x <- RunPCA(x, npcs = dim.all)
		})
	}
}
features <- SelectIntegrationFeatures(object.list = all_rds, nfeatures = nf.usage, fvf.nfeatures = nf.usage)
if(args$transform == "SCT"){
	all_rds <- PrepSCTIntegration(object.list = all_rds, anchor.features = features)
#	all_rds <- lapply(X = all_rds, FUN = RunPCA, features = features)
}


if(is.null(args$ref_based)){
	anchors <- FindIntegrationAnchors(object.list = all_rds, anchor.features = features, normalization.method = args$transform, dims=1:dim.usage, reduction = args$workflow, k.anchor = args$kanchor)
}else{
	args$ref_based <- lapply(strsplit(args$ref_based, ','), as.numeric)[[1]]
	anchors <- FindIntegrationAnchors(object.list = all_rds, reference = c(args$ref_based), anchor.features = features, normalization.method = args$transform, dims=1:dim.usage, reduction = args$workflow, k.anchor = args$kanchor)
}
assay_name <- paste(args$sample, "integrated", sep="_")
data_integrated <- IntegrateData(anchorset = anchors, normalization.method = args$transform, dims = 1:dim.usage, k.weight = 50, new.assay.name = assay_name)
DefaultAssay(data_integrated) <- assay_name

if(args$transform != "SCT"){
	if(!is.null(args$cc_gene)){
		data_integrated <- CellCycleScoring(data_integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
		data_integrated <- ScaleData(data_integrated, vars.to.regress = c('S.Score', 'G2M.Score'), features = rownames(data_integrated))
	}else{
		data_integrated <- ScaleData(data_integrated, features = rownames(data_integrated))
	}
}


data_integrated <- RunPCA(data_integrated, npcs = dim.all, verbose = FALSE)
dim.usage2 = dim.all - 5
pdf(paste0(args$out,"/PCA_DimHeatmap.pdf"), height = 20, width = 20)
DimHeatmap(data_integrated, dims = 1:dim.usage, cells = 100, balanced = TRUE)
dev.off()
#if(args$transform != "SCT"){
#	data_integrated <- JackStraw(data_integrated, dims=dim.usage2, num.replicate = 100)
#	data_integrated <- ScoreJackStraw(data_integrated, dims = 1:dim.usage)
#	pdf(paste0(args$out, "/PCA_JackStrawPlot.pdf"))
#	JackStrawPlot(data_integrated, dims = 1:dim.usage)
#	dev.off()
#}
pdf(paste0(args$out, "/PCA_ElbowPlot.pdf"))
ElbowPlot(data_integrated)
dev.off()

data_integrated <- RunTSNE(data_integrated, reduction = "pca", dims = 1:dim.usage)
data_integrated <- RunUMAP(data_integrated, reduction = "pca", dims = 1:dim.usage)
data_integrated <- FindNeighbors(data_integrated, reduction = "pca", dims = 1:dim.usage)

if(args$algorithm == "leiden"){
    data_integrated <- FindClusters(data_integrated, resolution = res.usage, algorithm=4, method='igraph')
}else{
    data_integrated <- FindClusters(data_integrated, resolution = res.usage)
}


p1 <- DimPlot(data_integrated, reduction = "umap", label=T, raster=F)
p2 <- DimPlot(data_integrated, reduction = "tsne", label=T, raster=F)
p <- p1|p2
ggsave(filename = paste0(args$out, "/Umap_TSNE.pdf"), plot = p, device = "pdf", width = 16, height = 7)

nfeature_quantile <- quantile(data_integrated$nFeature_RNA, probs = c(.2, .8))
p1 <- FeaturePlot(data_integrated, "nFeature_RNA", cols=c("blue","yellow","red"), min.cutoff=nfeature_quantile[[1]], max.cutoff=nfeature_quantile[[2]], reduction = "umap")
p2 <- FeaturePlot(data_integrated, "nFeature_RNA", cols=c("blue","yellow","red"), min.cutoff=nfeature_quantile[[1]], max.cutoff=nfeature_quantile[[2]], reduction = "tsne")
p <- p1|p2
ggsave(filename = paste0(args$out, "/Umap_TSNE.nfeature.pdf"), plot = p, device = "pdf", width = 15, height = 7)

p1 <- DimPlot(data_integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(data_integrated, reduction = "tsne", group.by = "orig.ident")
p <- p1|p2
ggsave(filename = paste0(args$out, "/Umap_TSNE.sample.pdf"), plot = p, device = "pdf", width = 15, height = 7)


if(!is.null(args$marker)){
	sample_marker <- read.table(args$marker,header=F)
	colnames(sample_marker) = c("marker", "anno")
	plot_theme<-theme(panel.background = element_blank(),axis.line = element_line(size=0.1),axis.ticks = element_line(size=0.1),axis.text = element_text(size=7),axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),axis.ticks.length = unit(1, "pt"),plot.title = element_text(size = 7, face = "bold",margin=margin(0,0,4,0)),plot.margin = unit(c(0,0,0,0),"cm"),legend.key.size = unit(0.1, 'cm'),legend.key.height = unit(0.1, 'cm'),legend.key.width = unit(0.1, 'cm'),legend.title = element_text(size=6),legend.text = element_text(size=6))

	dot_p <- DotPlot(data_integrated, assay="RNA", features=sample_marker$marker,cols = c("grey","red"),col.min= -0.5,dot.scale=3)+labs(x="",y="")+plot_theme
	ggsave(filename = paste0(args$out, "/Marker.DotPlot.pdf"), plot = dot_p, device = "pdf", width = round(nrow(sample_marker)/8+2), height = 4)
}

write.table(data_integrated@meta.data,paste0(args$out,"/",args$sample,".cell_info.txt"),sep="\t", quote=FALSE)
saveRDS(data_integrated,paste(args$out,"/",args$sample,".integrated.RDS",sep=""))

markers <- FindAllMarkers(data_integrated, only.pos=TRUE)
write.table(markers,paste0(args$out,"/",args$sample,".marker_info.txt"),sep="\t", quote=FALSE)

library(clustree)
data_integrated <- FindClusters(data_integrated, resolution = seq(0.2,1.5,0.1), n.start=10)
pdf(paste0(args$out,"/",args$sample,".clustree.pdf"), width = 15, height = 15)
assay_name <- paste(assay_name, "snn", "res", sep="_")
assay_name <- paste(assay_name, ".", sep="")
clustree(data_integrated, prefix = assay_name)

dev.off()
