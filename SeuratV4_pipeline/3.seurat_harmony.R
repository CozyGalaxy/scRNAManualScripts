## Rcode description: scRNA-seq integration by Seurat using harmony
##
## Author: Zhang Pei (zhangpei@genomics.cn)
##
## Date: 2022.7.13

parser = argparse::ArgumentParser(description="scRNA-seq integration by Seurat using harmony")
parser$add_argument('-I','--input', help = 'seurat object name and directory, without header')
parser$add_argument('-G','--mregex', help = 'regex for mitochondria')
parser$add_argument('-T','--transform', help = 'transform method, LogNormalize or SCT (default)')
parser$add_argument('-F','--nf',help = 'number of variable features')
parser$add_argument('-D','--dim',help = 'PCA dim')
parser$add_argument('-X','--dusage1',help = 'PCA dim usage for integration')
parser$add_argument('-Y','--dusage2',help = 'PCA dim usage for clustering, default=20')
parser$add_argument('-R','--res',help = 'Map resolution usage, 1 as default')
parser$add_argument('-O','--output',help = 'out directory, ./ as default')
parser$add_argument('-M','--marker',help= 'chosen marker list, NULL as default, file containing chosen marker genes for DotPlot, without header')
parser$add_argument('-C','--cc_gene',help= 'input list for cell cycle genes, NULL as default, 1 for cc genes provided by seurat, or should be file including gene name and type (s.genes and g2m.genes) and without header')
parser$add_argument('-S','--sample',help = 'Sample name, also be used for slot name. There must be no replicate names for the total project')

parser$add_argument('--algorithm',help='clustering algorithm, leiden (default) or louvain')
parser$add_argument('--python',help='python3 bin path')

##harmony arguments
parser$add_argument('--theta',help = 'Diversity clustering penalty parameter. Specify for each variable in group.by.vars. Default theta=2. theta=0 does not encourage any diversity. Larger values of theta result in more diverse clusters')
parser$add_argument('--lambda',help= 'Ridge regression penalty parameter. Specify for each variable in group.by.vars. Default lambda=1. Lambda must be strictly positive. Smaller values result in more aggressive correction.')
parser$add_argument('--sigma',help='Width of soft kmeans clusters. Default sigma=0.1. Sigma scales the distance from a cell to cluster centroids. Larger values of sigma result in cells assigned to more clusters. Smaller values of sigma make soft kmeans cluster approach hard clustering.')
parser$add_argument('--nclust',help='Number of clusters in model. nclust=1 equivalent to simple linear regression.')
parser$add_argument('--blocksize',help='What proportion of cells to update during clustering. Between 0 to 1, default 0.01. Larger values may be faster but less accurate')



args = parser$parse_args()
args$mregex <- if(!is.null(args$mregex)) args$mregex else "Mitochondria"
args$mregex <- paste("^",args$mregex, sep="")
nf.usage <- as.numeric(if(!is.null(args$nf)) args$nf else 4000)
dim.all <- as.numeric(if(!is.null(args$dim)) args$dim else 50)
dim.usage1 <- as.numeric(if(!is.null(args$dusage1)) args$dusage1 else round(0.9*dim.all))
dim.usage2 <- as.numeric(if(!is.null(args$dusage2)) args$dusage2 else round(0.6*dim.all))
res.usage <- as.numeric(if(!is.null(args$res)) args$res else 1.5)
args$output <- if(!is.null(args$output)) args$output else "./"
args$transform <- if(!is.null(args$transform)) args$transform else "SCT"

args$algorithm <- if(!is.null(args$algorithm)) args$algorithm else "leiden"
args$python <- if(!is.null(args$python)) args$python else "/ldfssz1/ST_DIVERSITY/PUB/USER/zhangpei/bin/python3.7.13/bin/python3"


args$theta <- as.numeric(if(!is.null(args$theta)) args$theta else 2)
args$lambda <- as.numeric(if(!is.null(args$lambda)) args$lambda else 1)
args$sigma <- as.numeric(if(!is.null(args$sigma)) args$sigma else 0.1)
args$nclust <- if(!is.null(args$nclust)) as.numeric(args$nclust) else NULL
args$blocksize <- as.numeric(if(!is.null(args$blocksize)) args$blocksize else 0.01)


library(harmony)
library(Seurat)
library(ggplot2)
library(plyr)

if(args$transform == "SCT"){
	library(sctransform)
	library(glmGamPoi)
}

if(args$algorithm == "leiden"){
	library(reticulate)
	use_python(args$python, required = T)
	py_config()
	library(igraph)
	library(leiden)
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

if(args$transform == "SCT"){
	all_rds <- lapply(X = all_rds, FUN = function(x){
		DefaultAssay(x) <- "RNA"
		for (i in Assays(x)){ if(i != "RNA") {x[[i]] <- NULL}}
		x <- PercentageFeatureSet(x, pattern = args$mregex, col.name = "percent.mt")
		if(!is.null(args$cc_gene)){
			x <- SCTransform(x, method = "glmGamPoi", vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA'), verbose = FALSE)
			x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, assay = 'SCT', set.ident = TRUE)
			x <- SCTransform(x, method = "glmGamPoi", vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'), verbose = FALSE)
		}else{
			x <- SCTransform(x, method = "glmGamPoi", vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA'), verbose = FALSE)
		}
        x <- RunPCA(x, npcs = dim.all)
	})
}else{
	all_rds <- lapply(X = all_rds, FUN = function(x) {
		DefaultAssay(x) <- "RNA"
		for (i in Assays(x)){ if(i != "RNA") {x[[i]] <- NULL}}
		x <- PercentageFeatureSet(x, pattern = args$mregex, col.name = "percent.mt")
        x <- NormalizeData(x)
		if(!is.null(args$cc_gene)){
			x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
		}
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nf.usage)
		if(!is.null(args$cc_gene)){
			x <- ScaleData(x, vars.to.regress = c("percent.mt", "nCount_RNA", "nFeature_RNA",'S.Score', 'G2M.Score'))
		}else{
        	x <- ScaleData(x, vars.to.regress = c("percent.mt", "nCount_RNA", "nFeature_RNA"))
		}
        x <- RunPCA(x, npcs = dim.all)
	})
}

for(i in 1:length(all_rds)) all_rds[[i]]@meta.data$temp_sample <- names(all_rds[i])
data_integrated <- merge(x=all_rds[[1]], y=all_rds[2:length(all_rds)], merge.data = TRUE)

if(args$transform == "SCT"){
	if(!is.null(args$cc_gene)){
		data_integrated <- SCTransform(data_integrated, method = "glmGamPoi", vars.to.regress = c("temp_sample","percent.mt", "S.Score", "G2M.Score"), variable.features.n=nf.usage, verbose=TRUE)
	}else{
		data_integrated <- SCTransform(data_integrated, method = "glmGamPoi", vars.to.regress = c("temp_sample","percent.mt"), variable.features.n=nf.usage, verbose=TRUE)
	}
}else{
	data_integrated <- NormalizeData(data_integrated)
	data_integrated <- FindVariableFeatures(data_integrated,selection.method = "vst", nfeatures = nf.usage)
	if(!is.null(args$cc_gene)){
		data_integrated <- ScaleData(data_integrated, vars.to.regress = c("temp_sample","percent.mt", "S.Score", "G2M.Score"))
	}else{
		data_integrated <- ScaleData(data_integrated, vars.to.regress = c("temp_sample","percent.mt"))
	}
}
data_integrated <- RunPCA(data_integrated, npcs = dim.all)


if(args$transform == "SCT"){
	data_integrated <- RunHarmony(data_integrated, group.by.vars = "temp_sample", assay.use="SCT", reduction = "pca", dims.use=1:dim.usage1, reduction.save = "harmony", theta = args$theta, lambda = args$lambda, sigma = args$sigma, nclust=args$nclust, block.size = args$blocksize, max.iter.harmony = 30, max.iter.cluster = 40, epsilon.cluster = 1e-06, epsilon.harmony = 1e-05, plot_convergence = TRUE)
}else{
	data_integrated <- RunHarmony(data_integrated, group.by.vars = "temp_sample", reduction = "pca", dims.use=1:dim.usage1, reduction.save = "harmony", theta = args$theta, lambda = args$lambda, sigma = args$sigma, nclust=args$nclust, block.size = args$blocksize, max.iter.harmony = 30, max.iter.cluster = 40, epsilon.cluster = 1e-06, epsilon.harmony = 1e-05, plot_convergence = TRUE)
}

data_integrated <- RunUMAP(data_integrated, reduction = "harmony", dims=1:dim.usage2)
data_integrated <- RunTSNE(data_integrated, reduction = "harmony", dims=1:dim.usage2)
data_integrated <- FindNeighbors(data_integrated, reduction='harmony', dims=1:dim.usage2)
if(args$algorithm == "leiden"){
	data_integrated <- FindClusters(data_integrated, resolution = res.usage, algorithm=4, method='igraph')
}else{
	data_integrated <- FindClusters(data_integrated, resolution = res.usage)
}

write.table(data_integrated@meta.data,paste0(args$output,"/",args$sample,".cell_info.txt"),sep="\t", quote=FALSE)
saveRDS(data_integrated,paste(args$output,"/",args$sample,".integrated.RDS",sep=""))


p1 <- DimPlot(data_integrated, reduction = "umap", label=T, raster=F)
p2 <- DimPlot(data_integrated, reduction = "tsne", label=T, raster=F)
p <- p1|p2
ggsave(filename = paste0(args$output, "/", args$sample,".Umap_TSNE.pdf"), plot = p, device = "pdf", width = 16, height = 7)

nfeature_quantile <- quantile(data_integrated$nFeature_RNA, probs = c(.2, .8))
p1 <- FeaturePlot(data_integrated, "nFeature_RNA", cols=c("blue","yellow","red"), min.cutoff=nfeature_quantile[[1]], max.cutoff=nfeature_quantile[[2]], reduction = "umap")
p2 <- FeaturePlot(data_integrated, "nFeature_RNA", cols=c("blue","yellow","red"), min.cutoff=nfeature_quantile[[1]], max.cutoff=nfeature_quantile[[2]], reduction = "tsne")
p <- p1|p2
ggsave(filename = paste0(args$output, "/", args$sample, ".Umap_TSNE.nfeature.pdf"), plot = p, device = "pdf", width = 15, height = 7)

p <- VlnPlot(data_integrated, features=c("nCount_RNA","nFeature_RNA"), pt.size=0, ncol=1)
ggsave(filename = paste0(args$out, "/", args$sample, ".nFeature.vlnplot.pdf"), plot = p, device = "pdf", width = length(levels(data_integrated$seurat_clusters))/4, height = 8)

p1 <- DimPlot(data_integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(data_integrated, reduction = "tsne", group.by = "orig.ident")
p <- p1|p2
ggsave(filename = paste0(args$output, "/", args$sample, ".Umap_TSNE.orig_ident.pdf"), plot = p, device = "pdf", width = 15, height = 7)

p1 <- DimPlot(data_integrated, reduction = "umap", group.by = "temp_sample")
p2 <- DimPlot(data_integrated, reduction = "tsne", group.by = "temp_sample")
p <- p1|p2
ggsave(filename = paste0(args$output, "/", args$sample, ".Umap_TSNE.sample.pdf"), plot = p, device = "pdf", width = 15, height = 7)


if(!is.null(args$marker)){
	sample_marker <- read.table(args$marker,header=F)
	colnames(sample_marker) = c("marker", "anno")
	plot_theme<-theme(panel.background = element_blank(),axis.line = element_line(size=0.1),axis.ticks = element_line(size=0.1),axis.text = element_text(size=7),axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),axis.ticks.length = unit(1, "pt"),plot.title = element_text(size = 7, face = "bold",margin=margin(0,0,4,0)),plot.margin = unit(c(0,0,0,0),"cm"),legend.key.size = unit(0.1, 'cm'),legend.key.height = unit(0.1, 'cm'),legend.key.width = unit(0.1, 'cm'),legend.title = element_text(size=6),legend.text = element_text(size=6))

	dot_p <- DotPlot(data_integrated, features=sample_marker$marker,cols = c("grey","red"),col.min= -0.5,dot.scale=3)+labs(x="",y="")+plot_theme
	ggsave(filename = paste0(args$output, "/", args$sample, ".Marker.DotPlot.pdf"), plot = dot_p, device = "pdf", width = round(nrow(sample_marker)/8+2), height = 4)
}

markers <- FindAllMarkers(data_integrated, only.pos=TRUE)
write.table(markers,paste0(args$output,"/",args$sample,".marker_info.txt"),sep="\t", quote=FALSE)


