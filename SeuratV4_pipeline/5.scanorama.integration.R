## Rcode description: scRNA-seq integration by Seurat using Scanorama
##
## Author: Zhang Pei (zhangpei@genomics.cn)
##
## Date: 2022.8.16

parser = argparse::ArgumentParser(description="scRNA-seq integration by Seurat using Scanorama")
parser$add_argument('-I','--input', help = 'seurat object name and directory, without header')
parser$add_argument('-G','--mregex', help = 'regex for mitochondria')
parser$add_argument('-T','--transform', help = 'transform method, RNA or LogNormalize or SCT (default)')
parser$add_argument('-F','--nf',help = 'number of variable features, default 5000, useless for "-T RNA" ')
parser$add_argument('-D','--dim',help = 'PCA dim, default 100')
parser$add_argument('-U','--dusage',help = 'PCA dim usage, default 0.8*dim')
parser$add_argument('-R','--res',help = 'Map resolution usage, 1 as default')
parser$add_argument('-O','--out',help = 'out directory, ./ as default')
parser$add_argument("-M", '--marker',help= 'chosen marker list, NULL as default, file containing chosen marker genes for DotPlot, without header')
parser$add_argument("-C",'--cc_gene',help= 'input list for cell cycle genes, NULL as default, 1 for cc genes provided by seurat, or should be file including gene name and type (s.genes and g2m.genes) and without header')
parser$add_argument('-S','--sample',help = 'Sample name, also be used for slot name. There must be no replicate names for the total project')
parser$add_argument('--algorithm',help='clustering algorithm, leiden (default) or louvain')
parser$add_argument('--python',help='python3 bin path')


args = parser$parse_args()
args$mregex <- if(!is.null(args$mregex)) args$mregex else "Mitochondria"
args$mregex <- paste("^",args$mregex, sep="")
nf.usage <- as.numeric(if(!is.null(args$nf)) args$nf else 5000)
dim.all <- as.numeric(if(!is.null(args$dim)) args$dim else 100)
dim.usage <- as.numeric(if(!is.null(args$dusage)) args$dusage else dim.all)
res.usage <- as.numeric(if(!is.null(args$res)) args$res else 1)
args$out <- if(!is.null(args$out)) args$out else "./"
args$transform <- if(!is.null(args$transform)) args$transform else "SCT"
args$algorithm <- if(!is.null(args$algorithm)) args$algorithm else "leiden"
args$python <- if(!is.null(args$python)) args$python else "/ldfssz1/ST_DIVERSITY/PUB/USER/zhangpei/bin/python3.7.13/bin/python3"


library(Seurat)
library(ggplot2)

library(reticulate)
use_python(args$python, required = T)
py_config()
scanorama <- import('scanorama')

if(args$transform == "SCT"){
    library(sctransform)
    library(glmGamPoi)
}

if(args$algorithm == "leiden"){
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

for(i in 1:length(all_rds)) all_rds[[i]]@meta.data$temp_sample <- names(all_rds[i])
data_merge <- merge(x=all_rds[[1]], y=all_rds[2:length(all_rds)], merge.data = TRUE)


rm(all_rds)

all_rds <- SplitObject(data_merge, split.by = "temp_sample")


if(args$transform == "SCT"){
    all_rds <- lapply(X = all_rds, FUN = function(x){
	DefaultAssay(x) <- "RNA"
        x <- PercentageFeatureSet(x, pattern = args$mregex, col.name = "percent.mt")
	x <- NormalizeData(x)
        if(!is.null(args$cc_gene)){
            x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, assay = 'SCT', set.ident = TRUE)
            x <- SCTransform(x, method = "glmGamPoi", vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'), variable.features.n = nf.usage, verbose = FALSE)
        }else{
			x <- SCTransform(x, method = "glmGamPoi", vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA'), variable.features.n=nf.usage, verbose = FALSE)
		}
	})
}else if (args$transform == "LogNormalize"){
    all_rds <- lapply(X = all_rds, FUN = function(x) {
	DefaultAssay(x) <- "RNA"
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
  	})
}else{
	all_rds <- lapply(X = all_rds, FUN = function(x) {
		DefaultAssay(x) <- "RNA";
		x <- PercentageFeatureSet(x, pattern = args$mregex, col.name = "percent.mt")
		x <- NormalizeData(x)
	})
}


datasets <- list()
genes_list <- list()
for (i in 1:length(all_rds)) {
	if(args$transform == "SCT"){
		datasets[[i]] <- t(as.matrix(GetAssayData(all_rds[[i]], "scale.data", assay = "SCT")))	
	}else if(args$transform == "LogNormalize"){
		datasets[[i]] <- t(as.matrix(GetAssayData(all_rds[[i]], "scale.data", assay = "RNA")))
	}else{
		datasets[[i]] <- t(as.matrix(GetAssayData(all_rds[[i]], "counts", assay = "RNA")))
	}
	genes_list[[i]] <- colnames(datasets[[i]])
}

#integrated.corrected.data <- scanorama$integrate(datasets, genes_list, dimred= as.integer(dim.all))
integrated.corrected.data <- scanorama$correct(datasets, genes_list, return_dimred=TRUE, return_dense=TRUE, dimred=as.integer(dim.all))


intdata <- lapply(integrated.corrected.data[[2]], t)
corrected.matrix <- do.call(cbind, intdata)
rownames(corrected.matrix) <- as.character(integrated.corrected.data[[3]])
colnames(corrected.matrix) <- unlist(sapply(datasets, rownames))

intdimred <- do.call(rbind, integrated.corrected.data[[1]])
colnames(intdimred) <- paste0("PC_", 1:dim.all)
stdevs <- apply(intdimred, MARGIN = 2, FUN = sd)
data_merge@assays[["scanorama"]] <- CreateAssayObject(data = corrected.matrix, min.cells=0, min.features = 0)
data_merge@assays[["scanorama"]]@key <- "scanorama_"
rownames(intdimred) <- colnames(data_merge)
DefaultAssay(data_merge) <- "scanorama"
data_merge[["pca"]] <- CreateDimReducObject(embeddings = intdimred, stdev = stdevs, key = "PC_", assay = "scanorama")

for (i in Assays(data_merge)){ if(i != "RNA" & i != "scanorama") {data_merge[[i]] <- NULL}}


data_merge <- RunTSNE(data_merge, reduction = "pca", dims = 1:dim.usage)
data_merge <- RunUMAP(data_merge, reduction = "pca", dims = 1:dim.usage)
data_merge <- FindNeighbors(data_merge, reduction = "pca", dims = 1:dim.usage)

if(args$algorithm == "leiden"){
	data_merge <- FindClusters(data_merge, resolution = res.usage, algorithm=4, method='igraph')
}else{
	data_merge <- FindClusters(data_merge, resolution = res.usage)
}

write.table(data_merge@meta.data,paste0(args$out,"/",args$sample,".cell_info.txt"),sep="\t", quote=FALSE)

DefaultAssay(data_merge) <- "RNA"
data_merge <- NormalizeData(data_merge)
DefaultAssay(data_merge) <- "scanorama"

saveRDS(data_merge,paste(args$out,"/",args$sample,".integrated.RDS",sep=""))
			
p <- VlnPlot(data_merge, features=c("nCount_RNA","nFeature_RNA"), pt.size=0, ncol=1)
ggsave(filename = paste0(args$out, "/", args$sample, ".nFeature.vlnplot.pdf"), plot = p, device = "pdf", width = length(levels(data_merge$seurat_clusters))/4, height = 8)

p1 <- DimPlot(data_merge, reduction = "umap", label=T, raster=F)
p2 <- DimPlot(data_merge, reduction = "tsne", label=T, raster=F)
p <- p1|p2
ggsave(filename = paste0(args$out, "/", args$sample, ".Umap_TSNE.pdf"), plot = p, device = "pdf", width = 16, height = 7)

p1 <- DimPlot(data_merge, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(data_merge,  reduction = "tsne", group.by = "orig.ident")
p <- p1|p2
ggsave(filename = paste0(args$out, "/", args$sample, ".Umap_TSNE.orig_ident.pdf"), plot = p, device = "pdf", width = 15, height = 7)

if(!is.null(args$marker)){
	sample_marker <- read.table(args$marker,header=F)
	colnames(sample_marker) = c("marker", "anno")
	plot_theme<-theme(panel.background = element_blank(),axis.line = element_line(size=0.1),axis.ticks = element_line(size=0.1),axis.text = element_text(size=7),axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),axis.ticks.length = unit(1, "pt"),plot.title = element_text(size = 7, face = "bold",margin=margin(0,0,4,0)),plot.margin = unit(c(0,0,0,0),"cm"),legend.key.size = unit(0.1, 'cm'),legend.key.height = unit(0.1, 'cm'),legend.key.width = unit(0.1, 'cm'),legend.title = element_text(size=6),legend.text = element_text(size=6))

	dot_p <- DotPlot(data_merge, features=sample_marker$marker,cols = c("grey","red"),col.min= -0.5,dot.scale=3)+labs(x="",y="")+plot_theme
	ggsave(filename = paste0(args$out,"/", args$sample, ".Marker.DotPlot.pdf"), plot = dot_p, device = "pdf", width = round(nrow(sample_marker)/8+2), height = 4)
i}

markers <- FindAllMarkers(data_merge, assay= "RNA", only.pos=TRUE)
write.table(markers,paste0(args$out,"/",args$sample,".marker_info.txt"),sep="\t", quote=FALSE)

