# Data is available in the GEO database under accession number GSE274338.
# Short summary of experimentation:
# PBMCs were isolated from four healthy human blood donors.
# PBMCs were stained with 5OPRU-hMR1-Tet and positively enriched before sorting of MAIT and non-MAIT cells.
# Libraries were generated following the 10X Genomics Next-GEM v2 5' Kit with Feature Barcode Technology and TCR sequencing.

# Setup -------------------------------------------------------------------

devtools::install:github('plger/scDblFinder') # DOI: 10.12688/f1000research.73600.2
devtools::install:github('powellgenomicslab/Nebulosa') # DOI: 10.1093/bioinformatics/btab003
libraries <- c('scDblFinder', 'Nebulosa')
required_packages <- c("readxl", "reshape2", "see", "purrr", "ggprism", "rstatix",
                       "scales", "stats", "RColorBrewer", "ggpubr",
                       "glmGamPoi", "UCell", "org.Hs.eg.db", "Seurat", "HGNChelper",
                       "ReactomePA", #Guangchuang Yu, Qing-Yu He. ReactomePA: Molecular BioSystems 2016, 12(2):477-479
                       "msigdbr", "DOSE", "fgsea", "clusterProfiler",  "MAST", "tidyverse")
for(package in required_packages){if(!require(package, character.only=TRUE))
{install.packages(package, dependencies=TRUE)
  library(package, character.only=TRUE) }}
lapply(libraries, library, character.only=TRUE)
rm(required_packages, libraries)

script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)

# Cellranger ouput processing ------------------------------------------------------

# cellranger-identified cell barcodes
cells_D1 <- Read10X_h5('../D1_multi_2023_06_13/count/sample_filtered_feature_bc_matrix.h5')
cells_D2 <- Read10X_h5('../D2_multi_2023_06_13/count/sample_filtered_feature_bc_matrix.h5')
cells_D3 <- Read10X_h5('../D3_multi_2023_06_13/count/sample_filtered_feature_bc_matrix.h5')
cells_D4 <- Read10X_h5('../D4_multi_2023_06_13/count/sample_filtered_feature_bc_matrix.h5')
cells_list <- list(D1=cells_D1, D2=cells_D2, D3=cells_D3, D4=cells_D4)

adt    <- c('CD4_TotalSeqC', 'CD8_TotalSeqC', 'CD16_TotalSeqC', 'CD27_TotalSeqC', 'CD28_TotalSeqC', 'CD45RA_TotalSeqC', 'CD45RO_TotalSeqC', 'CD56_TotalSeqC', 'CD62L_TotalSeqC', 'CD161_TotalSeqC', 'CCR7_TotalSeqC', 'Va72_TotalSeqC')
hto    <- c('Hashtag1', 'Hashtag3') # #1 == MR1-reactive, #3 == eluted MR1-non-reactive fraction
donors <- c('D1', 'D2', 'D3', 'D4')

# separate HTO from ADT, create SeuratObject for each donor
seurat_list <- list()
for (i in donors) {
  hto_data     <- get(paste0('cells_', i))$`Antibody Capture`[hto, ]
  adt_data     <- get(paste0('cells_', i))$`Antibody Capture`[adt, ]
  rna_data <- get(paste0('cells_', i))$`Gene Expression`
  seurat_obj          <- CreateSeuratObject(counts=rna_data, project=i)
  seurat_obj          <- AddMetaData(seurat_obj, col.name='donor', metadata=i)
  seurat_obj[['ADT']] <- CreateAssayObject(counts=adt_data)
  seurat_obj[['HTO']] <- CreateAssayObject(counts=hto_data)
  seurat_list[[i]]<- seurat_obj
  rm(hto_data, adt_data, seurat_obj)
}
rm(cells_D1, cells_D2, cells_D3, cells_D4)

# Add Metadata
seurat_D1 <- seurat_list[['D1']] %>% AddMetaData(metadata=list(donor='D1', batch='Batch1', age=57, AB0='B', gender='female', rhd='RhDpos'))
seurat_D2 <- seurat_list[['D2']] %>% AddMetaData(metadata=list(donor='D2', batch='Batch2', age=58, AB0='A', gender='male', rhd='RhDpos'))
seurat_D3 <- seurat_list[['D3']] %>% AddMetaData(metadata=list(donor='D3', batch='Batch2', age=53, AB0='0', gender='male', rhd='RhDpos'))
seurat_D4 <- seurat_list[['D4']] %>% AddMetaData(metadata=list(donor='D4', batch='Batch2', age=55, AB0='0', gender='male', rhd='RhDpos'))
rm(seurat_list)

# data deposition
seurat_D1_RNA <- GetAssayData(seurat_D1, assay='RNA', slot='counts')
seurat_D1_ADT <- GetAssayData(seurat_D1, assay='ADT', slot='counts')
seurat_D1_HTO <- GetAssayData(seurat_D1, assay='HTO', slot='counts')

seurat_D2_RNA <- GetAssayData(seurat_D2, assay='RNA', slot='counts')
seurat_D2_ADT <- GetAssayData(seurat_D2, assay='ADT', slot='counts')
seurat_D2_HTO <- GetAssayData(seurat_D2, assay='HTO', slot='counts')

seurat_D3_RNA <- GetAssayData(seurat_D3, assay='RNA', slot='counts')
seurat_D3_ADT <- GetAssayData(seurat_D3, assay='ADT', slot='counts')
seurat_D3_HTO <- GetAssayData(seurat_D3, assay='HTO', slot='counts')

seurat_D4_RNA <- GetAssayData(seurat_D4, assay='RNA', slot='counts')
seurat_D4_ADT <- GetAssayData(seurat_D4, assay='ADT', slot='counts')
seurat_D4_HTO <- GetAssayData(seurat_D4, assay='HTO', slot='counts')

# save as rds files
saveRDS(seurat_D1_RNA, file='./Data deposition/Donor1_RNA_Kammann_et_al.rds')
saveRDS(seurat_D2_RNA, file='./Data deposition/Donor2_RNA_Kammann_et_al.rds')
saveRDS(seurat_D3_RNA, file='./Data deposition/Donor3_RNA_Kammann_et_al.rds')
saveRDS(seurat_D4_RNA, file='./Data deposition/Donor4_RNA_Kammann_et_al.rds')

saveRDS(seurat_D1_ADT, file='./Data deposition/Donor1_ADT_Kammann_et_al.rds')
saveRDS(seurat_D2_ADT, file='./Data deposition/Donor2_ADT_Kammann_et_al.rds')
saveRDS(seurat_D3_ADT, file='./Data deposition/Donor3_ADT_Kammann_et_al.rds')
saveRDS(seurat_D4_ADT, file='./Data deposition/Donor4_ADT_Kammann_et_al.rds')

saveRDS(seurat_D1_HTO, file='./Data deposition/Donor1_HTO_Kammann_et_al.rds')
saveRDS(seurat_D2_HTO, file='./Data deposition/Donor2_HTO_Kammann_et_al.rds')
saveRDS(seurat_D3_HTO, file='./Data deposition/Donor3_HTO_Kammann_et_al.rds')
saveRDS(seurat_D4_HTO, file='./Data deposition/Donor4_HTO_Kammann_et_al.rds')

rm(seurat_D1_RNA, seurat_D2_RNA, seurat_D3_RNA, seurat_D4_RNA,
   seurat_D1_ADT, seurat_D2_ADT, seurat_D3_ADT, seurat_D4_ADT,
   seurat_D1_HTO, seurat_D2_HTO, seurat_D3_HTO, seurat_D4_HTO)

# HTODemux ----------------------------------------------------------------

# developed by Stoeckius et al., DOI: 10.1186/s13059-018-1603-1

seurat_list <- list(D1=seurat_D1, D2=seurat_D2, D3=seurat_D3, D4=seurat_D4)
seurat_list_post_HTODemux <- list()

for (i in donors) {
  seurat_obj <- seurat_list[[i]]
  m <- GetAssayData(seurat_obj, assay='HTO')
  m <- m+1
  HTO <- CreateAssayObject(counts=m)
  seurat_obj[['HTO']] <- HTO
  seurat_obj <- NormalizeData(seurat_obj, assay='HTO', normalization.method='CLR')
  seurat_obj <- HTODemux(seurat_obj, assay='HTO', positive.quantile=0.99, seed=42)
  seurat_list_post_HTODemux[[i]] <- seurat_obj
  rm(m, seurat_obj, HTO)
}

# subset singlets
table(seurat_list_post_HTODemux[["D1"]]@meta.data[["HTO_classification.global"]])
table(seurat_list_post_HTODemux[["D2"]]@meta.data[["HTO_classification.global"]])
table(seurat_list_post_HTODemux[["D3"]]@meta.data[["HTO_classification.global"]])
table(seurat_list_post_HTODemux[["D4"]]@meta.data[["HTO_classification.global"]])

cells_singlets <- list()
for (i in donors){
  seurat_obj          <- seurat_list_post_HTODemux[[i]]
  seurat_singlet      <- subset(seurat_obj, HTO_classification.global=='Singlet')
  cells_singlets[[i]] <- seurat_singlet
  rm(seurat_obj, seurat_singlet)
}

rm(seurat_list, seurat_list_post_HTODemux, hto, seurat_D1, seurat_D2, seurat_D3, seurat_D4)

# scDblFinder -------------------------------------------------------------

# developed by Germain et al. (2022), DOI: 10.12688/f1000research.73600.2

singlets <- list() # singlets post Cellranger & post HTODemux & post scDblFinder
for (i in donors){
  X             <- cells_singlets[[i]]
  X_singlets    <- scDblFinder(GetAssayData(X, slot='counts'))
  X$scDblFinder.class <- X_singlets$scDblFinder.class
  X             <- subset(X, scDblFinder.class=='singlet')
  singlets[[i]] <- X
  rm(X, X_singlets)
}

# QC ----------------------------------------------------------------------

cells <- merge(singlets[['D1']], c(singlets[['D2']], singlets[['D3']], singlets[['D4']]), add.cell.ids=c('D1', 'D2', 'D3', 'D4'))

# calculate quality control metrics
cells <- PercentageFeatureSet(cells, '^MT-', col.name='percent_mito')
cells <- PercentageFeatureSet(cells, '^RP[SL]', col.name='percent_ribo')

feature_qc <- c('nFeature_RNA', 'nCount_RNA', 'percent_mito', 'percent_ribo')
VlnPlot(cells, group.by='donor', features=feature_qc, pt.size=0.05, ncol=2)

# keep high quality cells
cells <- subset(cells, nFeature_RNA>200 & nCount_RNA>1000 & nCount_RNA<8000 & percent_mito<7.5 & percent_ribo>10)

# adjust ADT tag names
cells@assays[["ADT"]]@data@Dimnames[[1]]   <- sub('-TotalSeqC', '', cells@assays[['ADT']]@data@Dimnames[[1]])
cells@assays[["ADT"]]@counts@Dimnames[[1]] <- sub('-TotalSeqC', '', cells@assays[['ADT']]@counts@Dimnames[[1]])

# Cells: Initial Integration ------------------------------------------------

# RNA (standard Seurat pipeline)
cells.list     <- SplitObject(cells, split.by='donor')
cells.list     <- lapply(cells.list, FUN=function(x) {x <- NormalizeData(object=x, assay='RNA') } )
cells.list     <- lapply(cells.list, FUN=function(x) {x <- FindVariableFeatures(object=x, assay='RNA', nfeatures=2000, selection.method='vst') } )
cells.list     <- lapply(cells.list, FUN=function(x) {x <- ScaleData(object=x, assay='RNA') } )
cells.features <- SelectIntegrationFeatures(object.list=cells.list)
cells.anchors  <- FindIntegrationAnchors(object.list=cells.list, anchor.features=cells.features)
int.cells      <- IntegrateData(anchorset=cells.anchors, normalization.method='LogNormalize', new.assay.name='RNAint')

int.cells <- NormalizeData(int.cells, assay='ADT', normalization.method='CLR', margin=2)
int.cells <- ScaleData(int.cells, assay='ADT')
rm(cells.list, cells.features, cells.anchors, cells)

# Cells: Initial PCA and UMAP ---------------------------------------------------------------------

DefaultAssay(int.cells) <- 'RNAint'
int.cells <- ScaleData(int.cells)
int.cells <- RunPCA(int.cells, npcs=50, reduction.name='rna_pca', reduction.key='rnaPC_')
int.cells <- ProjectDim(int.cells, reduction='rna_pca')
ElbowPlot(int.cells, ndims=50, reduction='rna_pca')
int.cells <- RunUMAP(int.cells, reduction='rna_pca', reduction.name='rna_umap', reduction.key='rnaUMAP_', dims=1:10)
int.cells <- FindNeighbors(int.cells, dims=1:10, reduction='rna_pca')
int.cells <- FindClusters (int.cells, resolution=1)

Idents(int.cells) <- int.cells@meta.data[['RNAint_snn_res.1']]
DimPlot(int.cells, pt.size=0.2, reduction='rna_umap', label=TRUE)+NoLegend()
ggsave(filename='UMAP/DimPlot_rnaUMAP_cells.png', path=plotdir, width=100, height=100, units='mm')

# Cells: Celltype annotation with sc-type ----------------------------------------

# Ianevski et al: https://doi.org/10.1038/s41467-022-28803-w
source('https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R') # load gene set preparation function
source('https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R') # load cell type annotation function
db_ <- './sc-type annotation/ScType_withMAIT.xlsx' # Load reference dataset with MAIT cell customization
tissue        <- 'Immune system' # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus
gs_list       <- gene_sets_prepare(db_, tissue)
es.max        <- sctype_score(scRNAseqData=int.cells[['RNA']]@scale.data, scaled=TRUE, gs=gs_list$gs_positive, gs2=gs_list$gs_negative)

cL_results    <- do.call('rbind', lapply(unique(int.cells@meta.data[['RNAint_snn_res.1']]), function(cl){
  es.max.cl=sort(rowSums(es.max[ ,rownames(int.cells@meta.data[int.cells@meta.data[['RNAint_snn_res.1']]==cl, ])]), decreasing=!0)
  head(data.frame(cluster=cl, type=names(es.max.cl), scores=es.max.cl, ncells=sum(int.cells@meta.data[['RNAint_snn_res.1']]==cl)), 10) }))

sctype_scores <- cL_results %>% group_by(cluster) %>% top_n(n=1, wt=scores)
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4]='Unknown'

# annotate sc-type result in metadata
int.cells@meta.data$sctype <-  ''
for(j in unique(sctype_scores$cluster)){
  cl_type=sctype_scores[sctype_scores$cluster==j,];
  int.cells@meta.data$sctype[int.cells@meta.data[['RNAint_snn_res.1']]==j] <- as.character(cl_type$type[1])}

DimPlot(int.cells, pt.size=0.2, reduction='rna_umap', group.by='sctype', label=TRUE)+NoLegend()

# remove contaminating monocytes and B cells
TNK <- subset(int.cells, sctype!='Non-classical monocytes' & sctype!='Immature B cells')
rm(int.cells, cl_type, cL_results, es.max, gs_list, tissue, db_, sctype_scores, sctype_score, gene_sets_prepare)

# TNK: Integration post filtering ------------------------------------------------

# RNA (standard Seurat pipeline)
TNK.list     <- SplitObject(TNK, split.by='donor')
TNK.list     <- lapply(TNK.list, FUN=function(x) {x <- NormalizeData(object=x, assay='RNA') } )
TNK.list     <- lapply(TNK.list, FUN=function(x) {x <- FindVariableFeatures(object=x, assay='RNA', nfeatures=5000, selection.method='vst') } )
TNK.list     <- lapply(TNK.list, FUN=function(x) {x <- ScaleData(object=x, assay='RNA') } )
TNK.features <- SelectIntegrationFeatures(object.list=TNK.list)
TNK.anchors  <- FindIntegrationAnchors(object.list=TNK.list, anchor.features=TNK.features)
int.TNK      <- IntegrateData(anchorset=TNK.anchors, normalization.method='LogNormalize', new.assay.name='RNAint')

int.TNK <- NormalizeData(int.TNK, assay='ADT', normalization.method='CLR', margin=2)
int.TNK <- ScaleData(int.TNK, assay='ADT')

# TNK: PCA and UMAP post filtering ---------------------------------------------------------------------

DefaultAssay(int.TNK) <- 'RNAint'
int.TNK <- ScaleData(int.TNK)
int.TNK <- RunPCA(int.TNK, npcs=50, reduction.name='rna_pca', reduction.key='rnaPC_')
int.TNK <- ProjectDim(int.TNK, reduction='rna_pca')
int.TNK <- RunUMAP(int.TNK, reduction='rna_pca', reduction.name='rna_umap', reduction.key='rnaUMAP_', dims=1:10)
int.TNK <- FindNeighbors(int.TNK, dims=1:10, reduction='rna_pca')
int.TNK <- FindClusters (int.TNK, resolution=1)
Idents(int.TNK) <- int.TNK@meta.data[['RNAint_snn_res.1']]

# TNK: Cluster analysis / DEG --------------------------------------------------------

# DEG analysis
int.TNK_rna_markers <- FindAllMarkers(int.TNK, assay='RNA', test.use='wilcox', only.pos=TRUE, min.pct=0.3, logfc.threshold=0.3)
int.TNK_rna_markers <- int.TNK_rna_markers %>% arrange(cluster, desc(avg_log2FC))

Idents(int.TNK) <- int.TNK@meta.data[["RNAint_snn_res.1"]]
int.TNK <- RenameIdents(int.TNK,
                        '0'='MAIT cells', '1'='MAIT cells', '5'='MAIT cells', '6'='MAIT cells', '7'='MAIT cells',
                        '12'='NK cells',
                        '13'='regulatory CD4 T cells',
                        '2'='naive CD4 T cells',
                        '3'='memory CD4 T cells',
                        '11'='naive CD8 T cells',
                        '4'='cytotoxic CD8 T cells',
                        '9'='memory CD8 T cells', '14'='memory CD8 T cells',
                        '10'='gd T cells', '8'='gd T cells')
int.TNK@meta.data$celltype <- Idents(int.TNK)

# Fig. 5A
colors_clusters <- c('MAIT cells'='#608080', 'NK cells'='#472D7BFF',
                     'regulatory CD4 T cells'='#3E356BFF','naive CD4 T cells'='#8BBEDA', 'memory CD4 T cells'='#3B5698',
                     'naive CD8 T cells'='#F6BB97FF', 'memory CD8 T cells'='#F4875EFF', 'cytotoxic CD8 T cells'='#CB1B4FFF',
                     'gd T cells'='#96DDB5FF')
DimPlot(int.TNK, pt.size=0.2, reduction='rna_umap', cols=colors_clusters, label=FALSE, shuffle=TRUE)
plot_density(int.TNK, features=c('adt_Va72', 'adt_CD161'), reduction='rna_umap', combine=TRUE, joint=FALSE)

# Fig. S9A
top10_TNK_rna_markers <- int.TNK_rna_markers %>% filter(p_val_adj<0.05) %>% group_by(cluster) %>% top_n(10, avg_log2FC) %>% mutate(rank=base::rank(dplyr::desc(avg_log2FC), ties.method='random')) %>% arrange(cluster, rank, avg_log2FC)
int.TNK %>% subset(downsample=200) %>% DoHeatmap(features=top10_TNK_rna_markers$gene, assay='RNA', slot='scale.data', size=4, group.bar=TRUE, group.colors=colors_clusters) & NoLegend()

# MAIT cell CD56 phenotyping ------------------------------------------------------

maits      <- subset(int.TNK, celltype=='MAIT cells')
maits.list <- SplitObject(maits, split.by='donor')

# check equal anti-human CD56 antibody signal
D1 <- maits.list[['D1']]
D2 <- maits.list[['D2']]
D3 <- maits.list[['D3']]
D4 <- maits.list[['D4']]
FeatureScatter(D1, feature1='adt_CD56', feature2='adt_CD8', slot='scale.data', pt.size=1, raster=FALSE, jitter=TRUE, plot.cor=FALSE)+geom_vline(xintercept=0)+geom_vline(xintercept=0.5)
FeatureScatter(D2, feature1='adt_CD56', feature2='adt_CD8', slot='scale.data', pt.size=1, raster=FALSE, jitter=TRUE, plot.cor=FALSE)+geom_vline(xintercept=0)+geom_vline(xintercept=0.5)
FeatureScatter(D3, feature1='adt_CD56', feature2='adt_CD8', slot='scale.data', pt.size=1, raster=FALSE, jitter=TRUE, plot.cor=FALSE)+geom_vline(xintercept=0)+geom_vline(xintercept=0.5)
FeatureScatter(D4, feature1='adt_CD56', feature2='adt_CD8', slot='scale.data', pt.size=1, raster=FALSE, jitter=TRUE, plot.cor=FALSE)+geom_vline(xintercept=0)+geom_vline(xintercept=0.5)
rm(D1, D2, D3, D4)

# gate CD56 for each donor
for (i in donors){
  X <- maits.list[[i]]
  CD56pos_cells <- WhichCells(X, slot='scale.data', expression=adt_CD56>=0.5 | rna_NCAM1>0) # CD56+ MAIT cells express NCAM1 or have above threshold adt_CD56 signal
  CD56dim_cells <- WhichCells(X, slot='scale.data', expression=adt_CD56<0.5 & adt_CD56>=0)
  X$CD56status  <- ifelse(colnames(X) %in% CD56pos_cells, yes='CD56pos', no=ifelse(colnames(X) %in% CD56dim_cells, yes='CD56dim', no='CD56neg'))
  maits.list[[i]] <- X
  rm(X, CD56pos_cells, CD56dim_cells)
}

Idents(maits) <- 'CD56status'
CD56pos_cells  <- WhichCells(maits, idents='CD56pos')
CD56neg_cells  <- WhichCells(maits, idents='CD56neg')
CD56dim_cells  <- WhichCells(maits, idents='CD56dim')
int.TNK$CD56status <- ifelse(colnames(int.TNK) %in% CD56pos_cells, yes='CD56pos', no=ifelse(colnames(int.TNK) %in% CD56dim_cells, yes='CD56dim', no=ifelse(colnames(int.TNK) %in% CD56neg_cells, yes='CD56neg', no=NA)) )

# Fig. 5B
maits %>% subset(downsample=2000) %>% FeatureScatter(feature1='adt_CD8', feature2='adt_CD56', slot='scale.data', group.by='CD56status', col=mycolors3, pt.size=0.2, plot.cor=FALSE)+geom_point(shape=21, alpha=0.005)+scale_y_continuous(limits=c(-1, 9), breaks=c(0,4,8))+scale_x_continuous(limits=c(-3,2), breaks=c(-2,0,2))+theme(panel.grid.major.x=element_line(linewidth=0.35, linetype='dotted', color='grey'))+NoLegend()
maits %>% subset(downsample=2000) %>% FeatureScatter(feature1='adt_CD8', feature2='rna_NCAM1', slot='scale.data', group.by='CD56status', col=mycolors3, pt.size=0.2, plot.cor=FALSE)+geom_point(shape=21, alpha=0.005)+scale_y_continuous(limits=c(-0.5, 10.5), breaks=c(0,5,10))+scale_x_continuous(limits=c(-3,2), breaks=c(-2,0,2))+theme(panel.grid.major.x=element_line(linewidth=0.35, linetype='dotted', color='grey'))+NoLegend()

# data deposition
RNA_matrix <- int.TNK@assays$RNA@counts
ADT_matrix <- int.TNK@assays$ADT@counts
HTO_matrix <- GetAssayData(int.TNK, assay='HTO', slot='counts')
metadata   <- int.TNK@meta.data

# write matrix to .rds files
saveRDS(RNA_matrix, file='./Data deposition/RNA_Kammann_et_al.rds')
saveRDS(ADT_matrix, file='./Data deposition/ADT_Kammann_et_al.rds')
saveRDS(HTO_matrix, file='./Data deposition/HTO_Kammann_et_al.rds')
saveRDS(metadata,   file='./Data deposition/Metadata_Kammann_et_al.rds')

rm(CD56pos_cells, CD56neg_cells, CD56dim_cells, int.TNK, RNA_matrix, ADT_matrix, HTO_matrix, metadata)

# MAIT cell clustering ----------------------------------------------------

DefaultAssay(maits.list) <- 'RNA'
maits.list     <- lapply(maits.list, FUN=function(x) {x <- NormalizeData(object=x, assay='RNA') } )
maits.list     <- lapply(maits.list, FUN=function(x) {x <- FindVariableFeatures(object=x, assay='RNA', nfeatures=5000, selection.method='vst') } )
maits.list     <- lapply(maits.list, FUN=function(x) {x <- ScaleData(object=x, assay='RNA') } )
maits.features <- SelectIntegrationFeatures(object.list=maits.list, nfeatures=5000)
maits.anchors  <- FindIntegrationAnchors(object.list=maits.list, anchor.features=maits.features)
maits          <- IntegrateData(anchorset=maits.anchors, normalization.method='LogNormalize', new.assay.name='RNAint')
rm(maits.anchors, maits.list)

DefaultAssay(maits) <- 'RNAint'
maits <- Trex::quietTCRgenes(maits) # silence TCR chain genes from VariableFeatures present in RNAint assay
maits <- ScaleData(maits)
maits <- RunPCA(maits, npcs=50, reduction.name='rna_pca', reduction.key='rnaPC_')
maits <- ProjectDim(maits, reduction='rna_pca')
ElbowPlot(maits, ndims=50, reduction='rna_pca')
maits <- RunUMAP(maits, reduction='rna_pca', reduction.name='rna_umap', reduction.key='rnaUMAP_', dims=1:8)
maits <- FindNeighbors(maits, dims=1:8, reduction='rna_pca')
maits <- FindClusters (maits, resolution=0.5)
DimPlot(maits, pt.size=0.5, reduction='rna_umap', label=TRUE)+NoLegend()

# WNN clustering
DefaultAssay(maits)     <- 'ADT'
adt_features            <- sub('_TotalSeqC', '', adt)
VariableFeatures(maits) <-rownames(maits@assays[['ADT']]@meta.features)

# initialize pseudo-PCA slot for ADT assay, scale ADT data, load ADT scale.data into PCA dimensions
maits <- ScaleData(maits, assay='ADT')
maits <- RunPCA(maits, reduction.name='adt_pca', reduction.key='adtPC_', features=VariableFeatures(maits), verbose=FALSE)
pseudo_adt <- t(GetAssayData(maits, assay='ADT', slot='scale.data'))
colnames(pseudo_adt) <- paste('pseudo', 1:length(adt_features), sep='_')
maits@reductions$adt_pca@cell.embeddings <- pseudo_adt

maits <- FindMultiModalNeighbors(maits, reduction.list=list('adt_pca', 'rna_pca'), dims.list=list(1:12, 1:8), modality.weight.name=c('ADT.weight', 'RNA.weight'))
maits <- RunUMAP(maits, nn.name='weighted.nn', reduction.name='wnn_umap', reduction.key='wnnUMAP_')
maits <- FindClusters(maits, graph.name='wsnn', algorithm=3, resolution=0.2)
Idents(maits) <- maits@meta.data[['wsnn_res.0.2']]

# Fig. 5F
DimPlot(maits, pt.size=0.2, reduction='wnn_umap', label=TRUE, group.by='wsnn_res.0.2')+NoLegend()
rm(adt_features, adt, pseudo_adt)

# MAIT cell CD56 DEG analysis with MAST test ---------------------------------------------

# Compare CD56-MAIT cells with CD56+ MAIT cells
Idents(maits)       <- 'CD56status'
DefaultAssay(maits) <- 'RNA'
cd56maits           <- subset(maits, CD56status!='CD56dim')

cd56_rna_markers_mast <- FindAllMarkers(cd56maits, assay='RNA', test.use='MAST', latent.vars='donor',  only.pos=FALSE, min.pct=0.01, logfc.threshold=0)
cd56_rna_markers_mast <- cd56_rna_markers_mast %>% adjust_pvalue(method='BH', p.col='p_val', output.col='adjust_p') %>% arrange(cluster, desc(avg_log2FC))
cd56_rna_markers_mast_CD56pos  <- cd56_rna_markers_mast %>% filter(cluster=='CD56pos') %>% View()

# gene annotations for volcano plot
x              <- cd56_rna_markers_mast %>% filter(cluster=='CD56pos')
xlabels        <- x %>% filter(avg_log2FC>0.2 | avg_log2FC< (-0.2) | p_val_adj<1e-30) %>% pull(gene)
alllabels      <- x %>% filter(avg_log2FC>0.2 | avg_log2FC< (-0.2) | p_val_adj<1e-30) %>% pull(gene)
labels_to_keep <- c('NCAM1', 'NKG7', 'GNLY', 'AOAH', 'KLRD1', 'TMIGD2', 'PRF1', 'CST7', 'CCR7', 'S100B', 'S100A', 'IL2RB', 'ID2', 'CEBPD', 'TYROBP', 'PRDM13', 'IL1R1', 'KLRF1', 'CD160', 'IFITM3', 'SCGB3A1', 'MAL', 'LEF1', 'CD4', 'S100A9', 'AC114930.1', 'PI16')
xlabels        <- xlabels[xlabels %in% labels_to_keep]

# Fig. 5C
# DEG Volcano CD56-MAIT vs CD56+MAIT
ggplot(data=x, aes(x=avg_log2FC, y=-log10(adjust_p)))+
  geom_vline(xintercept=c(-0.2, 0.2), linetype="dotted", color='black')+
  ggrepel::geom_label_repel(data=subset(x, gene %in% xlabels), aes(label=gene), box.padding=0.2, label.padding=0.15, size=2, show.legend=FALSE, max.overlaps=50)+
  geom_point(aes(fill=avg_log2FC), shape=21, size=2, alpha=1, show.legend=TRUE)+
  scale_x_continuous(limits=c(-4,4), breaks=c(-4,-2,0,2,4),  guide='prism_offset', expand=c(0,0) )+
  scale_y_continuous(limits=c(0,100), breaks=c(0,20,40,60,80,100), guide='prism_offset', expand=c(0.02, 0))+
  scale_fill_gradientn(colors=c('#C0C0C0','white','#008080'), na.value='transparent', breaks=c(-4, 4), limits=c(-4, 4), labels=c('CD56-\nMAIT', 'CD56+\nMAIT') )+
  theme(panel.grid.major.y=element_blank())+
  labs(x='Log2FC (CD56+ MAIT / CD56- MAIT)', y='-log10(adjusted p-value)', fill='', title='')

#saveRDS(cd56maits, file='cd56maits.rds')

# Pseudobulk generation for gene set and pathway enrichment analysis -----------------------------------------------------

maits_agg              <- AggregateExpression(cd56maits, return.seurat=TRUE, group.by='CD56status')
maits_agg$CD56status   <- Idents(maits_agg)
maits_agg_rna <- GetAssayData(maits_agg, assay='RNA', slot='data')
maits_agg_rna <- as.data.frame(maits_agg_rna)
maits_agg_rna$gene <- rownames(maits_agg_rna)
maits_agg_rna <- maits_agg_rna %>% mutate(diff=CD56pos-CD56neg)
xlabels <- xlabels[xlabels!='']

ggplot(maits_agg_rna, aes(CD56pos, CD56neg))+
  geom_point(aes(fill=diff),shape=21, color='black', alpha=0.8, size=2, show.legend=FALSE)+
  scale_fill_gradientn(colors=c('#C0C0C0','white','#008080'), na.value='transparent', breaks=c(-max(abs(maits_ave_rna$diff)), max(abs(maits_ave_rna$diff))), limits=c(-max(abs(maits_ave_rna$diff)), max(abs(maits_ave_rna$diff))), labels=c('CD56-\nMAIT', 'CD56+\nMAIT'))+
  ggrepel::geom_label_repel(data=subset(maits_ave_rna, gene %in% xlabels), aes(label=gene), box.padding=0.2, label.padding=0.15, size=2, show.legend=FALSE)+
  labs(x='mean CD56+ MAIT gene expression', y='mean CD56- MAIT gene expression')+
  scale_x_continuous(limits=c(0,6), guide='prism_offset', expand=c(0.02,0))+
  scale_y_continuous(limits=c(0,6), guide='prism_offset', expand=c(0.02,0))

# Gene ontology (GO) overrepresentation analysis (ORA) on MAST results -----------------------------------------------------------

# annotate DEGs with EntrezID
markers                <- cd56_rna_markers_mast_CD56pos %>% filter(avg_log2FC<(-0.2) | avg_log2FC>0.2)
markers_list           <- markers
markers_list$SYMBOL    <- markers_list$gene
markers_list           <- bitr(markers_list$SYMBOL, fromType='SYMBOL', toType='ENTREZID', OrgDb='org.Hs.eg.db', drop=T)
markers_list$gene      <- markers_list$SYMBOL
markers                <- left_join(markers, markers_list)

# create universe for ORA GO analysis (= all genes expressed by MAIT cells)
universe_background <- maits_agg_rna$gene
ensembl             <- useEnsembl(biomart='genes', dataset='hsapiens_gene_ensembl')
symbol_to_entrez    <- getBM(values=universe_background, attributes=c('external_gene_name', 'ensembl_gene_id', 'entrezgene_id'), filters='external_gene_name', mart=ensembl)
universe_background <- drop_na(symbol_to_entrez)
colnames(universe_background) <- c('SYMBOL', 'ENSEMBL', 'ENTREZID')

# take 100 most important DEGs
top100  <- markers %>% filter(adjust_p<0.001 & avg_log2FC>0.2 & cluster=='CD56pos') %>% top_n(n=100, wt=avg_log2FC) %>% arrange(desc(avg_log2FC))
top100  <- split(top100$ENTREZID, top100$cluster)

# GO ORA on biological processes (BP)
CD56pos_go        <- enrichGO(ont='BP', gene=top100[['CD56pos']], universe=universe_background, OrgDb=org.Hs.eg.db, pAdjustMethod='BH', pvalueCutoff=0.05, qvalueCutoff=0.05, readable=TRUE)
CD56pos_go_result <- CD56pos_go@result %>% as.data.frame() %>% filter(p.adjust<0.05) %>% View()

# Fig. 5D
selected_keys <- c('cell killing', 'positive regulation of cytokine production', 'leukocyte mediated immunity', 'leukocyte cell-cell adhesion','positive regulation of type II interferon production', 'regulation of inflammatory response', 'positive regulation of immune effector process')
CD56pos_go_selected <- CD56pos_go_result %>% filter(Description %in% selected_keys) %>%
ggplot(aes(y=reorder(Description, Count), x=1))+
  geom_point(aes(size=Count, fill=p.adjust), shape=21, color='black', show.legend=TRUE)+
  geom_label(aes(label=Count), nudge_x=0.02, label.padding=unit(0.1, 'lines'))+
  scale_x_continuous(limits=c(0.95,1.05), breaks=c(1), expand=c(0,0))+
  scale_y_discrete(expand=c(0.05,0.05), guide='prism_offset')+
  scale_fill_viridis()+
  labs(x='CD56+ MAIT cell GO ORA pathway', y='')+
  theme(panel.grid.major.y=element_line(linewidth=0.3, linetype='dotted'), axis.text.x=element_blank())

rm(CD56pos_go, CD56pos_go_result, CD56pos_go_selected, selected_keys)

# Reactome ORA analysis on MAST results ----------------------------

CD56pos_reactome        <- enrichPathway(gene=top100[['CD56pos']], pAdjustMethod='BH', pvalueCutoff=0.05, qvalueCutoff=0.2, readable=TRUE, universe=universe_background, organism='human', minGSSize=10, maxGSSize=103)
CD56pos_reactome_result <- CD56pos_reactome@result %>% as.data.frame() %>% filter(Count>=2) %>% adjust_pvalue(method='BH') %>% filter(pvalue.adj<0.05) %>% arrange(desc(Count)) %>% select(-p.adjust) %>% View()

# Fig. 5E
important <- c('Interleukin-12 family signaling', 'Costimulation by the CD28 family', 'Interleukin-10 signaling', 'Signal regulatory protein family interactions', 'FCGR activation')
CD56pos_reactome_result$Description <- factor(CD56pos_reactome_result$Description, levels=important)
ggplot(data=subset(CD56pos_reactome_result, Description%in%important), aes(y=reorder(Description, Count), x=1))+
  geom_point(aes(size=Count), shape=21, color='black', fill='grey', show.legend=FALSE)+
  geom_label(aes(label=Count), nudge_x=0.05)+
  scale_x_continuous(limits=c(0.9,1.1), breaks=c(1), expand=c(0,0))+
  scale_y_discrete(expand=c(0.05,0.05), guide='prism_offset')+
  labs(x='CD56+ MAIT cell reactome pathway', y='')+
  theme(panel.grid.major.y=element_line(linewidth=0.3, linetype='dotted'), axis.text.x=element_blank())

rm(CD56pos_reactome, CD56pos_reactome_result, important)

# GO gene set enrichment analysis (GSEA) ---------------------------------------

CD56pos_vs_CD56neg                 <- maits_agg_rna %>% mutate(SYMBOL=gene) %>% dplyr::select(-gene) %>% dplyr::mutate(delta_avg_exp=CD56pos-CD56neg) %>% arrange(desc(delta_avg_exp))
CD56pos_vs_CD56neg_delta           <- CD56pos_vs_CD56neg %>% dplyr::select(SYMBOL, delta_avg_exp) %>% arrange(desc(delta_avg_exp)) # control NES ranking
any(duplicated(CD56pos_vs_CD56neg_delta$SYMBOL)) # control for no duplicatess
CD56pos_vs_CD56neg_delta_FC        <- CD56pos_vs_CD56neg_delta$delta_avg_exp
names(CD56pos_vs_CD56neg_delta_FC) <- CD56pos_vs_CD56neg_delta$SYMBOL # fgsea requirement

# GO GSEA
gse <- gseGO(geneList=CD56pos_vs_CD56neg_delta_FC, ont='ALL', OrgDb='org.Hs.eg.db', keyType='SYMBOL', minGSSize=25, maxGSSize=250, pvalueCutoff=0.05, pAdjustMethod='BH', seed=42)
gse_result <- gse@result

# Visualization of GO GSEA
BPs <- c('cell killing', 'cytosolic ribosome', 'cytoplasmic translation', 'immunoglobulin complex')
ggplot(data=subset(gse_result, Description %in% BPs), aes(x=reorder(Description, NES), y=NES))+
  geom_col(aes(fill=p.adjust), color='black', show.legend=TRUE)+
  coord_flip()+
  scale_fill_viridis()+
  scale_y_continuous(limits=c(-abs(max(gse_result$NES)+0.2), abs(max(gse_result$NES)+0.2)))+
  labs(x='Gene Set', y='NES')

rm(CD56pos_vs_CD56neg, CD56pos_vs_CD56neg_delta, CD56pos_vs_CD56neg_delta_FC, BPs, gse, gse_result)

# Reactome GSEA on MAST results -------------------------

CD56pos_vs_CD56neg                 <- maits_agg_rna %>% mutate(SYMBOL=gene) %>% dplyr::select(-gene) %>% dplyr::mutate(delta_avg_exp=CD56pos-CD56neg) %>% arrange(desc(delta_avg_exp))
CD56pos_vs_CD56neg_2               <- left_join(CD56pos_vs_CD56neg, markers_list)
CD56pos_vs_CD56neg_delta           <- CD56pos_vs_CD56neg_2 %>% dplyr::select(ENTREZID, delta_avg_exp) %>% arrange(desc(delta_avg_exp)) %>% as.data.frame()
any(duplicated(CD56pos_vs_CD56neg_delta$ENTREZID)) # control for no duplicates
CD56pos_vs_CD56neg_delta           <- CD56pos_vs_CD56neg_delta %>% distinct(ENTREZID, .keep_all=TRUE) %>% drop_na()
CD56pos_vs_CD56neg_delta_FC        <- CD56pos_vs_CD56neg_delta$delta_avg_exp
names(CD56pos_vs_CD56neg_delta_FC) <- CD56pos_vs_CD56neg_delta$ENTREZID # fgsea requirement

# GSEA against Reactome database
gse        <- gsePathway(geneList=CD56pos_vs_CD56neg_delta_FC, organism='human', minGSSize=5, maxGSSize=800, pvalueCutoff=0.05, pAdjustMethod='BH', eps=0, seed=42)
gse_result <- gse@result

# Filtering of disease pathways
gse_filter <- c('Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell',
  'Signaling by Interleukins',
  'Cytokine Signalling in immune system',
  'Antimicrobial peptides',
  'Metabolism of amino acids and derivatives',
  'rRNA processing',
  'Cellular response to starvation',
  'Eukaryotic Translation Initiation',
  'Eukaryotic Translation Elongation',
  'Selenoamino acid metabolism'
)

# Fig. S9B
gse_result_filtered <- gse_result %>% filter(Description %in% gse_filter) %>%
ggplot(aes(x=reorder(Description, NES), y=NES))+
  geom_col(aes(fill=log10(p.adjust)),color='black', show.legend=TRUE)+
  coord_flip()+ scale_y_continuous(limits=c(-2.2,2.2))+ scale_fill_viridis()+ labs(fill='adjusted p')+ labs(x='Pathway', y='NES')

rm(list=ls())
gc()

