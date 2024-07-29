# Unsupervised analysis (UMAP, CytoTree) ---------------------------------------------------

# Setup
if (!require('devtools', quietly=TRUE)) install.packages('devtools')
if (!require('flowCore', quietly=TRUE)) devtools::install:github("RGLab/flowCore") # Hahne et al. (2009), DOI: 10.1186/1471-2105-10-106
if (!require('flowUtils', quietly=TRUE)) devtools::install_github("jspidlen/flowUtils") #Spidlen et al.(2021)
if (!require('CytoTree', quietly=TRUE)) devtools::install_github("JhuangLab/CytoTree") # Dai et al. (2020), DOI: 10.1186/s12859-021-04054-2

required_packages <- c("readxl", "reshape2", "see", "purrr", "ggprism", "rstatix", "scales", "stats", "pheatmap", "RColorBrewer", "Hmisc", "corrplot", "vtable", "ggpubr", "tidyverse")
for(package in required_packages){
  if(!require(package, character.only=TRUE)) {install.packages(package, dependencies=TRUE)
  library(package, character.only=TRUE) }}
umap_libraries <- c("flowUtils", "flowCore", "CytoTree")
lapply(umap_libraries, library, character.only=TRUE)
rm(required_packages)
mycolors <- c('Blood'='#E88984', 'Spleen'='#A9261F', 'Liver'='#813c5e', 'Ileum'='#8c510a', 'Caecum'='#bf812d', 'Colon'='#dfc27d', 'mLN'='#f6e8c3', 'Lung'='#c7eae5', 'LungLN'='#92c5de', 'Skin'='#01662c')

script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)

# Fig 2:
### (optional): loading concatenated MAIT cells as fcs for UMAP
fcs_file  <- './IHOPE_UMAP_MST_analysis_Fig2.R'
fcs_data <- runExprsExtract(fcs_file, comp=FALSE, transformMethod='none')
fcs_protein_data <- fcs_data[, 8:34] # get protein values
fcs_meta_data <- fcs_data[,c(35,37)] # cut flowjo's numerical distribution into discrete metadata
fcs_meta_data[,"Pseudo_Donor<NA>"]  <- as.numeric(as.character(cut(fcs_meta_data[,"Pseudo_Donor<NA>"],  breaks=c(5500,6500,7500,8500,9500,10500,11500,12500,13500,14500,15500,16500,17500,18500), labels=c('6','7','8','9','10','11','12','13','14','15','16','17','18'), right=FALSE)))
fcs_meta_data[,"Pseudo_Tissue<NA>"] <- as.numeric(as.character(cut(fcs_meta_data[,"Pseudo_Tissue<NA>"], breaks=c(500,1500,2500,3500,4500,5500,6500,7500,8500,9500,10500), labels=c(1,2,3,4,5,6,7,8,9,10))))
fcs_meta_data <- as.data.frame(fcs_meta_data)
colnames(fcs_meta_data) <- c('Pseudo_Donor', 'Pseudo_Tissue')
fcs_protein_data_columns <- c("FJComp-APC-A<NA>"='Va7.2', "FJComp-APC-H7-A<NA>"='IL-22', "FJComp-Alexa Fluor 700-A<NA>"='CD69', "FJComp-BB515-A<NA>"='PD-1', "FJComp-BB630-A<NA>"='Gnly', "FJComp-BB700-A<NA>"="CD4", "FJComp-BB755-A<NA>" ="Perforin", "FJComp-BB790-A<NA>"="GzmB",
                              "FJComp-BUV395-A<NA>"="CD103", "FJComp-BUV496-A<NA>"="CD39", "FJComp-BUV563-A<NA>"="HLA-DR", "FJComp-BUV615-A<NA>"="CD27", "FJComp-BUV661-A<NA>"="CD127",  "FJComp-BUV737-A<NA>"="CD56", "FJComp-BUV805-A<NA>"="CD45", "FJComp-BV421-A<NA>"="CXCR5",
                              "FJComp-BV510-A<NA>"="Dump",  "FJComp-BV570-A<NA>"="CD8",  "FJComp-BV605-A<NA>"="IL17A", "FJComp-BV650-A<NA>"="CD3", "FJComp-BV711-A<NA>"="TNF", "FJComp-BV750-A<NA>"="CD62L", "FJComp-BV786-A<NA>"="IFNg", "FJComp-PE-A<NA>"="5OPRU-Tet",
                              "FJComp-PE-CF594-A<NA>"="IL-10", "FJComp-PE-Cy5-A<NA>"="CD161", "FJComp-PE-Cy7-A<NA>"="CXCR3")
colnames(fcs_protein_data)[match(names(fcs_protein_data_columns), colnames(fcs_protein_data))] <- fcs_protein_data_columns

# create CYT object (S4 object)
cyt_ihope <- createCYT(raw.data=fcs_protein_data, normalization.method='log')
cyt_ihope@meta.data$PseudoDonor  <- fcs_meta_data$Pseudo_Donor
cyt_ihope@meta.data$PseudoTissue <- fcs_meta_data$Pseudo_Tissue
cyt_ihope@meta.data$Tissue       <- case_match(cyt_ihope@meta.data$PseudoTissue, 1~'Blood', 2~'Spleen', 3~'Liver', 4~'Ileum', 5~'Caecum', 6~'Colon', 7~'mLN', 8~'Lung', 9~'LungLN', 10~'Skin')
cyt_ihope@meta.data$Tissue       <- factor(cyt_ihope@meta.data$Tissue, levels=c('Blood', 'Spleen', 'Liver', 'Ileum', 'Caecum', 'Colon', 'mLN', 'Lung', 'LungLN', 'Skin'))
tissues <-  levels(cyt_ihope@meta.data$Tissue)
rm(fcs_data, fcs_meta_data, fcs_protein_data, fcs_protein_data_columns, fcs_file)

# select relevant proteins for UMAP, remove markers used for MAIT cell identification
clustering_markers <- c('CD4', 'CD8','CD27', 'CD39', 'CD56', 'CD62L' ,'CD69', 'CD103', 'CD127', 'CXCR3', 'CXCR5', 'PD-1', 'HLA-DR', 'GzmB', 'Perforin', 'Gnly', 'IL-10', 'IL17A', 'IL-22', 'TNF', 'IFNg')
cyt_ihope          <- changeMarker(cyt_ihope, markers=clustering_markers)

# Downsample MAIT cells per donor and tissue
cyt <- lapply(tissues, function(i){
  cells  <- cyt_ihope@meta.data[cyt_ihope@meta.data$Tissue==i,]
  donors <- unique(cells$PseudoDonor)
  sampled_cells_df <- data.frame() # initate dataframe
  # Downsampling of 500 MAIT cells per tissue of each donor
  for (j in donors) {
    donor_cells     <- cells[cells$PseudoDonor==j, ]
    cells_to_sample <- min(nrow(donor_cells), 500)
    if (cells_to_sample>0) {
      set.seed(1337)
      sampled_cells_df <- rbind(sampled_cells_df, donor_cells[sample(nrow(donor_cells), cells_to_sample), ])
    }
  }
  # Downsampling to 2000 cells per tissue from all donors if too many
  if (nrow(sampled_cells_df)>2000) {
    set.seed(1337)
    sampled_cells_df <- sampled_cells_df[sample(nrow(sampled_cells_df), 2000), ]
  }
  sampled_cells_df
})

cyt <- do.call(rbind, cyt)
sampled_cells <- rownames(cyt) # vector of cells to keep
cyt <- subsetCYT(cyt_ihope, cells=sampled_cells)

### k-means clustering, UMAP
cyt_kmeans <- cyt
cyt_kmeans <- runCluster(cyt_kmeans, cluster.method='kmeans', verbose=TRUE, iter.max=100) # computes kmeans, flowsom, phenograph or other cluster.methods.
cyt_kmeans <- processingCluster(cyt_kmeans, k=20, perplexity=5, downsampling.size=1, verbose=TRUE) # no further downsampling
cyt_kmeans <- runUMAP(cyt_kmeans, n_neighbors=20, verbose=TRUE)
cyt_kmeans <- buildTree(cyt_kmeans, dim.type='umap', dim.use=1:2)

# loading cyt_kmeans.rds because of stochasticity in process of UMAP and tree building
cyt_kmeans <- readRDS('../Data/Fig2_UMAP_data.rds')

# Fig. 2A, Fig. S5A, # UMAP for each tissue
for (i in tissues){
  umapplot <- plot2D(cyt_kmeans, item.use=c('UMAP_1', 'UMAP_2'), color.by=paste0(i), alpha=ifelse(cyt_kmeans@meta.data$Tissue==i, yes=1, no=0.2), category='categorical', size=ifelse(cyt_kmeans@meta.data$Tissue==i, 1,0.2))+theme_transparent()+theme(plot.background=element_rect(fill='transparent'), panel.background=element_rect(fill='transparent'))+
  scale_color_manual(values=c(mycolors[i], 'other'='lightgrey'))+theme(legend.position='none')
  print(umapplot)
}

# Fig. 2B, Fig. S5A, # UMAP for each marker
for (i in clustering_markers){
  umapplot<- plot2D(cyt_kmeans, item.use=c('UMAP_1', 'UMAP_2'), color.by=i, alpha=0.9, category='numeric', size=1, main=paste0(i))+
  theme_minimal()+theme(legend.position='none')+scale_color_gradientn(colors=c("#4575b4", "#fafac8", "#b2182b"))
  print(umapplot)
}

# Fig.S5B, tree colored by branch
plotTree(cyt_kmeans, color.by='branch.id', show.node.name=TRUE, cex.size=1)+ scale_color_gradientn(colors=c('black','darkgrey','#A9261F', '#813c5e', 'orange' ))

# Fig. S5B, UMAP colored by Tree
plot2D(cyt_kmeans, item.use=c('UMAP_1', 'UMAP_2'), alpha=1, category='categorical', color.by='branch.id', size=1, show.cluser.id=FALSE)+theme_transparent()+ scale_color_manual(values=c('black','darkgrey','#A9261F', '#813c5e', 'orange' ))+theme(legend.position='none')

# Fig. S5B Tree for each marker
for (i in clustering_markers){
  Treeplot <- plotTree(cyt_kmeans, color.by=i, show.node.name=FALSE, cex.size=0.8, as.tree=FALSE)+
    scale_colour_gradientn(colors=c("#4575b4", "#fafac8", "#b2182b"))+
    theme(legend.position='none')+
    labs(title=paste0(i))
  print(Treeplot)
}



