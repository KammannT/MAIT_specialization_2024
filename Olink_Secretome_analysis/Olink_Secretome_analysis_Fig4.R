# Secretome analysis of purified cell-free supernatants obtained from purified immune cell populations by Olink Target-96 Inflammation Panel
# Setup
required_packages <- c('tidyverse', 'vtable', 'rstatix', 'scales', 'RColorBrewer', 'ggpubr', 'readxl',
                       'OlinkAnalyze', 'FactoMineR', 'factoextra')
for(package in required_packages){
  if(!require(package, character.only=TRUE, quietly=TRUE))
  {install.packages(package, dependencies=TRUE)
    library(package, character.only=TRUE) }}

script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)

olink_data <- read_csv(file='../Data/Fig4_Olink_data.csv')

# Remove IFNg in IL-12+IL-18 stimulated assays because of hook effect observed in pilot run with dilution serious of supernatant
olink_data     <- olink_data %>% filter(!(Assay=='IFNg' & Stimulation=='IL-12 + IL-18'))
# Remove IIL-12B and IL-18 in IL-12+IL-18 stimulated assays because of experimental design
olink_data     <- olink_data %>% filter(!(Assay=='IL-12B' & Stimulation=='IL-12 + IL-18')) %>% filter(!(Assay=='IL-18' & Stimulation=='IL-12 + IL-18'))

# Check for positive protein secretion in at least 3 out of 7 supernatants per cell population and condition
olink_data <- olink_data  %>% group_by(Donor, Assay, Stimulation, Celltype_CD56status) %>% mutate(aboveT=NPX>=LOD) %>% ungroup
olink_data <- olink_data  %>% group_by(Assay, Stimulation, Celltype_CD56status) %>% mutate(detected3x=sum(aboveT==TRUE)>=3) %>% mutate(detectednumber=sum(aboveT==TRUE)) %>% ungroup

# Number of different proteins robustly detected
olink_data %>% select(Celltype_CD56status, Assay, Stimulation, detected3x) %>% filter(detected3x==TRUE) %>% group_by(Stimulation, Celltype_CD56status) %>% summarize(Numberofproteins_detected3x=length(unique(Assay))) %>% print

# Fig. 4B
# Normalized Protein eXpression (N) Heatmap
olinkdataheatmap        <- olink_data %>% group_by(Assay, Stimulation, Celltype_CD56status) %>% summarize(MeanNPX=mean(NPX), detected3x=mean(detected3x) )
olinkdataheatmap_nulled_3xdetected       <- olinkdataheatmap %>% filter(detected3x==1) %>% mutate(MeanNPX=ifelse(test=MeanNPX<0, yes=0, no=MeanNPX))
olinkdataheatmap_nulled_3xdetected$Assay <- fct_reorder(olinkdataheatmap_nulled_3xdetected$Assay, olinkdataheatmap_nulled_3xdetected$MeanNPX, .fun=max, .desc=FALSE)


heatmapdata <- olinkdataheatmap %>% filter(Assay %in% olinkdataheatmap_nulled_3xdetected$Assay) %>% filter(Stimulation!='control') %>% droplevels()
ggplot(data=heatmapdata, aes(x=Celltype_CD56status, y=reorder(Assay, MeanNPX), fill=ifelse(detected3x==1, yes=MeanNPX, no=NA)))+
  geom_tile(color='black', lwd=0.2)+
  scale_fill_gradientn(limits=c(-4, 18), colors=c('#fee6ce', '#fec44f','#e31a1c','#a50f15', '#4d004b'), values=scales::rescale(c(-4,0,4,8,16)), na.value='white')+
  scale_x_discrete(expand=c(0,0))+
  labs(x='Celltype', y='Assay', fill='NPX')+
  facet_grid(.~Stimulation, scales='free', drop=TRUE)+
  theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
        axis.text.y=element_text(size=8),
        axis.ticks=element_blank(),
        panel.grid.major.y=element_blank(),
        axis.line.y=element_blank(),
        axis.title=element_text(face='plain', size=10))

# Fig. 4C
# CD56+ MAIT, resting vs IL-12 + IL-18
volcano_CD56pos_cyto_vs_base  <- olink_data %>% filter(Celltype_CD56status=='CD56+ MAIT' & Stimulation %in% c('control', 'IL-12 + IL-18')) %>% group_by(Donor, Assay) %>% droplevels %>%
  group_by(Assay, Stimulation) %>% mutate(detected3xinstimulation=sum(aboveT==TRUE)>=3 ) %>% ungroup %>% # detected3xincelltype = assays detected in  CD56- MAIT subset or NK cell subset, respectively, from 3 or more donors
  group_by(Assay) %>% filter(any(detected3xinstimulation==TRUE)) %>%
  filter(!Assay %in% c('IFNg', 'IL-12B', 'IL-18') ) # remove assays that are defect in IL-12 + IL-18 stimulation group

ttest_volcano_CD56pos_cyto_vs_base <- olink_ttest(df=volcano_CD56pos_cyto_vs_base, variable='Stimulation', pair_id='Donor') #
ttest_volcano_CD56pos_cyto_vs_base <- ttest_volcano_CD56pos_cyto_vs_base %>% mutate(logFC=-estimate) # positive logFC == enriched in IL12+IL18 stimulation
ttest_NPX_CD56pos_cyto_vs_base     <- ttest_volcano_CD56pos_cyto_vs_base %>% filter(p.value<0.05 ) %>% pull(OlinkID)

ttest_volcano_CD56pos_cyto_vs_base_table   <- ttest_volcano_CD56pos_cyto_vs_base
ttest_volcano_CD56pos_cyto_vs_base_table$Condition <- 'IL-12+IL-18 - resting'
ttest_volcano_CD56pos_cyto_vs_base_table   <- ttest_volcano_CD56pos_cyto_vs_base_table %>% select(Panel, method, OlinkID, UniProt, Assay, Condition, logFC, p.value, Adjusted_pval) %>% arrange(desc(logFC)) %>% View()

ttest_volcano_CD56pos_cyto_vs_base %>% ggplot(aes(x=logFC, y=-log10(p.value)))+
  geom_hline(yintercept=-log10(0.05), linetype="dotted", color='grey')+
  ggrepel::geom_label_repel(data=subset(ttest_volcano_CD56pos_cyto_vs_base, OlinkID %in% ttest_NPX_CD56pos_cyto_vs_base ), aes(label=Assay), box.padding=0.1, show.legend=FALSE, size=2)+
  geom_point(aes(fill=logFC), shape=21, size=2.5, show.legend=FALSE)+
  scale_x_continuous(limits=c(-2,8), breaks=c(-2,0,2,4,6,8), guide='prism_offset', expand=c(0,0) )+
  scale_y_continuous(limits=c(0,8), breaks=c(0,1,2,3,4,5,6,7,8), guide='prism_offset', expand=c(0.02,0) )+
  scale_fill_gradientn(colors=c('#4d4d4d', 'white', '#b2182b'), na.value='transparent', breaks=c(-2,8), limits=c(-2,8), labels=c('control', 'IL-12+\nIL-18'), values=scales::rescale(c(-2,0,8)) )+
  theme(panel.grid.major.y=element_blank())+
  labs(x='Log2FC (IL-12 + IL-18 / control)', y='-log10(p-value)', fill='', title='CD56+ MAIT cell secretome')

# Fig. 4D
# CD56+ MAIT, resting vs PMA+Ionomycin
volcano_CD56pos_PMAI_vs_base  <- olink_data %>% filter(Celltype_CD56status=='CD56+ MAIT' & Stimulation %in% c('control', 'PMA + Ionomycin')) %>% group_by(Donor, Assay) %>% droplevels %>%
  group_by(Assay, Stimulation) %>% mutate(detected3xinstimulation=sum(aboveT==TRUE)>=3 ) %>% ungroup %>% # detected3xincelltype = assays detected in  CD56- MAIT subset or NK cell subset, respectively, from 3 or more donors
  group_by(Assay) %>% filter(any(detected3xinstimulation==TRUE))

ttest_volcano_CD56pos_PMAI_vs_base <- olink_ttest(df=volcano_CD56pos_PMAI_vs_base, variable='Stimulation', pair_id='Donor') #
ttest_volcano_CD56pos_PMAI_vs_base <- ttest_volcano_CD56pos_PMAI_vs_base %>% mutate(logFC=-estimate) # positive logFC == enriched in PMA+Ionomycin stimulation
ttest_NPX_CD56pos_PMAI_vs_base     <- ttest_volcano_CD56pos_PMAI_vs_base %>% filter(p.value<0.05) %>% pull(OlinkID)

ttest_volcano_CD56pos_PMAI_vs_base_table   <- ttest_volcano_CD56pos_PMAI_vs_base
ttest_volcano_CD56pos_PMAI_vs_base_table$Condition <- 'PMA+Ionomycin - resting'
ttest_volcano_CD56pos_PMAI_vs_base_table   <- ttest_volcano_CD56pos_PMAI_vs_base_table %>% select(Panel, method, OlinkID, UniProt, Assay, Condition, logFC, p.value, Adjusted_pval) %>% arrange(desc(logFC)) %>% View()

ttest_volcano_CD56pos_PMAI_vs_base %>% ggplot(aes(x=logFC, y=-log10(p.value)))+
  geom_hline(yintercept=-log10(0.05), linetype="dotted", color='grey')+
  ggrepel::geom_label_repel(data=subset(ttest_volcano_CD56pos_PMAI_vs_base, OlinkID %in% ttest_NPX_CD56pos_PMAI_vs_base), aes(label=Assay), box.padding=0.1, size=2, show.legend=FALSE)+
  geom_point(aes(fill=logFC), shape=21, size=2.5, show.legend=FALSE)+
  scale_x_continuous(limits=c(-2,10), breaks=c(-2,0,2,4,6,8,10),guide='prism_offset', expand=c(0,0) )+
  scale_y_continuous(limits=c(0, 7), breaks=c(0,1,2,3,4,5,6,7), guide='prism_offset', expand=c(0.02,0) )+
  scale_fill_gradientn(colors=c('#4d4d4d', 'white', '#2166ac'), na.value='transparent', breaks=c(-2,10), limits=c(-2,10), labels=c('base', 'PMA+\nIonomycin'), values=scales::rescale(c(-2,0,10)) )+
  theme(panel.grid.major.y=element_blank())+
  labs(x='Log2FC (PMA + Ionomycin / control)', y='-log10(p-value)', fill='', title='CD56+ MAIT cell secretome')

# Fig. 4E
# PCA of CD56+ MAIT vs NK cells
pca_data_cyto_mait_nk <- olink_data %>% filter(Stimulation=='IL-12 + IL-18', Celltype_CD56status=='CD56+ MAIT' | Celltype_CD56status=='NK') %>%
  select(-CD56status, -Normalization, -Adj_factor, -Celltype, -detected3x, -detectednumber) %>%
  group_by(Assay, Celltype_CD56status) %>% mutate(detected3xincelltype=sum(aboveT==TRUE)>=3 ) %>% ungroup # detected3xincelltype = assays detected in  CD56 MAIT subset or NK cell subset, respectively, from 3 or more donors

pca_data_cyto_mait_nk_imputed <- pca_data_cyto_mait_nk %>% mutate(NPX=ifelse(NPX<LOD, yes=LOD, no=NPX)) %>% # impute missing data if NPX below LOD with LOD
  group_by(Assay) %>% filter(any(detected3xincelltype==TRUE)) %>% ungroup %>% droplevels()

pca_data_cyto_mait_nk <- pca_data_cyto_mait_nk_imputed %>% mutate(Sample=paste0(Donor,' ',Celltype_CD56status)) %>%
  select(Sample, Celltype_CD56status, Assay, NPX) %>% pivot_wider(names_from='Assay', values_from='NPX') %>%
  column_to_rownames(var='Sample')

pca_data_cyto_mait_nk <- FactoMineR::PCA(pca_data_cyto_mait_nk, graph=FALSE, scale.unit=TRUE, quali.sup='Celltype_CD56status')

fviz_screeplot(pca_data_cyto_mait_nk, choice='variance', addlabels=TRUE, ylim=c(0, 50), ggtheme=theme_tk, barcolor='black', barfill='darkgrey')+labs(title='', x='PC', y='Explained variance (%)')+scale_y_continuous(expand=c(0,0))

x <- fviz_pca_biplot(pca_data_cyto_mait_nk, geom.ind='point', geom.var='', title='')
y <- x + geom_point(aes(fill=pca_data_cyto_mait_nk$call$X$Celltype_CD56status), color='black', shape=21, size=2, show.legend=FALSE)+
  labs(fill='Celltype',)+
  scale_fill_manual(values=celltype_colors)+
  scale_x_continuous(limits=c(-9,9), breaks=c(-9,0,9), expand=c(0,0), guide='prism_offset')+
  scale_y_continuous(limits=c(-8,8), breaks=c(-8,0,8), expand=c(0,0), guide='prism_offset')+
  theme_tk+
  theme(panel.grid.major.y=element_blank())
print(y)

# PCA of CD56+ MAIT vs CD45RO+ CD8+ T cells
pca_data_pmai_mait_cd8t <- olink_data %>% filter(Stimulation=='PMA + Ionomycin', Celltype_CD56status=='CD56+ MAIT' | Celltype_CD56status=='CD8T') %>%
  select(-CD56status, -Normalization, -Adj_factor, -Celltype, -detected3x, -detectednumber) %>%
  group_by(Assay, Celltype_CD56status) %>% mutate(detected3xincelltype=sum(aboveT==TRUE)>=3 ) %>% ungroup # detected3xincelltype = assays detected in  CD56 MAIT subset or NK cell subset, respectively, from 3 or more donors

pca_data_pmai_mait_cd8t_imputed <- pca_data_pmai_mait_cd8t %>% mutate(NPX=ifelse(NPX<LOD, yes=LOD, no=NPX)) %>% # impute missing data if NPX below LOD with LOD
  group_by(Assay) %>% filter(any(detected3xincelltype==TRUE)) %>% ungroup %>% droplevels()

pca_data_pmai_mait_cd8t <- pca_data_pmai_mait_cd8t_imputed %>% mutate(Sample=paste0(Donor,' ',Celltype_CD56status)) %>%
  select(Sample, Celltype_CD56status, Assay, NPX) %>% pivot_wider(names_from='Assay', values_from='NPX') %>%
  column_to_rownames(var='Sample')

pca_data_pmai_mait_cd8t  <- FactoMineR::PCA(pca_data_pmai_mait_cd8t, graph=FALSE, scale.unit=TRUE, quali.sup='Celltype_CD56status')

fviz_screeplot(pca_data_pmai_mait_cd8t, choice='variance', addlabels=TRUE, ylim=c(0, 50), ggtheme=theme_tk, barcolor='black', barfill='darkgrey')+labs(title='', x='PC', y='Explained variance (%)')+scale_y_continuous(expand=c(0,0))

x <- fviz_pca_biplot(pca_data_pmai_mait_cd8t, geom.ind='point', geom.var='', title='')
y <- x + geom_point(aes(fill=pca_data_pmai_mait_cd8t$call$X$Celltype_CD56status), color='black', shape=21, size=2, show.legend=FALSE)+
  labs(fill='Celltype', color='Contribution')+
  scale_fill_manual(values=celltype_colors)+
  scale_x_continuous(limits=c(-8,8), breaks=c(-8,0,8), expand=c(0,0), guide='prism_offset')+
  scale_y_continuous(limits=c(-9,9), breaks=c(-9,0,9), expand=c(0,0))+
  theme_tk+
  theme(panel.grid.major.y=element_blank())
print(y)

# Fig. 4F:
volcano_NK_vs_CD56pos_cyto <- olink_data %>% filter(Celltype_CD56status %in% c('CD56+ MAIT', 'NK') & Stimulation=='IL-12 + IL-18') %>% group_by(Donor, Assay) %>% droplevels %>%
  group_by(Assay, Celltype_CD56status) %>% mutate(detected3xincelltype=sum(aboveT==TRUE)>=3 ) %>% ungroup %>% # detected3xincelltype = assays detected in  CD56- MAIT subset or NK cell subset, respectively, from 3 or more donors
  group_by(Assay) %>% filter(any(detected3xincelltype==TRUE))

ttest_volcano_NK_vs_CD56pos_cyto<- olink_ttest(df=volcano_NK_vs_CD56pos_cyto, variable='Celltype_CD56status', pair_id='Donor')
ttest_volcano_NK_vs_CD56pos_cyto<- ttest_volcano_NK_vs_CD56pos_cyto %>% mutate(logFC=estimate) # enriched in CD56+ MAIT = positive logFC
ttest_NPX_cyto_NK_vs_CD56pos    <- ttest_volcano_NK_vs_CD56pos_cyto %>% filter(Assay %in% c('CCL20', 'CD5', 'CSF-1', 'VEGFA' ,'IL-6', 'IL-8', 'TRAIL', 'TRANCE', 'TWEAK', 'OSM') ) %>% pull(OlinkID)
ttest_volcano_NK_vs_CD56pos_cyto_table   <- ttest_volcano_NK_vs_CD56pos_cyto
ttest_volcano_NK_vs_CD56pos_cyto_table$Condition <- 'IL-12 + IL-18'
ttest_volcano_NK_vs_CD56pos_cyto_table <- ttest_volcano_NK_vs_CD56pos_cyto_table %>% select(Panel, method, OlinkID, UniProt, Assay, Condition, logFC, p.value, Adjusted_pval) %>% arrange(desc(logFC)) %>% View()


ttest_volcano_NK_vs_CD56pos_cyto %>% ggplot(aes(x=logFC, y=-log10(p.value)))+
  geom_hline(yintercept=-log10(0.05), linetype="dotted", color='grey')+
  ggrepel::geom_label_repel(data=subset(ttest_volcano_NK_vs_CD56pos_cyto, OlinkID %in% ttest_NPX_cyto_NK_vs_CD56pos), aes(label=Assay), box.padding=0.1, size=2, min.segment.length=0.2, show.legend=FALSE)+
  geom_point(aes(fill=logFC), shape=21, size=2.5, show.legend=FALSE)+
  scale_x_continuous(limits=c(-8.2,8.2), breaks=c(-8,-4,0,4,8), guide='prism_offset', expand=c(0,0) )+
  scale_y_continuous(limits=c(0,3), guide='prism_offset', expand=c(0.02,0.02) )+
  scale_fill_gradientn(colors=c(celltype_colors[3], 'white', celltype_colors[2]), na.value='transparent', breaks=c(-8,8), limits=c(-8,8), labels=c('NK', 'CD56+\nMAIT') )+
  theme(panel.grid.major.y=element_blank())+
  labs(x='Log2FC (CD56+ MAIT / NK)', y='-log10(p-value)', fill='', title='CD56+ MAIT vs NK cell secretome')

# Fig. 4G:
# CD56+ MAIT vs CD8T cells after PMA + Ionomycin stimulation
volcano_CD8T_vs_CD56pos_PMAI <- olink_data %>% filter(Celltype_CD56status %in% c('CD56+ MAIT', 'CD8T') & Stimulation=='PMA + Ionomycin') %>% group_by(Donor, Assay) %>% droplevels %>%
  group_by(Assay, Celltype_CD56status) %>% mutate(detected3xincelltype=sum(aboveT==TRUE)>=3 ) %>% ungroup %>% # detected3xincelltype = assays detected in  CD56- MAIT subset or NK cell subset, respectively, from 3 or more donors
  group_by(Assay) %>% filter(any(detected3xincelltype==TRUE))
ttest_volcano_CD8T_vs_CD56pos_PMAI<- olink_ttest(df=volcano_CD8T_vs_CD56pos_PMAI, variable='Celltype_CD56status', pair_id='Donor')
ttest_volcano_CD8T_vs_CD56pos_PMAI<- ttest_volcano_CD8T_vs_CD56pos_PMAI %>% mutate(logFC=estimate) # enriched in CD56- MAIT = positive logFC
ttest_NPX_CD8T_vs_CD56pos_PMAI    <- ttest_volcano_CD8T_vs_CD56pos_PMAI %>% filter(Assay %in% c('CCL20', 'IL-17A', 'VEGFA', 'TNFRSF9', 'TNFSF14', 'IL-18R1', 'IL-8', 'uPA', 'TWEAK') ) %>% pull(OlinkID)

ttest_volcano_CD8T_vs_CD56pos_PMAI_table   <- ttest_volcano_CD8T_vs_CD56pos_PMAI
ttest_volcano_CD8T_vs_CD56pos_PMAI_table$Condition <- 'PMA + Ionomycin'
ttest_volcano_CD8T_vs_CD56pos_PMAI_table <- ttest_volcano_CD8T_vs_CD56pos_PMAI_table %>% select(Panel, method, OlinkID, UniProt, Assay, Condition, logFC, p.value, Adjusted_pval) %>% arrange(desc(logFC)) %>% View()

ttest_volcano_CD8T_vs_CD56pos_PMAI %>% ggplot(aes(x=logFC, y=-log10(p.value)))+
  geom_hline(yintercept=-log10(0.05), linetype="dotted", color='grey')+
  ggrepel::geom_label_repel(data=subset(ttest_volcano_CD8T_vs_CD56pos_PMAI, OlinkID %in% ttest_NPX_CD8T_vs_CD56pos_PMAI), aes(label=Assay), box.padding=0.1, size=2, min.segment.length=0.4, show.legend=FALSE)+
  geom_point(aes(fill=logFC), shape=21, size=2.5, show.legend=FALSE)+
  scale_x_continuous(limits=c(-8,8), breaks=c(-8,-4,0,4,8), guide='prism_offset', expand=c(0,0) )+
  scale_y_continuous(limits=c(0, 7), breaks=c(0,1,2,3,4,5,6,7), guide='prism_offset', expand=c(0.02,0) )+
  scale_fill_gradientn(colors=c(celltype_colors[4], 'white', celltype_colors[2]), na.value='transparent', breaks=c(-7,7), limits=c(-7,7), labels=c('CD8 T', 'CD56+\nMAIT'))+
  theme(panel.grid.major.y=element_blank())+
  labs(x='Log2FC (CD56+ MAIT / CD8 T)', y='-log10(p-value)', fill='', title='CD56+ MAIT vs CD8 T cell secretome')

