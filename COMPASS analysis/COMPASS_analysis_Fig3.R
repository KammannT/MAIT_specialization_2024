# COMPASS analysis
# COMPASS was developed by Lin et al. (2015), Nat. Biotechnol, DOI: 10.1038/nbt.3187
# Setup
required_packages1 <- c("readxl", "reshape2", "see", "purrr", "ggprism", "rstatix", "scales", "stats", "RColorBrewer", "Hmisc", "corrplot", "vtable", "ggpubr", "tidyverse", "dplyr")
for(package in required_packages1){if(!require(package, character.only=TRUE))
{install.packages(package, dependencies=TRUE)
  library(package, character.only=TRUE) }}
if (!require('BiocManager', quietly=TRUE)) install.packages('BiocManager')
if (!require('COMPASS', quietly=TRUE)) devtools::install_github('RGLab/COMPASS')
if (!require('cytoqc', quietly=TRUE)) devtools::install_github('RGLab/cytoqc')
required_packages2 <- c('flowWorkspace', 'flowWorkspaceData', 'CytoML')
lapply(required_packages2, FUN=function(x) {BiocManager::install(x) })
libraries <- c('COMPASS', 'cytoqc', required_packages1, required_packages2)
lapply(libraries, FUN=function(x) {library(x, character.only=TRUE) })
mycolors <- c('Blood'='#E88984', 'Spleen'='#A9261F', 'Liver'='#813c5e', 'Ileum'='#8c510a', 'Caecum'='#bf812d', 'Colon'='#dfc27d', 'mLN'='#f6e8c3', 'Lung'='#c7eae5', 'LungLN'='#92c5de', 'Skin'='#01662c')

script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)

# Fig. 3H
## load Count data
cd56compass             <- read_csv('../Data/Fig3_COMPASS_data.csv')
cd56compass$Stimulation <- recode_factor(cd56compass$Stimulation, none='unstimulated', PMA_I='PMA + Ionomycin')
cd56compass             <- cd56compass %>% select(-CountCD56)

## run COMPASS
metadata_PI <- cd56compass %>% filter(Stimulation %in% c('unstimulated', 'PMA + Ionomycin'))  %>% select(1:5) %>% mutate(ID=paste(`Sample:`, Tissue, CD56status, sep='-')) %>% as.data.frame()
count_base  <- cd56compass %>% filter(Tissue=='Liver') %>% filter(Stimulation=='unstimulated')    %>% arrange(CD56status, Donor) %>% select(-(1:5)) %>% as.matrix()
count_PI    <- cd56compass %>% filter(Tissue=='Liver') %>% filter(Stimulation=='PMA + Ionomycin') %>% arrange(CD56status, Donor) %>% select(-(1:5)) %>% as.matrix()

rownames(count_base)      <- subset(metadata_PI$ID, metadata_PI$Stimulation=='unstimulated')
rownames(count_PI)        <- subset(metadata_PI$ID, metadata_PI$Stimulation=='PMA + Ionomycin')
metadata_PI               <- subset(metadata_PI,    metadata_PI$Stimulation=="PMA + Ionomycin")

## Change name of gates into COMPASS-suitable format.
nms <- basename(colnames(count_PI))
nms <- gsub(nms, pattern=' | ', replacement=',', fixed=TRUE) # Flowjo V10 adds " | " instead of ","
nms <- COMPASS::translate_marker_names(nms)
colnames(count_PI) <- nms

nms <- basename(colnames(count_base))
nms <- gsub(nms, pattern=' | ', replacement=',', fixed=TRUE) # Flowjo V10 adds " | " instead of ","
nms <- COMPASS::translate_marker_names(nms)
colnames(count_base) <- nms
rm(nms)

## COMPASS modelling
# calculating model is computationally expensive
compass_liver_fit_PI <- COMPASS::SimpleCOMPASS(n_s=count_PI, n_u=count_base,
                                               meta=metadata_PI,
                                               individual_id='ID',
                                               iterations=40000,
                                               replications=8, seed=1337)
compass_liver_fit_PI_scores  <- COMPASS::scores(compass_liver_fit_PI)

## COMPASS heatmap
plot(compass_liver_fit_PI, row_annotation='CD56status',
     border_color='black',
     legend_breaks=c(0, 0.2,0.4,0.6,0.8),
     palette=colorRampPalette(brewer.pal(10, "Purples"))(20),
     cellwidth=4, cellheight=4)

# Fig. 3I
richness_data  <- cd56compass %>% filter(Stimulation=='PMA + Ionomycin') %>% select(-`Sample:`) %>% droplevels
n_combinations <- (dim(richness_data)[-1]-4) # number of effector molecule combinations

richness1 <- richness_data %>% mutate(Richness=rowSums(.[, 5:260] !=0, na.rm=TRUE)) %>% mutate(totalcells=rowSums(.[, 5:260], na.rm=TRUE)) %>% select(1:4, Richness, totalcells)

sign_test_richness <- richness1 %>% arrange(Donor, CD56status) %>% wilcox_test(Richness~CD56status, detailed=TRUE, paired=TRUE) %>%
  add_significance(p.col='p', cutpoints=c(0,0.001,0.01,0.05,1), symbols=c('***','**','*','ns')) %>% add_xy_position(x='CD56status', fun='max')
ggplot(richness1, aes(y=Richness, x=CD56status))+
  geom_boxplot(aes(fill=CD56status), color='black', show.legend=FALSE, outlier.shape=NA)+
  geom_line(aes(group=interaction(Donor)), position=position_jitter(width=0.1, height=0, seed=42), color='black', linewidth=0.5)+
  geom_point(fill='white', position=position_jitter(width=0.1, height=0, seed=42), shape=21, size=2, show.legend=FALSE)+
  stat_summary(fun=mean, geom='point', shape=15, size=1, color='black', show.legend=FALSE)+
  scale_fill_manual(values=c('#813c5e', '#813c5e'))+
  add_pvalue(sign_test_richness, tip.length=0, label.size=3, bracket.size=0.35, y.position=80*1.02)+
  scale_y_continuous(limits=c(0,80*1.05), breaks=c(0,20,40,60,80), expand=c(0.02,0), guide='prism_offset')+
  scale_x_discrete(labels=c('-', '+'))+
  labs(x='', y='Richness (# of effector combinations)')+
  theme_classic()+
  theme(axis.text.x=element_text(size=10))

# Fig. 3J
sign_test_compass_CD56_FS <- compass_liver_fit_PI_scores %>% group_by(Tissue) %>% arrange(Donor, Tissue) %>% wilcox_test(FS~CD56status, detailed=TRUE, paired=TRUE, p.adjust.method='none') %>%
  add_significance(p.col='p', cutpoints=c(0,0.001,0.01,0.05,1), symbols=c('***','**','*','ns')) %>%
  filter(p<0.05) %>% add_xy_position(x='CD56status', group='Tissue', fun='max')

ggplot(compass_liver_fit_PI_scores, aes(y=FS, x=CD56status))+
  geom_boxplot(aes(fill=Tissue), color='black', show.legend=FALSE, outlier.shape=NA)+
  geom_line(aes(group=interaction(Tissue, Donor)), position=position_jitter(width=0.1, height=0, seed=42), color='black', linewidth=0.5)+
  geom_point(fill='white', shape=21, size=2, show.legend=FALSE, position=position_jitter(width=0.1, height=0, seed=42))+
  stat_summary(fun=mean, geom='point', shape=15, size=1, color='black', show.legend=FALSE)+
  scale_fill_manual(values=mycolors)+
  add_pvalue(sign_test_compass_CD56_FS, tip.length=0, label.size=3, bracket.size=0.35, y.position=0.08*1.02)+
  scale_y_continuous(limits=c(0,0.08*1.05), breaks=c(0,0.02,0.04,0.06,0.08), expand=c(0.02,0), guide='prism_offset')+
  scale_x_discrete(labels=c('-', '+'))+
  labs(x='', y='Functionality Score')+
  theme_classic()