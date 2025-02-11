# Data is available in the GEO database under accession number GSE274338.
# Short summary of experimentation:
# PBMCs were isolated from four healthy human blood donors.
# PBMCs were stained with 5OPRU-hMR1-Tet and positively enriched before sorting of MAIT and non-MAIT cells.
# Libraries were generated following the 10X Genomics Next-GEM v2 5' Kit with Feature Barcode Technology and TCR sequencing.

# Setup -------------------------------------------------------------------

devtools::install:github('ncborcherding/scRepertoire') #  DOI: 10.12688/f1000research.22139.2, https://www.borch.dev/uploads/screpertoire/
library(scRepertoire)
required_packages <- c("readxl", "reshape2", "see", "purrr", "ggprism", "rstatix", "viridis",
                       "scales", "stats", "RColorBrewer", "ggpubr", "Seurat", "HGNChelper", "tidyverse")
for(package in required_packages){if(!require(package, character.only=TRUE))
{install.packages(package, dependencies=TRUE)
  library(package, character.only=TRUE) }}
rm(required_packages)

mycolors2 <- c('CD56neg'='#C0C0C0','CD56pos'='#008080')

script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)

# Import TCRseq data -----------------------------------------------------------

# load CD56-annotated MAIT cell seurat object
cd56maits <- readRDS('cd56maits.rds')
cd56maits$Donor_CD56status <- paste0(cd56maits$donor, '_', cd56maits$CD56status)

# load VDJ annotations from cellranger output
VDJ1 <- read.csv('../D1_multi_2023_06_13/vdj_t/filtered_contig_annotations.csv')
VDJ2 <- read.csv('../D2_multi_2023_06_13/vdj_t/filtered_contig_annotations.csv')
VDJ3 <- read.csv('../D3_multi_2023_06_13/vdj_t/filtered_contig_annotations.csv')
VDJ4 <- read.csv('../D4_multi_2023_06_13/vdj_t/filtered_contig_annotations.csv')

# subset VDJ information to match MAIT cell seurat object
mait_ids_D1 <- colnames(cd56maits)[cd56maits$donor=='D1']
mait_ids_D1 <- sub(".*_", "", mait_ids_D1)
mait_ids_D2 <- colnames(cd56maits)[cd56maits$donor=='D2']
mait_ids_D2 <- sub(".*_", "", mait_ids_D2)
mait_ids_D3 <- colnames(cd56maits)[cd56maits$donor=='D3']
mait_ids_D3 <- sub(".*_", "", mait_ids_D3)
mait_ids_D4 <- colnames(cd56maits)[cd56maits$donor=='D4']
mait_ids_D4 <- sub(".*_", "", mait_ids_D4)

VDJ1mait <- VDJ1 %>% filter(barcode %in% mait_ids_D1)
VDJ2mait <- VDJ2 %>% filter(barcode %in% mait_ids_D2)
VDJ3mait <- VDJ3 %>% filter(barcode %in% mait_ids_D3)
VDJ4mait <- VDJ4 %>% filter(barcode %in% mait_ids_D4)

# TCR object
contig_list <- list(VDJ1mait, VDJ2mait, VDJ3mait, VDJ4mait)
TCR <- combineTCR(contig_list, samples=c('D1', 'D2', 'D3', 'D4'), removeNA=TRUE, filterMulti=TRUE)
TCR <- addVariable(TCR, name='donor', variables=c('D1', 'D2', 'D3', 'D4'))

# Data deposition
TCR_file <- do.call(rbind.data.frame, TCR)
TCR_file <- TCR_file %>% dplyr::select(-sample) %>% rename(CellBarcode=barcode, Donor=donor, TCRa=TCR1, TCRb=TCR2)
saveRDS(TCR_file, file='./Data deposition/TCR_MAIT_Kammann_et_al.rds')

# MAIT cell TCR clonotype annotation ----------------------------------------

# paste donor metadata to cell ID
donors <- c('D1', 'D2', 'D3', 'D4')
for (i in donors){
  TCR_info <- TCR[[i]]
  cd56maits_info  <- cd56maits@meta.data %>% filter(donor==i) %>% rownames_to_column(var='barcode') %>% filter(barcode %in% TCR_info$barcode)
  TCR_info$CD56status <- cd56maits_info$CD56status
  TCR_info$Donor_CD56status <- cd56maits_info$Donor_CD56status
  TCR_info$gender <- cd56maits_info$gender
  TCR[[i]] <- TCR_info
  rm(TCR_info, cd56maits_info)
}
rm(donors)

# combine TCR information and seurat object
fullTCR            <- combineExpression(TCR, cd56maits, cloneCall='CTstrict', chain='both', group.by='CD56status', proportion=TRUE, filterNA=FALSE, cloneSize=c(Small=0.001, Medium=0.01, Large=0.1, Hyperexpanded=1), addLabel=FALSE)
fullTCR_cd56       <- scRepertoire:::.expression2List(fullTCR, split.by='CD56status')
fullTCR_cd56_donor <- scRepertoire:::.expression2List(fullTCR, split.by='Donor_CD56status')

# TRAV and TRBV usage ---------------------------------------------------

# Fig. S9C
# TRAV comparison
TRAV_cd56_plot      <- vizGenes(fullTCR_cd56, x.axis='TRAV', plot='barplot', group.by='Donor_CD56status', order='variance', scale=TRUE)+scale_y_continuous(limits=c(0,1), breaks=c(0,0.25,0.50,0.75,1.0), expand=c(0,0), guide='prism_offset')
TRAV_cd56_plot_data <- TRAV_cd56_plot[["data"]] %>% arrange(desc(mean))
TRAV_chains         <- TRAV_cd56_plot_data %>% distinct(x.axis) %>% head(5) %>% pull(x.axis) %>% droplevels()
top10_TRAV          <- TRAV_cd56_plot_data %>% separate(y.axis, into=c('Donor', 'CD56status'), sep='_', remove=FALSE) %>% mutate(TRAV=x.axis, Donor_CD56status=y.axis, Abundance=mean*100) %>% filter(x.axis %in% TRAV_chains)

ggplot(top10_TRAV, aes(x=CD56status, y=Abundance, fill=CD56status))+
  geom_boxplot(position=position_dodge(width=0.8), color='black', show.legend=FALSE)+
  geom_line(aes(group=Donor), position='identity')+
  geom_point(aes(group=CD56status), shape=21, size=1, fill='white', position=position_dodge(width=0.8))+
  scale_fill_manual(values=mycolors2, labels=c('CD56-','CD56+'))+
  scale_x_discrete()+
  scale_y_continuous(limits=c(0,100), breaks=c(0,20,40,60,80,100), expand=c(0.05,0), guide='prism_offset')+
  labs(x='', y='chain abundance (%)')+
  theme_classic()+
  theme(axis.text.x=element_blank(), panel.grid.major.y=element_line(color='grey',linewidth=0.35, linetype='dotted'), strip.background=element_blank(), strip.placement='outside')+
  facet_grid(.~TRAV, switch='x')

# Fig. S9C
# TRBV comparison
TRBV_cd56_plot      <- vizGenes(fullTCR_cd56, x.axis='TRBV', plot='barplot', group.by='Donor_CD56status', order='variance', scale=TRUE)+scale_y_continuous(limits=c(0,1), breaks=c(0,0.25,0.50,0.75,1.0), expand=c(0,0), guide='prism_offset')
TRBV_cd56_plot_data <- write_delim(as.data.frame(TRBV_cd56_plot[["data"]]), file='./Tables/MAIT_fullTCR_cd56_TRBV_combined.csv')
TRBV_cd56_plot_data <- TRBV_cd56_plot[["data"]] %>% arrange(desc(mean))
TRBV_chains         <- TRBV_cd56_plot_data %>% distinct(x.axis) %>% head(5) %>% pull(x.axis) %>% droplevels()
top10_TRBV          <- TRBV_cd56_plot_data %>% separate(y.axis, into=c('Donor', 'CD56status'), sep='_', remove=FALSE) %>% mutate(TRBV=x.axis, Donor_CD56status=y.axis, Abundance=mean*100) %>% filter(x.axis %in% TRBV_chains)

ggplot(top10_TRBV, aes(x=CD56status, y=Abundance, fill=CD56status))+
  geom_boxplot(position=position_dodge(width=0.8), color='black', show.legend=FALSE)+
  geom_line(aes(group=Donor), position='identity')+
  geom_point(aes(group=CD56status), shape=21, size=2, fill='white', position=position_dodge(width=0.8))+
  scale_fill_manual(values=mycolors2, labels=c('CD56-','CD56+'))+
  scale_x_discrete()+
  scale_y_continuous(limits=c(0,50), breaks=c(0,10,20,30,40,50), expand=c(0.05,0), guide='prism_offset')+
  labs(x='', y='chain abundance (%)')+
  theme_classic()+
  theme(axis.text.x=element_blank(), panel.grid.major.y=element_line(color='grey',linewidth=0.35, linetype='dotted'), strip.background=element_blank(), strip.placement='outside')+
  facet_grid(.~TRBV, switch='x')

# VDJ clonotype comparison (only TRAV1-2+ MAIT cells) ----------------------------------------------------

# filter for TRAV1-2+ cells
TRAV12_cd56pos <- fullTCR_cd56[['CD56pos']] %>% separate(CTgene, into='MAIT_TRAV', sep='\\.', remove=FALSE) %>% filter(MAIT_TRAV=='TRAV1-2')
TRAV12_cd56neg <- fullTCR_cd56[['CD56neg']] %>% separate(CTgene, into='MAIT_TRAV', sep='\\.', remove=FALSE) %>% filter(MAIT_TRAV=='TRAV1-2')
TRAV12 <- rbind(TRAV12_cd56pos, TRAV12_cd56neg)

# check size of (smaller) CD56+ MAIT cell population in each donor for sampling
table(TRAV12_cd56pos$Donor_CD56status)
table(TRAV12_cd56neg$Donor_CD56status)

# Downsampling of CD56-MAIT cells to donor-specific equal number of CD56+MAIT cells for diversity comparison.
set.seed(42)
sampled_D1 <- TRAV12_cd56neg %>% filter(Donor_CD56status=='D1_CD56neg') %>% sample_n(159, replace=FALSE, seed=42)
sampled_D2 <- TRAV12_cd56neg %>% filter(Donor_CD56status=='D2_CD56neg') %>% sample_n(523, replace=FALSE, seed=42)
sampled_D3 <- TRAV12_cd56neg %>% filter(Donor_CD56status=='D3_CD56neg') %>% sample_n(562, replace=FALSE, seed=42)
sampled_D4 <- TRAV12_cd56neg %>% filter(Donor_CD56status=='D4_CD56neg') %>% sample_n(92,  replace=FALSE, seed=42)
TRAV12_cd56neg_sampled <- rbind(sampled_D1, sampled_D2, sampled_D3, sampled_D4)

list_TRAV12_sampled <- list(TRAV12_cd56pos, TRAV12_cd56neg_sampled)
TRAV12_sampled      <- rbind(TRAV12_cd56pos, TRAV12_cd56neg_sampled)
rm(sampled_D1, sampled_D2, sampled_D3, sampled_D4, TRAV12_cd56neg_sampled, TRAV12_cd56neg, TRAV12_cd56pos)

# Fig. 5G
# TCR repertoire diversity
TCR_diversity_plot <- clonalDiversity(list_TRAV12_sampled, cloneCall='CTstrict', group.by='Donor_CD56status', x.axis='donor', n.boots=500)
TCR_diversity_plot
TCR_diversity_plot_data <- TCR_diversity_plot[['data']] %>% filter(variable=='inv.simpson') %>% pivot_wider(names_from='variable', values_from='value') %>% separate(Donor_CD56status, into=c('Donor', 'CD56status'))

ttest <- TCR_diversity_plot_data %>% arrange(donor) %>% t_test(inv.simpson~CD56status, paired=TRUE) %>% add_significance(p.col='p', cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c('***', '**', '*', 'ns')) %>% add_xy_position(x='CD56status', fun='max')
ggplot(TCR_diversity_plot_data, aes(x=CD56status, y=inv.simpson))+
  geom_boxplot(aes(fill=CD56status), color='black', show.legend=FALSE)+
  scale_fill_manual(values=mycolors2)+
  geom_line(aes(group=donor), linewidth=0.35, color='black')+
  geom_point(size=2.5, color='black', shape=21, fill='white', show.legend=FALSE)+
  add_pvalue(ttest, tip.length=0, bracket.size=0.35, label.size=5, y.position=30)+
  scale_x_discrete(labels=c('CD56-', 'CD56+'))+
  scale_y_continuous(limits=c(0, 33), expand=c(0,0), guide='prism_offset')+
  labs(x='MAIT cell', y='Inverse Simpson Diversity Index')+
  theme_classic()
  theme(panel.grid.major.y=element_line(linetype='dotted', linewidth=0.35, color='grey'))
rm(ttest, TCR_diversity_plot_data, TCR_diversity_plot)

# Fig. 5H
# Occupied TCR repertoire space
TCR_proportion_plot           <- clonalProportion(list_TRAV12_sampled, cloneCall='CTstrict', group.by='Donor_CD56status', clonalSplit=c(5,10,20,50,550))
TCR_proportion_plot_data      <- TCR_proportion_plot[['data']] %>% separate(Var1, into=c('Donor', 'CD56status'), sep='_', remove=TRUE)
TCR_proportion_plot_data$Var2 <- fct_recode(TCR_proportion_plot_data$Var2, '1-5'='[1:5]', '6-10'='[6:10]', '11-20'='[11:20]', '21-50'='[21:50]', '>50'='[51:550]')

ggplot(TCR_proportion_plot_data, aes(x=CD56status, y=value, fill=Var2))+
  geom_col(position='fill', color='black')+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1), labels=scales::percent_format(scale=100), expand=c(0,0), guide='prism_offset')+
  scale_fill_brewer(palette='Greys', direction=-1)+
  labs(x='MAIT cell population', y='Occupied repertoire space', fill='Clones')+
  scale_x_discrete(labels=c('CD56-', 'CD56+'))+
  theme_classic()+
  theme(axis.text.x=element_text(angle=45, hjust=1), strip.background=element_blank())+
  facet_grid(.~Donor)

# Fig. S9D
# only focus on VDJ gene information without CDR3 sequence details
TRAV12$Cell_id <- rownames(TRAV12)
TRAV12_seurat <- subset(fullTCR, cells=TRAV12$Cell_id)

circles   <- getCirclize(TRAV12_seurat, cloneCall='CTgene', group.by='Donor_CD56status')
grid.cols <- viridis(length(unique(fullTCR$Donor_CD56status)))
names(grid.cols) <- levels(fullTCR@meta.data$Donor_CD56status)
circlize::chordDiagram(circles, self.link=1, grid.col=grid.cols, link.lwd=0.35)
dev.off()
rm(circles, grid.cols)

rm(list=ls())
gc()
