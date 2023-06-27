
library(Seurat)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(cowplot)
library(base2grob)
library(gridExtra)
library(stringr)

setwd('~/Desktop/fantom/')
magma <- read.table('~/Desktop/fantom/state-signature genes_id_symbol.overlap.zscore.50%.fdr(1).txt',sep='\t',header = T)
colnames(magma)[5] <- 'gene'
magma <-subset(magma, gene!='' & state_name%in%c('Astrocyte - cerebellum','Astrocyte - cerebral cortex','Neurons',"Oligodendrocyte - precursors"))

# ensg <- read.table('~/Desktop/fantom/ens.txt',sep='\t',header = T)
# as.data.frame(cbind(str_split_fixed(ensg$X10000.101A1,'\\|',2)[,1],str_split_fixed(ensg$X10000.101A1,'\\|',2)[,2])) -> gene_info
# rownames(gene_info) <- gene_info[,1]
# 
# magma$gene <- ''
# for(i in 1:nrow(magma)){
#   magma$gene[i] <- paste(gene_info[unlist(strsplit( magma$signature.genes[i],' ')),2],collapse = ' ')
# }


stat <- magma

unique(stat$state_name)

load('~/Desktop/fantom/sub_pfc_p6_14.RData')
h.combined <- sub_pfc_p6_14
meta_data <- sub_pfc_p6_14@meta.data

setwd('~/Desktop/fantom')
cell_type_info <- read.table('./cell_types_v2.txt',row.names = 1, sep='\t',header = T)
period_infos <- read.table('./Period.txt',header = T,sep='\t',row.names = 1)
Area_infos <- read.table('./Area.txt',header = T,sep='\t',row.names = 1)
o.ids <- c(paste('InN', 1:10, sep=''), paste('ExN', 1:15,sep=''))
n.ids <- c(paste('InN', c(3, 6, 5, '1a', '4a', '4_5', '1b', '4b', '5_6', '7_8' ), sep=''),
           paste('ExN', c('1_4', '1a', 3, '6a', '2', '9', 10, '6b', 11, '8', '4_6', '1b', 5, '1c', 4 ), sep=''))
meta_d <- as.data.frame(table(meta_data[,c('Area','Period','cluster1')]))
meta_d <- meta_d[meta_d$Freq!=0,]

#Neurons   
options(stringsAsFactors = F)
i <- 1
stat_neuron <-  subset(stat, state_name=='Neurons')
h.combined_sub <- subset(h.combined, cluster%in%c('ExN','InN' ) & Period%in%c('P6','P14'))
dd_Neurons <- c()
for(disease in stat_neuron$gwas_file){
  
  gene <- unlist(strsplit(stat_neuron$gene[i],' '))
  d <- as.matrix(h.combined_sub[['RNA']]@data[ intersect(gene, rownames(h.combined_sub[['RNA']]@data)), ])
  dd <- melt(d, id = row.names)
  dd <- dd %>% dplyr::rename(gene = Var1, cell = Var2)
  dd$Area <- h.combined_sub$Area[dd$cell]
  dd$Period <- h.combined_sub$Period[dd$cell]
  
  #str(dd$tree.ident)
  #dd$gene <- factor(dd$gene, levels = intersect(gene, rownames(t[['RNA']]@data)))
  dd$disease <- disease
  dd_Neurons <- rbind(dd_Neurons, dd)
  i <- i+1
  
}


#Astrocyte - cerebellum  
#Astrocyte - cerebral cortex  

options(stringsAsFactors = F)
i <- 1
stat_Astro <-  subset(stat, state_name=='Astrocyte - cerebral cortex')
h.combined_sub <- subset(h.combined, cluster%in%c('Astro' ) & Period%in%c('P6','P14') )
dd_Astro <- c()
for(disease in stat_neuron$gwas_file){
  
  gene <- unlist(strsplit(stat_Astro$gene[i],' '))
  d <- as.matrix(h.combined_sub[['RNA']]@data[ intersect(gene, rownames(h.combined_sub[['RNA']]@data)), ])
  dd <- melt(d, id = row.names)
  dd <- dd %>% dplyr::rename(gene = Var1, cell = Var2)
  dd$Area <- h.combined_sub$Area[dd$cell]
  dd$Period <- h.combined_sub$Period[dd$cell]
  
  #str(dd$tree.ident)
  #dd$gene <- factor(dd$gene, levels = intersect(gene, rownames(t[['RNA']]@data)))
  dd$disease <- disease
  dd_Astro <- rbind(dd_Astro, dd)
  i <- i+1
  
}



#Oligodendrocyte - precursors
options(stringsAsFactors = F)
i <- 1
stat_neuron <-  subset(stat, state_name=='Oligodendrocyte - precursors')
h.combined_sub <- subset(h.combined, cluster%in%c('OPC' ) & Period%in%c('P6'))
dd_OPC <- c()
for(disease in stat_neuron$gwas_file){
  
  gene <- unlist(strsplit(stat_neuron$gene[i],' '))
  d <- as.matrix(h.combined_sub[['RNA']]@data[ intersect(gene, rownames(h.combined_sub[['RNA']]@data)), ])
  dd <- melt(d, id = row.names)
  dd <- dd %>% dplyr::rename(gene = Var1, cell = Var2)
  dd$Area <- h.combined_sub$Area[dd$cell]
  dd$Period <- h.combined_sub$Period[dd$cell]
  
  #str(dd$tree.ident)
  #dd$gene <- factor(dd$gene, levels = intersect(gene, rownames(t[['RNA']]@data)))
  dd$disease <- disease
  dd_OPC <- rbind(dd_OPC, dd)
  i <- i+1
  
}


dd_Astro_cellt <- subset(dd_Astro,Area=='PFC' )
dd_Astro_cellt$celltype <- 'Astro'
dd_Neurons_cellt <- subset(dd_Neurons,Area=='PFC' )
dd_Neurons_cellt$celltype <- 'Neuron'
dd_OPC_cellt <- subset(dd_OPC,Area=='PFC' )
dd_OPC_cellt$celltype <- 'OPC'

dd_cellt <- rbind(dd_Astro_cellt,dd_Neurons_cellt,dd_OPC_cellt)
dd_cellt <- subset(dd_cellt,!disease%in%c('AD_Kunkle','Cannabis','MDD','ALCDEP'))

dd_cellt <- within(dd_cellt,{
  disease[disease=='INSOMNIA'] <- 'Insomnia'
  disease[disease=='INTELLIGENCE'] <- 'Intelligence'
  disease[disease=='Neuroticism'] <- 'Neuroticism'
  disease[disease=='RISK_behavior'] <- 'Risk behavior'
  disease[disease=='EPILEPSY'] <- 'Epilepsy'
  disease[disease=='AD_Jansen'] <- 'AD'
  
  disease_cat <- 'Trait'
  disease_cat[disease%in%c('AD','Epilepsy','MS','PD','ALS')] <- 'Neurological'
  disease_cat[disease%in%c('MDD','SCZ','BIP','ASD','ADHD','AUD')] <- 'Psychiatric'
  disease_cat[is.na(disease_cat)] <- 'Trait'
}) 



dd_cellt_mean <- aggregate(value ~ gene+disease+Area+Period+ disease_cat+celltype, data = dd_cellt, mean)

ggboxplot(subset(dd_cellt_mean, Period==c('P6','P14')[2]), x = "disease", y = "value",
          color = "celltype", palette = "jco",size=0.5,add = "none",
          outlier.shape = NA)+
  facet_grid(Period~disease_cat,scales = 'free_x')+xlab('')+ylab('Normalized Expression')+
  stat_compare_means(aes(group = celltype), label = "p.signif")+
  theme(axis.text.x = element_text(angle = 44,vjust = 1,hjust = 1))+
  ylim(c(0,1))

dev.print(pdf,'~/Desktop/dd_cellt_mean_P14_update.pdf')

ggboxplot(subset(dd_cellt_mean, Period==c('P6','P14')[1]), x = "disease", y = "value",
          color = "celltype", palette = "jco",size=0.5,add = "none",
          outlier.shape = NA)+
  facet_grid(Period~disease_cat,scales = 'free_x')+xlab('')+ylab('Normalized Expression')+
  stat_compare_means(aes(group = celltype), label = "p.signif")+
  theme(axis.text.x = element_text(angle = 44,vjust = 1,hjust = 1))+
  ylim(c(0,1))

dev.print(pdf,'~/Desktop/dd_cellt_mean_p6_update.pdf')


