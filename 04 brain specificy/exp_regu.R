library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggsci)

#Boxplots showing expression level and specificity of promoters with or without a regulating enhancer in each FANTOM5 
sample.options(stringsAsFactors = F)

# e-p link for each individual
e_p <- read.table('/Users/songlt/Desktop/fantom/link.single_sample/stats.txt')
e_p$edge <- paste(e_p$V1, e_p$V2, sep='--')
e_p_sum_gene <- as.data.frame(xtabs(~V4+V1,e_p))
rownames(e_p_sum_gene) <- paste(e_p_sum_gene$V1,e_p_sum_gene$V4, sep='--')

meta<- read.table('/Users/songlt/Desktop/fantom/metadata.FANTOM/metadata.FANTOM.biosample.txt', header = T,sep='\t',quote = '')

########## The distribution of enhancer promoter link across different tissues
active_m <- c()
for(sta in unique(e_p$V4)){
  active_m <- rbind(active_m,c(length(unique(subset(e_p,V4==sta)[,'V1'])), length(unique(subset(e_p,V4==sta)[,'V2'])), length(unique(subset(e_p,V4==sta)[,'edge']))))
}
active_m <- as.data.frame(active_m)
rownames(active_m) <- unique(states$V4)
colnames(active_m) <- c('promoter','enhancer','edge')
rownames(meta) <- meta$BiosampleID
active_m$source <- meta[rownames(active_m),'c']
active_m_mlt <- melt(active_m)

ggplot(active_m_mlt,aes(x = reorder(source,value,median),y=value,color=source,fill=source))+
  geom_boxplot(alpha=0)+ylab('')+xlab('')+
  facet_wrap(.~variable,nrow = 3,scales="free_y")+
  theme_bw()+
  scale_color_manual(breaks = rownames(color_code),
                     values = color_code$ColorCode) +
  scale_fill_manual(breaks = rownames(color_code),
                    values = color_code$ColorCode) +
  theme(axis.text.x =element_text(angle = 45, hjust = 1,vjust = 1), legend.position = '')

dev.print(pdf, file='~/Desktop/fantom/figures/active_elemet_box2_individual.pdf')

# expression level
exp_level <- (read.table('/Users/songlt/Desktop/fantom/promoter_signal.matrix', row.names = 1, header = T,check.names = F))
exp_level_m <- melt(as.matrix(exp_level))
exp_level_m <- subset(exp_level_m, value>1)
rownames(exp_level_m) <-  paste(exp_level_m$Var1,exp_level_m$Var2, sep='--')


# Expression level of promoters with or without a regulating enhancer in each FANTOM5 sample.
e_p_sum_gene <- e_p_sum_gene[intersect(rownames(e_p_sum_gene),rownames(exp_level_m)),]

exp_level_m$Freq[!rownames(exp_level_m)%in%rownames(e_p_sum_gene)] <- 0
exp_level_m$Freq[is.na(exp_level_m$Freq)] <- e_p_sum_gene[ rownames(exp_level_m)[is.na(exp_level_m$Freq)], 'Freq']

exp_level_m$reg <- 'Nonregulated' 
exp_level_m$reg[exp_level_m$Freq!=0] <- 'Regulated' 

exp_level_m$Freq[exp_level_m$Freq>=5] <- '5+'

ggplot(exp_level_m,aes(x = reg,y=value,color=reg,fill=reg))+
  geom_boxplot(alpha=0,width=0.7)+ylab('Expression level')+xlab('')+
  theme_bw()+theme( legend.position = '')+ylim(0,30)+
  geom_signif(test='wilcox.test', map_signif_level=TRUE,y_position=28,
              comparisons = list(c("Nonregulated", "Regulated") ))+scale_fill_jco()+scale_color_jco()

dev.print(pdf, file='~/Desktop/fantom/figures/box_reg_exp.pdf')

ggplot(exp_level_m,aes(x = Freq,y=value,color=Freq,fill=Freq))+
  geom_boxplot(alpha=0, width=0.7)+ylab('Expression level')+xlab('Number of regulating enhancers')+
  theme_bw()+theme( legend.position = '')+
  geom_signif(test='wilcox.test', map_signif_level=TRUE,y_position=c(40,41,42,43),
              comparisons = list(c("1", "2"),c('2','3'),c('3','4' ),c('4','5+'))) +scale_fill_jco() +scale_color_jco()+ylim(0,45)

dev.print(pdf, file='~/Desktop/fantom/figures/box_numreg_exp.pdf')

# sample-specific expression specificity of promoters with or without a regulating enhancer in each FANTOM5 sample.
RowS <- rowSums(exp_level)
div <- exp_level/RowS
div2 <- div*log2(div)
div3 <- rowSums(div2)
exp_spec <- log2(ncol(exp_level))+div3+log2(div)

exp_spec_m <- melt(as.matrix(exp_spec))
rownames(exp_spec_m) <-  paste(exp_spec_m$Var1,exp_spec_m$Var2, sep='--')

exp_spec_m$Freq[!rownames(exp_spec_m)%in%rownames(e_p_sum_gene)] <- 0
exp_spec_m$Freq[is.na(exp_spec_m$Freq)] <- e_p_sum_gene[ rownames(exp_spec_m)[is.na(exp_spec_m$Freq)], 'Freq']

exp_spec_m$reg <- 'Nonregulated' 
exp_spec_m$reg[exp_spec_m$Freq!=0] <- 'Regulated' 

ggplot(exp_spec_m,aes(x = reg,y=value,color=reg,fill=reg))+
  geom_boxplot(alpha=0)+ylab('Expression specificity')+xlab('')+
  theme_bw()+theme( legend.position = '')+
  geom_signif(test='wilcox.test', map_signif_level=TRUE,y_position=-6,size = 0.4, textsize = 2,
              comparisons = list(c("Nonregulated", "Regulated") ))+ylim(-15,-5)+
  scale_fill_jco()+scale_color_jco()#

dev.print(pdf, file='~/Desktop/fantom/figures/box_reg_speci.pdf')

exp_spec_m$Freq[exp_spec_m$Freq>=5] <- '5+'

ggplot(exp_spec_m,aes(x = Freq,y=value,color=Freq))+
  geom_boxplot(alpha=0)+ylab('Expression specificity')+xlab('Number of regulating enhancers')+
  theme_bw()+theme( legend.position = '')  +
  geom_signif(test='wilcox.test', map_signif_level=TRUE,y_position=c(-6,-5,-4,-3,-2), textsize = 2,
              comparisons = list(c("1", "2"),c('2','3'),c('3','4' ),c('4','5+'))) +scale_fill_jco() +scale_color_jco()+ylim(-15,-2)
 dev.print(pdf, file='~/Desktop/fantom/figures/box_numreg_speci.pdf')


