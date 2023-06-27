# 2021.06.28
# pca promoter
# SONG Liting

library(ggplot2)
library(tsne)
library(Rtsne)
library(reshape2)

# meta data
metadata_dataset <- read.table('/home1/yangyc/share/data_resource/FANTOM/metadata/metadata.FANTOM.dataset.txt',sep='\t',header = T)
metadata_biosample <- read.table('/home1/yangyc/share/data_resource/FANTOM/metadata/metadata.FANTOM.biosample.txt',sep='\t',header = T,quote = '')
colorcode <- read.table('/home1/GENE_proc/SONGLITING/FANTOM/figures/colorcode.txt',sep='\t',header = T,comment.char = '',stringsAsFactors = F)
rownames(metadata_biosample) <- metadata_biosample$BiosampleID

# promoter coverage matrix
# sample: sample name; v14: gene; depth: depth
promoter <- read.table('/home1/GENE_proc/SONGLITING/FANTOM/promoter/all_promoter.txt',stringsAsFactors = F)[,c(1,4,13,17)]
colnames(promoter) <- c('chr','gene','depth','sample')
promoter$sample_gene <- paste(promoter$gene,promoter$sample,sep='_')
promoter <- promoter[!duplicated(promoter$sample_gene),]
promoter <- subset(promoter, !chr%in%c("chrX",  "chrY" , "chrM"))  
promoter_m <- acast(promoter,sample~gene,value.var=c('depth'))

#promoter_m <- colSums(promoter_m>1)

# all mapped reads
depth <- read.table('/home1/GENE_proc/SONGLITING/FANTOM/depth/mapped_seq.txt')
depth$n_million <- depth$V1/1000000
rownames(depth) <- metadata_dataset$BiosampleID
depth <- depth[rownames(promoter_m),]

# promoter_tpm
promoter_tpm <- promoter_m/depth$n_million
hist(promoter_tpm[1,],xlim = c(1,100),breaks = 100000)

metadata_biosample <- metadata_biosample[rownames(promoter_tpm),]

# PCA all gene promoter
pca_all <- prcomp(t(promoter_tpm), scale. = T, rank. = 2) 
pca_all_meta <- cbind(pca_all$rotation,metadata_biosample)


ggplot(pca_all_meta, aes(x=PC1,y=PC2,color=BiosampleGroup,shape=BiosampleGroup))+
  geom_point()+
  scale_color_manual(breaks = colorcode$BiosampleGroup,
                     values = colorcode$ColorCode) + theme_bw()+
  scale_shape_manual(values = c(0:25,0:3))
  
ggplot(pca_all_meta, aes(x=PC1,y=PC2,color=BiosampleGroup))+
  geom_point()+
  scale_color_manual(breaks = colorcode$BiosampleGroup,
                     values = colorcode$ColorCode) + theme_bw()

dev.print(pdf, file='/home1/GENE_proc/SONGLITING/FANTOM/figures/pca_promoter_all.pdf')

# tsne all
set.seed(1234)
tsne_out <- Rtsne(promoter_tpm,pca=FALSE,dims=2,
                  perplexity=30,theta=0.0) # Run TSNE

tsne_res <- as.data.frame(tsne_out$Y)
colnames(tsne_res) <- c("tSNE1","tSNE2")
tsne_res$BiosampleGroup <- metadata_biosample$BiosampleGroup
head(tsne_res)

ggplot(tsne_res,aes(tSNE1,tSNE2,color=BiosampleGroup))+
  geom_point()+
  scale_color_manual(breaks = colorcode$BiosampleGroup,
                     values = colorcode$ColorCode) + theme_bw()
dev.print(pdf, file='/home1/GENE_proc/SONGLITING/FANTOM/figures/tsne_promoter_all.pdf')

ggplot(tsne_res,aes(tSNE1,tSNE2,color=BiosampleGroup,shape=BiosampleGroup))+
  geom_point()+
  scale_color_manual(breaks = colorcode$BiosampleGroup,
                     values = colorcode$ColorCode) + theme_bw()+
  scale_shape_manual(values = c(0:25,0:3))

dev.print(pdf, file='/home1/GENE_proc/SONGLITING/FANTOM/figures/tsne_promoter_all_shape.pdf')


# tsne top
set.seed(1234)
promoter_tpm_top <- promoter_tpm[, order(apply(promoter_tpm,2,mad), decreasing = T)[1:2000]]

tsne_out_top <- Rtsne(promoter_tpm_top,pca=FALSE,dims=2,
                  perplexity=30,theta=0.0) # Run TSNE

tsne_res_top <- as.data.frame(tsne_out_top$Y)
colnames(tsne_res_top) <- c("tSNE1","tSNE2")
tsne_res_top$BiosampleGroup <- metadata_biosample$BiosampleGroup

ggplot(tsne_res_top,aes(tSNE1,tSNE2,color=BiosampleGroup))+
  geom_point()+
  scale_color_manual(breaks = colorcode$BiosampleGroup,
                     values = colorcode$ColorCode) + theme_bw()

dev.print(pdf, file='/home1/GENE_proc/SONGLITING/FANTOM/figures/tsne_promoter_top.pdf')

# PCA top 2000 promoter
pca_top <- prcomp(t(promoter_tpm_top), scale. = T, rank. = 2) 
pca_top_meta <- cbind(pca_top$rotation,metadata_biosample)

ggplot(pca_top_meta, aes(x=PC1,y=PC2,color=metadata_biosample$BiosampleGroup))+
  geom_point()+
  scale_color_manual(breaks = colorcode$BiosampleGroup,
                     values = colorcode$ColorCode) + theme_bw()

dev.print(pdf, file='/home1/GENE_proc/SONGLITING/FANTOM/figures/pca_promoter_top.pdf')

save.image(file='/home1/GENE_proc/SONGLITING/FANTOM/promoter/promoter.RData')
write.table(t(promoter_tpm), file='/home1/GENE_proc/SONGLITING/FANTOM/promoter/combine_promoter_tpm.txt',quote = F,sep='\t')
write.table(t(promoter_m), file='/home1/GENE_proc/SONGLITING/FANTOM/promoter/combine_promoter_count.txt',quote = F,sep='\t')

