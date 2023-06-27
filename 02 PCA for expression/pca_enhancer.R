# 2021.06.30
# pca enhancer
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

# enhancer coverage matrix
# sample: sample name; v14: encode_acc; depth: depth
enhancer <- read.table('/home1/GENE_proc/SONGLITING/FANTOM/enhancer/filter_enhancer.txt',stringsAsFactors = F,sep='\t')
colnames(enhancer) <- c('chr','encode_acc','depth','sample')
#enhancer$sample_encode_acc <- paste(enhancer$encode_acc,enhancer$sample,sep='_')
#enhancer <- enhancer[!duplicated(enhancer$sample_encode_acc),]
enhancer <- subset(enhancer, !chr%in%c("chrX",  "chrY" , "chrM"))  
enhancer_m <- acast(enhancer,sample~encode_acc,value.var=c('depth'))
enhancer_m[is.na(enhancer_m)] <- 0
#enhancer_m1 <- enhancer_m[,colSums(enhancer_m>5) > 100]

# all mapped reads
depth <- read.table('/home1/GENE_proc/SONGLITING/FANTOM/depth/mapped_seq.txt')
depth$n_million <- depth$V1/1000000
rownames(depth) <- metadata_dataset$BiosampleID
depth <- depth[rownames(enhancer_m),]

# enhancer_tpm
enhancer_tpm <- enhancer_m/depth$n_million
hist(enhancer_tpm[1,],xlim = c(1,100),breaks = 100000)

metadata_biosample <- metadata_biosample[rownames(enhancer_tpm),]

# Cumulative Distribution for samples

set.seed(1234)
ecdf <- melt(t(enhancer_m[sample(1:nrow(enhancer_m),10),]))
ggplot(ecdf,aes(x = value,color=Var2))+
  stat_ecdf()+
  labs(title = "",x='depth',y='Cumulative Distribution')+xlim(c(1,20))
dev.print(pdf, file='/home1/GENE_proc/SONGLITING/FANTOM/figures/sample_cdf.pdf')

# Cumulative Distribution for enhancers

set.seed(1234)
ecdf_enhancer <- melt(t(enhancer_m[,sample(1:nrow(enhancer_m),20)]))
ggplot(ecdf_enhancer,aes(x = value,color=Var1))+
  stat_ecdf()+
  labs(title = "",x='depth',y='Cumulative Distribution')+xlim(c(1,20))

dev.print(pdf, file='/home1/GENE_proc/SONGLITING/FANTOM/figures/enhancer_cdf.pdf')


# PCA all encode_acc enhancer
pca_all <- prcomp(t(enhancer_tpm), scale. = T, rank. = 2) 
pca_all_meta <- cbind(pca_all$rotation,metadata_biosample)

ggplot(pca_all_meta, aes(x=PC1,y=PC2,color=BiosampleGroup))+
  geom_point()+
  scale_color_manual(breaks = colorcode$BiosampleGroup,
                     values = colorcode$ColorCode) + theme_bw()

dev.print(pdf, file='/home1/GENE_proc/SONGLITING/FANTOM/figures/pca_enhancer_all.pdf')

# tsne all
set.seed(1234)
tsne_out <- Rtsne(enhancer_tpm,pca=FALSE,dims=2,
                  perplexity=30,theta=0.0) # Run TSNE

tsne_res <- as.data.frame(tsne_out$Y)
colnames(tsne_res) <- c("tSNE1","tSNE2")
tsne_res$BiosampleGroup <- metadata_biosample$BiosampleGroup
head(tsne_res)

ggplot(tsne_res,aes(tSNE1,tSNE2,color=BiosampleGroup))+
  geom_point()+
  scale_color_manual(breaks = colorcode$BiosampleGroup,
                     values = colorcode$ColorCode) + theme_bw()
dev.print(pdf, file='/home1/GENE_proc/SONGLITING/FANTOM/figures/tsne_enhancer_all.pdf')
#dev.print(pdf, file='/home1/GENE_proc/SONGLITING/FANTOM/figures/tsne_enhancer_exp.pdf')

ggplot(tsne_res,aes(tSNE1,tSNE2,color=BiosampleGroup,shape=BiosampleGroup))+
  geom_point()+
  scale_color_manual(breaks = colorcode$BiosampleGroup,
                     values = colorcode$ColorCode) + theme_bw()+
  scale_shape_manual(values = c(0:25,0:3))

dev.print(pdf, file='/home1/GENE_proc/SONGLITING/FANTOM/figures/tsne_enhancer_all_shape.pdf')


# tsne top
set.seed(1234)
enhancer_tpm_top <- enhancer_tpm[, order(apply(enhancer_tpm,2,mad), decreasing = T)[1:2000]]

tsne_out_top <- Rtsne(enhancer_tpm_top,pca=FALSE,dims=2,
                      perplexity=30,theta=0.0) # Run TSNE

tsne_res_top <- as.data.frame(tsne_out_top$Y)
colnames(tsne_res_top) <- c("tSNE1","tSNE2")
tsne_res_top$BiosampleGroup <- metadata_biosample$BiosampleGroup

ggplot(tsne_res_top,aes(tSNE1,tSNE2,color=BiosampleGroup))+
  geom_point()+
  scale_color_manual(breaks = colorcode$BiosampleGroup,
                     values = colorcode$ColorCode) + theme_bw()

dev.print(pdf, file='/home1/GENE_proc/SONGLITING/FANTOM/figures/tsne_enhancer_top.pdf')

# PCA top 2000 enhancer
pca_top <- prcomp(t(enhancer_tpm_top), scale. = T, rank. = 2) 
pca_top_meta <- cbind(pca_top$rotation,metadata_biosample)

ggplot(pca_top_meta, aes(x=PC1,y=PC2,color=metadata_biosample$BiosampleGroup))+
  geom_point()+
  scale_color_manual(breaks = colorcode$BiosampleGroup,
                     values = colorcode$ColorCode) + theme_bw()

dev.print(pdf, file='/home1/GENE_proc/SONGLITING/FANTOM/figures/pca_enhancer_top.pdf')

save.image(file='/home1/GENE_proc/SONGLITING/FANTOM/enhancer/enhancer.RData')
write.table(t(enhancer_tpm), file='/home1/GENE_proc/SONGLITING/FANTOM/enhancer/combine_enhancer_tpm.txt',quote = F,sep='\t')
write.table(t(enhancer_m), file='/home1/GENE_proc/SONGLITING/FANTOM/enhancer/combine_enhancer_count.txt',quote = F,sep='\t')

