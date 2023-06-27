# Prioritizing genes associated with brain disorders by leveraging enhancer-promoter interactions in diverse neural cells and tissues

These scripts are part of codes for Yucheng Yang's project entitled "Prioritizing genes associated with brain disorders by leveraging enhancer-promoter interactions in diverse neural cells and tissues". 

## 01 make expression matrix

1. get_enhancer_cov2.sh & 1.get_promoter_cov2.sh: We get the expression (counts) using bedtools for promoter and enhancer regions, where the promoter regions were defined as the regions raging from up and down-stream 500bp at TSS (0.get_promoter_region.sh), and the enhancer regions were extracted form fantom dataset.


## 02 specificity of the expression of enhancers and promoters

We visulized the tissue specificity of expression of enhancers and promoters using principle compotent analysis (PCA) and TSNE.


## 03 specificity of the reconstructed EPIs (figure 2)

1. cosine_sim.R: Clustering of FANTOM5 samples based on cosine similarities of their enhancer-promoter networks
2. network_specific_X.R: Visualization of the enhancer-promoter networks active in only a single group of FANTOM5 tissues.Here we only presented the subnetworks on Chr1.
3. exp_regu.R: Boxplots showing expression level and specificity of promoters with or without a regulating enhancer in each 
FANTOM5 sample, and with various numbers of regulating enhancers in each FANTOM5 sample.

## 04 Specificity and functions of CREs in human brain (figure 3)
1. brain_enhancer_category_count.R: To identify brain-elevated and brain period-elevated enhancers for each stage
2. brain_promoter_category_count.R: To identify brain-elevated and brain period-elevated promoters for each stage
3. brain_path.R: biological process enrichment for the genes with distinct tissue- and stage-specific activity patterns
4. homer.sh and TF_enrich.R:The enrichment of TF binding in the enhancers with distinct tissue- and stage-specific activity patterns
5. enhancer_conserv.R: Conservation of enhancers with distinct tissu and stage-specific activity patterns


## 05 Expression profiles of the genes associated with diverse disorders and behavioral-cognitive phenotypes from different brain single cell types at the fetal and adult stage (figure 5)
cell_exp_2211.R: The gene expression data of brain single cells were collected from our STAB database. We used single cell datasets from 19-26 PCW and 40-60 Years in the analysis at the fetal and adult stage, respectively.

## Help
You can visit 'https://github.com/Soulnature/brainepl/' for more codes in this project. If you have any questions, please feel free to contact me (ltsong18@.fudan.edu.cn).
