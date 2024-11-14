rm(list=ls())

install.packages('pak')
pak::pkg_install("r-lib/rlang")
devtools::install_github("igordot/msigdbr")
BiocManager::install("fgsea")

install.packages('ggplot2')

library(data.table)
library(fgsea)
library(ggplot2)
library(BiocParallel)
library(msigdbr)
library(xlsx)
library(dplyr)


setwd("D:/KP/mouse_brain_project/fGSEA")
######################################################################################################
####################### fgsea using Wilcoxon gene rank list from Scanpy ##############################
######################################################################################################

###################### KO to others ################################################################

### KEGG ###

DEG <- read.csv("GC1_DEG_by_tissue_for_fgesa.csv", sep = ',')
DEG_1 <- DEG %>% select(KO_names, KO_score, KO_pvals)
DEG_1 %>% filter(KO_pvals < 0.05)
DEG_1$gene <- DEG_1$KO_names
row.names(DEG_1) <- DEG_1$gene

##### rank with log2FC may not correct ####
#### I used new rank (avg_log2FC * -log10(p.adj)) ####calculate and make new column in Excel ####
ranks_1 <- DEG_1$KO_score
#ranks_1 <- DEG_1$rank
names(ranks_1) <- row.names(DEG_1)
head(ranks_1)
fgsea.p1<-barplot(sort(ranks_1, decreasing = TRUE))

msigdbr_show_species()
genesets_1 = msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")
msigdbr_list = split(x = genesets_1$gene_symbol, f = genesets_1$gs_name)
fgseaRes_1 <- fgsea(msigdbr_list, ranks_1, minSize=15, maxSize = 500)
head(fgseaRes_1[order(padj, -abs(NES)), ], n=10)

topUp <- fgseaRes_1 %>% filter(ES > 0) %>% top_n(5, wt=-padj)
topDown <- fgseaRes_1 %>% filter(ES < 0) %>% top_n(5, wt=-padj)
topPathways <- bind_rows(topUp, topDown) %>% arrange(-ES)
fgsea.p2<-plotGseaTable(msigdbr_list[topPathways$pathway], ranks_1, fgseaRes_1, gseaParam = 0.5)

fgseaResTidy <- fgseaRes_1 %>%  as_tibble() %>% arrange(desc(NES))

fwrite(fgseaResTidy, file="GC1_KO_GSEA_results_KEGG.csv", sep=",", sep2=c("", " ", ""))


ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

#### REACTOME  ####

genesets_2 = msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME")
msigdbr_list_2 = split(x = genesets_2$gene_symbol, f = genesets_2$gs_name)
msigdbr_list_2
fgseaRes_2 <- fgsea(pathways=msigdbr_list_2, ranks_1, minSize=5, maxSize = 500)
head(fgseaRes_2[order(padj, -abs(NES)), ], n=10)


topUp_2 <- fgseaRes_2 %>% filter(ES > 0) %>% top_n(5, wt=-padj)
topDown_2 <- fgseaRes_2 %>% filter(ES < 0) %>% top_n(5, wt=-padj)
topPathways_2 <- bind_rows(topUp_2, topDown_2) %>% arrange(-ES)


fgsea.p3<-plotGseaTable(msigdbr_list_2[topPathways_2$pathway], ranks_1, fgseaRes_2, gseaParam = 0.5)


fgseaResTidy_2 <- fgseaRes_2 %>%  as_tibble() %>% arrange(desc(NES))
fgseaResTidy_2 %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

fwrite(fgseaResTidy_2, file="GC1_KO_GSEA_results_REACTOME.csv", sep=",", sep2=c("", " ", ""))
fgseaResTidy_2 = read.csv('GC1_KO_GSEA_results_REACTOME.csv', sep=",")


ggplot(fgseaResTidy_2, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="REACTOME Hallmark pathways") + 
  theme_minimal()


####### GOBP
genesets_5 = msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")
msigdbr_list_5 = split(x = genesets_5$gene_symbol, f = genesets_5$gs_name)
msigdbr_list_5
fgseaRes_5 <- fgsea(pathways=msigdbr_list_5, ranks_1, minSize=5, maxSize = 500)
head(fgseaRes_5[order(padj, -abs(NES)), ], n=10)


topUp_5 <- fgseaRes_5 %>% filter(ES > 0) %>% top_n(5, wt=-padj)
topDown_5 <- fgseaRes_5 %>% filter(ES < 0) %>% top_n(5, wt=-padj)
topPathways_5 <- bind_rows(topUp_5, topDown_5) %>% arrange(-ES)

dev.off() 
fgsea.p5<-plotGseaTable(msigdbr_list_5[topPathways_5$pathway], ranks_DEG, fgseaRes_5, gseaParam = 0.5)



fgseaResTidy_5 <- fgseaRes_5 %>%  as_tibble() %>% arrange(desc(NES))
fgseaResTidy_5 %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

fwrite(fgseaResTidy_5, file="GC1_KO_GSEA_results_GOBP.csv", sep=",", sep2=c("", " ", ""))
fgseaResTidy_5 = read.csv("GC1_KO_GSEA_results_GOBP.csv", sep = ",")


ggplot(fgseaResTidy_5, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="KOvs.others GOBP pathways") + 
  theme_minimal()

