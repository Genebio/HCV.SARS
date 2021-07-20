setwd("C:/Users/user/Downloads/HCV.SARS")
library(KEGGprofile)
library(openxlsx)
library(dplyr)
library(reshape)
library(ggplot2)
db <- read.csv("C:/Users/user/Downloads/human_gene_id.csv")
head(db)

#Import DE genes
HCV.NEW_up <- read.xlsx("HCV.NEW_up.xlsx", colNames = F)
HCV.NEW_down <- read.xlsx("HCV.NEW_down.xlsx", colNames = F)
SARS_up <- read.xlsx("SARS_up.xlsx", colNames = F)
SARS_down <- read.xlsx("SARS_down.xlsx", colNames = F)
SARS.HCV_up <- read.xlsx("SARS.HCV_up.xlsx", colNames = F)
SARS.HCV_down <- read.xlsx("SARS.HCV_down.xlsx", colNames = F)

#Summarize UP and DOWN
HCV.NEW_DE <- c(HCV.NEW_up$X1, HCV.NEW_down$X1)
SARS_DE <- c(SARS_up$X1, SARS_down$X1)
SARS.HCV_DE <- c(SARS.HCV_up$X1, SARS.HCV_down$X1)

HCV.NEW_db <- as.character(db$GeneID[na.omit(match(HCV.NEW_DE, toupper(db$Symbol)))])
KEGG_HCV.NEW <- find_enriched_pathway(HCV.NEW_db, species = "hsa", returned_adjpvalue = 0.01, download_latest=T)[[1]][,c(1,6)]
KEGG_HCV.NEW_diff_pathways <- as.character(KEGG_HCV.NEW$Pathway_Name)
KEGG_HCV.NEW_all_pathways <- find_enriched_pathway(HCV.NEW_db, returned_pvalue=1, returned_adjpvalue = 1, species = "hsa",download_latest=T)[[1]][,c(1,6)]

SARS_db <- as.character(db$GeneID[na.omit(match(SARS_DE, toupper(db$Symbol)))])
KEGG_SARS <- find_enriched_pathway(SARS_db, species = "hsa",returned_adjpvalue = 0.01, download_latest=T)[[1]][,c(1,6)]
KEGG_SARS_diff_pathways <- as.character(KEGG_SARS$Pathway_Name)
KEGG_SARS_all_pathways <- find_enriched_pathway(SARS_db, returned_pvalue=1, returned_adjpvalue = 1, species = "hsa",download_latest=T)[[1]][,c(1,6)]

SARS.HCV_db <- as.character(db$GeneID[na.omit(match(SARS.HCV_DE, toupper(db$Symbol)))])
KEGG_SARS.HCV <- find_enriched_pathway(SARS.HCV_db, species = "hsa",returned_adjpvalue = 0.01, download_latest=T)[[1]][,c(1,6)]
KEGG_SARS.HCV <- KEGG_SARS.HCV %>% arrange(pvalueAdj, .by_group = T)
KEGG_SARS.HCV_diff_pathways <- as.character(KEGG_SARS.HCV$Pathway_Name) #top 30 SARS.HCV pathways
KEGG_SARS.HCV_all_pathways <- find_enriched_pathway(SARS.HCV_db, returned_pvalue=1, returned_adjpvalue = 1, species = "hsa", download_latest=T)[[1]][,c(1,6)]
##################################################################################

#PLOT unique/union DE pathways
Unique_DE_pathways <- unique(c(KEGG_HCV.NEW_diff_pathways, KEGG_SARS_diff_pathways, KEGG_SARS.HCV_diff_pathways))[-7]
HCV.NEW_df <- KEGG_HCV.NEW_all_pathways[match(Unique_DE_pathways, KEGG_HCV.NEW_all_pathways$Pathway_Name),]
SARS_df <- KEGG_SARS_all_pathways[match(Unique_DE_pathways, KEGG_SARS_all_pathways$Pathway_Name),]
SARS.HCV_df <- KEGG_SARS.HCV_all_pathways[match(Unique_DE_pathways, KEGG_SARS.HCV_all_pathways$Pathway_Name),]
merged_pathway_df <- data.frame(Pathway_name=Unique_DE_pathways, SARS.HCV=SARS.HCV_df$pvalueAdj, HCV.NEW=HCV.NEW_df$pvalueAdj, SARS=SARS_df$pvalueAdj)
merged_pathway_df[,2:ncol(merged_pathway_df)] <- -log10(merged_pathway_df[,2:ncol(merged_pathway_df)])
merged_pathway_df[is.na(merged_pathway_df)] <- 0
merged_pathway_df[merged_pathway_df==Inf] <- max(as.numeric(merged_pathway_df[mapply(is.finite, merged_pathway_df)]))#max non-infinite value
merged_pathway_df <- merged_pathway_df %>% arrange(-(HCV.NEW), .by_group = T)

names(merged_pathway_df) <- c("Pathway_name", "SARS+HCV","HCV", "SARS")
melted_DE_pathways <- melt(merged_pathway_df, id.vars = "Pathway_name")
melted_DE_pathways$value[melted_DE_pathways$value>10] <- 10
str(melted_DE_pathways)
pdf("Unique_DE_pathways_Padj_0.01_REF.pdf", height = 14, width = 14)
ggplot(melted_DE_pathways, aes(x = Pathway_name, y = value, fill=variable)) + 
  geom_bar(position=position_dodge(width=0.85), stat = "identity", width=0.7) + coord_flip() + theme_bw()+
  scale_y_continuous(name=bquote(~-Log[10]~ '(adjusted p-value)')) +
  scale_x_discrete(name="", limits=rev(merged_pathway_df$Pathway_name)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank(), legend.position = c("top"),
        axis.text.x = element_text(face="bold", size=20, angle=0, color = 'black'),
        axis.text.y = element_text(face="bold", size=17, angle=0, color = 'black'),
        axis.title.x = element_text(size=17, face="bold"),
        legend.text = element_text(face="bold", size=17, angle=0, color = 'black'),
        panel.border = element_rect(colour = "black", fill=NA, size=2))
dev.off()

########################################################

#Write genes of DE pathways in .xlsx file
#HCV enriched
HCV_pathways_total <- find_enriched_pathway(HCV.NEW_db, species = "hsa",returned_adjpvalue = 0.0007, refGene=union(REF_db,HCV.NEW_db), download_latest=T)
Pathway_names <- HCV_pathways_total$stastic$Pathway_Name
names(HCV_pathways_total$detail) <- Pathway_names
for (i in 1:length(HCV_pathways_total$detail)){
  HCV_pathways_total$detail[[i]] <- db$Symbol[match(as.numeric(HCV_pathways_total$detail[[i]]), db$GeneID)]
}
xx <- lapply(HCV_pathways_total$detail, unlist)
max <- max(sapply(xx, length))
HCV_pathways_total <- do.call(cbind, lapply(xx, function(z)c(z, rep(NA, max-length(z)))))
HCV_pathways_total <- data.frame(HCV_pathways_total)
names(HCV_pathways_total) <- Pathway_names
write.xlsx(HCV_pathways_total, file="HCV_pathways_total_Padj.01_REF_3.xlsx")

#SARS enriched
SARS_pathways_total <- find_enriched_pathway(SARS_db, species = "hsa",returned_adjpvalue = 0.0007, refGene=union(REF_db,SARS_db), download_latest=T)
Pathway_names <- SARS_pathways_total$stastic$Pathway_Name
names(SARS_pathways_total$detail) <- Pathway_names
for (i in 1:length(SARS_pathways_total$detail)){
  SARS_pathways_total$detail[[i]] <- db$Symbol[match(as.numeric(SARS_pathways_total$detail[[i]]), db$GeneID)]
}
xx <- lapply(SARS_pathways_total$detail, unlist)
max <- max(sapply(xx, length))
SARS_pathways_total <- do.call(cbind, lapply(xx, function(z)c(z, rep(NA, max-length(z)))))
SARS_pathways_total <- data.frame(SARS_pathways_total)
names(SARS_pathways_total) <- Pathway_names
write.xlsx(SARS_pathways_total, file="SARS_pathways_total_Padj.01_REF_3.xlsx")

#SARS.HCV enriched
SARS.HCV_pathways_total <- find_enriched_pathway(SARS.HCV_db, species = "hsa",returned_adjpvalue = 0.0007, refGene=union(REF_db,SARS.HCV_db), download_latest=T)
Pathway_names <- SARS.HCV_pathways_total$stastic$Pathway_Name
names(SARS.HCV_pathways_total$detail) <- Pathway_names
for (i in 1:length(SARS.HCV_pathways_total$detail)){
  SARS.HCV_pathways_total$detail[[i]] <- db$Symbol[match(as.numeric(SARS.HCV_pathways_total$detail[[i]]), db$GeneID)]
}
xx <- lapply(SARS.HCV_pathways_total$detail, unlist)
max <- max(sapply(xx, length))
SARS.HCV_pathways_total <- do.call(cbind, lapply(xx, function(z)c(z, rep(NA, max-length(z)))))
SARS.HCV_pathways_total <- data.frame(SARS.HCV_pathways_total)
names(SARS.HCV_pathways_total) <- Pathway_names
write.xlsx(SARS.HCV_pathways_total, file="SARS.HCV_pathways_total_Padj.01_REF_3.xlsx")


#Heatmap of DE pathways enriched in any of groups
Unique_DE_pathways <- unique(c(KEGG_HCV.NEW_diff_pathways, KEGG_SARS_diff_pathways, KEGG_SARS.HCV_diff_pathways))
HCV.NEW_df <- KEGG_HCV.NEW_all_pathways[match(Unique_DE_pathways, KEGG_HCV.NEW_all_pathways$Pathway_Name),]
SARS_df <- KEGG_SARS_all_pathways[match(Unique_DE_pathways, KEGG_SARS_all_pathways$Pathway_Name),]
SARS.HCV_df <- KEGG_SARS.HCV_all_pathways[match(Unique_DE_pathways, KEGG_SARS.HCV_all_pathways$Pathway_Name),]
merged_pathway_df <- data.frame(Pathway_name=Unique_DE_pathways, SARS.HCV=SARS.HCV_df$pvalueAdj, HCV.NEW=HCV.NEW_df$pvalueAdj, SARS=SARS_df$pvalueAdj)
merged_pathway_df[,2:ncol(merged_pathway_df)] <- -log10(merged_pathway_df[,2:ncol(merged_pathway_df)])
merged_pathway_df[,2:ncol(merged_pathway_df)] <- log(merged_pathway_df[,2:ncol(merged_pathway_df)])
merged_pathway_df[is.na(merged_pathway_df)] <- 0
merged_pathway_df[mapply(is.infinite, merged_pathway_df)] <- max(as.numeric(merged_pathway_df[mapply(is.finite, merged_pathway_df)]))

numeric_df <- merged_pathway_df[,2:ncol(merged_pathway_df)]
numeric_df[numeric_df>2] <- 2
numeric_df[numeric_df<c(-2)] <- c(-2)
croped_df <- data.frame(numeric_df, row.names = merged_pathway_df$Pathway_name)
colnames(croped_df) <- c("SARS+HCV","HCV", "SARS")
library(pheatmap)
pdf("Heatmap_log_croped_Padj.01.pdf")#, width = 10, height = 14
pheatmap(croped_df, angle_col=0, scale = "none", fontsize_row=10, fontsize_col=10)#fontsize_row=10, fontsize_col=15, 
dev.off()

#Plot Venn diagram of Up- and Down-regulated
