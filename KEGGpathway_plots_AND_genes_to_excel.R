library(KEGGprofile)
library(openxlsx)
db <- read.csv("human_gene_id.csv")
head(db)

HCV.NEW_db <- as.character(db$GeneID[na.omit(match(HCV.NEW_DE, toupper(db$Symbol)))])
KEGG_HCV.NEW <- find_enriched_pathway(HCV.NEW_db, species = "hsa", returned_adjpvalue = 0.05, download_latest=T)[[1]][,c(1,6)]
KEGG_HCV.NEW_diff_pathways <- as.character(KEGG_HCV.NEW$Pathway_Name)
KEGG_HCV.NEW_all_pathways <- find_enriched_pathway(HCV.NEW_db, returned_pvalue=1, returned_adjpvalue = 1, species = "hsa",download_latest=T)[[1]][,c(1,6)]

SARS_db <- as.character(db$GeneID[na.omit(match(SARS_DE, toupper(db$Symbol)))])
KEGG_SARS <- find_enriched_pathway(SARS_db, species = "hsa",returned_adjpvalue = 0.05, download_latest=T)[[1]][,c(1,6)]
KEGG_SARS_diff_pathways <- as.character(KEGG_SARS$Pathway_Name)
KEGG_SARS_all_pathways <- find_enriched_pathway(SARS_db, returned_pvalue=1, returned_adjpvalue = 1, species = "hsa",download_latest=T)[[1]][,c(1,6)]

SARS.HCV_db <- as.character(db$GeneID[na.omit(match(SARS.HCV_DE, toupper(db$Symbol)))])
KEGG_SARS.HCV <- find_enriched_pathway(SARS.HCV_db, species = "hsa",returned_adjpvalue = 0.05, download_latest=T)[[1]][,c(1,6)]
KEGG_SARS.HCV <- KEGG_SARS.HCV %>% arrange(pvalueAdj, .by_group = T)
KEGG_SARS.HCV_diff_pathways <- as.character(KEGG_SARS.HCV$Pathway_Name)[-2] #top 30 SARS.HCV pathways
KEGG_SARS.HCV_all_pathways <- find_enriched_pathway(SARS.HCV_db, returned_pvalue=1, returned_adjpvalue = 1, species = "hsa",download_latest=T)[[1]][,c(1,6)]

#HCV pathways DE
Unique_DE_pathways <- KEGG_HCV.NEW_diff_pathways
HCV.NEW_df <- KEGG_HCV.NEW_all_pathways[match(Unique_DE_pathways, KEGG_HCV.NEW_all_pathways$Pathway_Name),]
SARS_df <- KEGG_SARS_all_pathways[match(Unique_DE_pathways, KEGG_SARS_all_pathways$Pathway_Name),]
SARS.HCV_df <- KEGG_SARS.HCV_all_pathways[match(Unique_DE_pathways, KEGG_SARS.HCV_all_pathways$Pathway_Name),]
merged_pathway_df <- data.frame(Pathway_name=Unique_DE_pathways, SARS.HCV=SARS.HCV_df$pvalueAdj, HCV.NEW=HCV.NEW_df$pvalueAdj, SARS=SARS_df$pvalueAdj)
merged_pathway_df[,2:ncol(merged_pathway_df)] <- -log10(merged_pathway_df[,2:ncol(merged_pathway_df)])
merged_pathway_df <- merged_pathway_df %>% arrange(-(HCV.NEW), .by_group = T)
merged_pathway_df[is.na(merged_pathway_df)] <- 0
names(merged_pathway_df) <- c("Pathway_name", "SARS+HCV","HCV", "SARS")
melted_DE_pathways <- melt(merged_pathway_df, id.vars = "Pathway_name")
pdf("HCV_DE_pathways_paper.pdf", height = 14, width = 14)
ggplot(melted_DE_pathways, aes(x = Pathway_name, y = value, fill=variable)) + 
  geom_bar(position=position_dodge(width=0.85), stat = "identity", width=0.7) + coord_flip() + theme_bw()+
  scale_y_continuous(name=bquote(~-Log[10]~ '(adjusted p-value)')) +
  scale_x_discrete(name="", limits=rev(merged_pathway_df$Pathway_name)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank(), legend.position = c("top"),
        axis.text.x = element_text(face="bold", size=20, angle=0, color = 'black'),
        axis.text.y = element_text(face="bold", size=17, angle=0, color = 'black'),
        axis.title.x = element_text(size=14, face="bold"),
        legend.text = element_text(face="bold", size=17, angle=0, color = 'black'),
        panel.border = element_rect(colour = "black", fill=NA, size=2))
dev.off()

HCV_pathways_total <- find_enriched_pathway(HCV.NEW_db, species = "hsa",returned_adjpvalue = 0.05, download_latest=T)
names(HCV_pathways_total$detail) <- HCV_pathways_total$stastic$Pathway_Name
for (i in 1:length(HCV_pathways_total$detail)){
  HCV_pathways_total$detail[[i]] <- db$Symbol[match(as.numeric(HCV_pathways_total$detail[[i]]), db$GeneID)]
}
# HCV_pathways_total <- HCV_pathways_total$detail
# save(HCV_pathways_total, file="HCV_pathways_total.Rdata")
xx <- lapply(HCV_pathways_total$detail, unlist)
max <- max(sapply(xx, length))
HCV_pathways_total <- do.call(cbind, lapply(xx, function(z)c(z, rep(NA, max-length(z)))))
HCV_pathways_total <- data.frame(HCV_pathways_total)
write.xlsx(HCV_pathways_total, file="HCV_pathways_total.xlsx")

#SARS pathways DE
Unique_DE_pathways <- KEGG_SARS_diff_pathways
HCV.NEW_df <- KEGG_HCV.NEW_all_pathways[match(Unique_DE_pathways, KEGG_HCV.NEW_all_pathways$Pathway_Name),]
SARS_df <- KEGG_SARS_all_pathways[match(Unique_DE_pathways, KEGG_SARS_all_pathways$Pathway_Name),]
SARS.HCV_df <- KEGG_SARS.HCV_all_pathways[match(Unique_DE_pathways, KEGG_SARS.HCV_all_pathways$Pathway_Name),]
merged_pathway_df <- data.frame(Pathway_name=Unique_DE_pathways, SARS.HCV=SARS.HCV_df$pvalueAdj, HCV.NEW=HCV.NEW_df$pvalueAdj, SARS=SARS_df$pvalueAdj)
merged_pathway_df[,2:ncol(merged_pathway_df)] <- -log10(merged_pathway_df[,2:ncol(merged_pathway_df)])
merged_pathway_df <- merged_pathway_df %>% arrange(-(SARS), .by_group = T)
merged_pathway_df[is.na(merged_pathway_df)] <- 0
names(merged_pathway_df) <- c("Pathway_name", "SARS+HCV","HCV", "SARS")
melted_DE_pathways <- melt(merged_pathway_df, id.vars = "Pathway_name")
pdf("SARS_DE_pathways_paper.pdf", height = 14, width = 14)
ggplot(melted_DE_pathways, aes(x = Pathway_name, y = value, fill=variable)) + 
  geom_bar(position=position_dodge(width=0.85), stat = "identity", width=0.7) + coord_flip() + theme_bw()+
  scale_y_continuous(name=bquote(~-Log[10]~ '(adjusted p-value)')) +
  scale_x_discrete(name="", limits=rev(merged_pathway_df$Pathway_name)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank(), legend.position = c("top"),
        axis.text.x = element_text(face="bold", size=20, angle=0, color = 'black'),
        axis.text.y = element_text(face="bold", size=17, angle=0, color = 'black'),
        axis.title.x = element_text(size=14, face="bold"),
        legend.text = element_text(face="bold", size=17, angle=0, color = 'black'),
        panel.border = element_rect(colour = "black", fill=NA, size=2))
dev.off()

SARS_pathways_total <- find_enriched_pathway(SARS_db, species = "hsa",returned_adjpvalue = 0.05, download_latest=T)
names(SARS_pathways_total$detail) <- SARS_pathways_total$stastic$Pathway_Name
for (i in 1:length(SARS_pathways_total$detail)){
  SARS_pathways_total$detail[[i]] <- db$Symbol[match(as.numeric(SARS_pathways_total$detail[[i]]), db$GeneID)]
}
# SARS_pathways_total <- SARS_pathways_total$detail
# save(SARS_pathways_total, file="SARS_pathways_total.Rdata")
xx <- lapply(SARS_pathways_total$detail, unlist)
max <- max(sapply(xx, length))
SARS_pathways_total <- do.call(cbind, lapply(xx, function(z)c(z, rep(NA, max-length(z)))))
SARS_pathways_total <- data.frame(SARS_pathways_total)
write.xlsx(SARS_pathways_total, file="SARS_pathways_total.xlsx")

#SARS.HCV pathways DE
Unique_DE_pathways <- KEGG_SARS.HCV_diff_pathways
HCV.NEW_df <- KEGG_HCV.NEW_all_pathways[match(Unique_DE_pathways, KEGG_HCV.NEW_all_pathways$Pathway_Name),]
SARS_df <- KEGG_SARS_all_pathways[match(Unique_DE_pathways, KEGG_SARS_all_pathways$Pathway_Name),]
SARS.HCV_df <- KEGG_SARS.HCV_all_pathways[match(Unique_DE_pathways, KEGG_SARS.HCV_all_pathways$Pathway_Name),]
merged_pathway_df <- data.frame(Pathway_name=Unique_DE_pathways, SARS.HCV=SARS.HCV_df$pvalueAdj, HCV.NEW=HCV.NEW_df$pvalueAdj, SARS=SARS_df$pvalueAdj)
merged_pathway_df[,2:ncol(merged_pathway_df)] <- -log10(merged_pathway_df[,2:ncol(merged_pathway_df)])
merged_pathway_df <- merged_pathway_df %>% arrange(-(SARS.HCV), .by_group = T)
merged_pathway_df[is.na(merged_pathway_df)] <- 0
names(merged_pathway_df) <- c("Pathway_name", "SARS+HCV","HCV", "SARS")
melted_DE_pathways <- melt(merged_pathway_df, id.vars = "Pathway_name")
pdf("SARS.HCV_DE_pathways_paper.pdf", height = 14, width = 14)
ggplot(melted_DE_pathways, aes(x = Pathway_name, y = value, fill=variable)) + 
  geom_bar(position=position_dodge(width=0.85), stat = "identity", width=0.7) + coord_flip() + theme_bw()+
  scale_y_continuous(name=bquote(~-Log[10]~ '(adjusted p-value)')) +
  scale_x_discrete(name="", limits=rev(merged_pathway_df$Pathway_name)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank(), legend.position = c("top"),
        axis.text.x = element_text(face="bold", size=20, angle=0, color = 'black'),
        axis.text.y = element_text(face="bold", size=17, angle=0, color = 'black'),
        axis.title.x = element_text(size=14, face="bold"),
        legend.text = element_text(face="bold", size=17, angle=0, color = 'black'),
        panel.border = element_rect(colour = "black", fill=NA, size=2))
dev.off()

SARS.HCV_pathways_total <- find_enriched_pathway(SARS.HCV_db, species = "hsa",returned_adjpvalue = 0.05, download_latest=T)
names(SARS.HCV_pathways_total$detail) <- SARS.HCV_pathways_total$stastic$Pathway_Name
for (i in 1:length(SARS.HCV_pathways_total$detail)){
  SARS.HCV_pathways_total$detail[[i]] <- db$Symbol[match(as.numeric(SARS.HCV_pathways_total$detail[[i]]), db$GeneID)]
}

xx <- lapply(SARS.HCV_pathways_total$detail, unlist)
max <- max(sapply(xx, length))
SARS.HCV_pathways_total <- do.call(cbind, lapply(xx, function(z)c(z, rep(NA, max-length(z)))))
SARS.HCV_pathways_total <- data.frame(SARS.HCV_pathways_total)
write.xlsx(SARS.HCV_pathways_total, file="SARS.HCV_pathways_total.xlsx")

#All unique pathways DE
Unique_DE_pathways <- unique(c(KEGG_HCV.NEW_diff_pathways, KEGG_SARS_diff_pathways, KEGG_SARS.HCV_diff_pathways))
HCV.NEW_df <- KEGG_HCV.NEW_all_pathways[match(Unique_DE_pathways, KEGG_HCV.NEW_all_pathways$Pathway_Name),]
SARS_df <- KEGG_SARS_all_pathways[match(Unique_DE_pathways, KEGG_SARS_all_pathways$Pathway_Name),]
SARS.HCV_df <- KEGG_SARS.HCV_all_pathways[match(Unique_DE_pathways, KEGG_SARS.HCV_all_pathways$Pathway_Name),]
merged_pathway_df <- data.frame(Pathway_name=Unique_DE_pathways, SARS.HCV=SARS.HCV_df$pvalueAdj, HCV.NEW=HCV.NEW_df$pvalueAdj, SARS=SARS_df$pvalueAdj)
merged_pathway_df[,2:ncol(merged_pathway_df)] <- -log10(merged_pathway_df[,2:ncol(merged_pathway_df)])
merged_pathway_df[,2:ncol(merged_pathway_df)] <- log(merged_pathway_df[,2:ncol(merged_pathway_df)])
merged_pathway_df[is.na(merged_pathway_df)] <- c(-2)
merged_pathway_df[mapply(is.infinite, merged_pathway_df)] <- c(-2)

numeric_df <- merged_pathway_df[,2:ncol(merged_pathway_df)]
numeric_df[numeric_df>2] <- 2
numeric_df[numeric_df<c(-2)] <- c(-2)
croped_df <- data.frame(numeric_df, row.names = merged_pathway_df$Pathway_name)
colnames(croped_df) <- c("SARS+HCV","HCV", "SARS")
pdf("Heatmap_log_croped.pdf", width = 10, height = 14)
pheatmap(croped_df, fontsize_row=10, fontsize_col=15, angle_col=0)
dev.off()
# names(merged_pathway_df) <- c("Pathway_name", "SARS+HCV","HCV", "SARS")
# merged_total_df <- data.frame(merged_pathway_df[,2:4], row.names = merged_pathway_df$Pathway_name)
library(pheatmap)
pdf("ALL_DE_total_heatmap_log_no_scaling.pdf", width = 10, height = 14)
pheatmap(merged_total_df, fontsize_row=10, fontsize_col=15, angle_col=0)#scale="row", 
dev.off()

log_merged_total_df <- log(merged_total_df)
pheatmap(log_merged_total_df, scale="row", fontsize_row=10, fontsize_col=15, angle_col=0)
# HCV_pathways_total
# SARS_pathways_total
# SARS.HCV_pathways_total
# unique_total_pathways <- unique(c(names(HCV_pathways_total), names(SARS_pathways_total), names(SARS.HCV_pathways_total)))
# total_list <- unlist(list(HCV_pathways_total), list(SARS_pathways_total), list(SARS.HCV_pathways_total))
# str(total_list)
# 
# test_HCV <- HCV_pathways_total$detail
# test_SARS <- SARS_pathways_total$detail
# test_SARS.HCV <- SARS.HCV_pathways_total$detail
# # lsNames <- c("test_HCV","test_SARS","test_SARS.HCV")
# # do.call(mapply, c(FUN=c, sapply(lsNames, as.symbol), SIMPLIFY=T))
# names(test_HCV)#47
# names(test_SARS)#34
# names(test_SARS.HCV)#38
# new_list <- list()
# new_list[1:47] <- test_HCV 
# names(new_list[1:47]) <- names(test_HCV)
# new_list[48:82] <- test_SARS 
# names(new_list[48:82]) <- names(test_SARS)
# 
# length(test_SARS.HCV)
# new_list[83:121] <- test_SARS.HCV 
# names(new_list[83:121]) <- names(test_SARS.HCV)
