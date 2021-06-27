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



# SARS.HCV_unique_db <- as.character(db$GeneID[na.omit(match(SARS.HCV_unique, toupper(db$Symbol)))])
# KEGG_SARS.HCV_unique <- find_enriched_pathway(SARS.HCV_unique_db, species = "hsa", download_latest=T)[[1]][,c(1,6)]
# KEGG_SARS.HCV_unique_diff_pathways <- as.character(KEGG_SARS.HCV_unique$Pathway_Name)
# KEGG_SARS.HCV_unique_all_pathways <- find_enriched_pathway(SARS.HCV_unique_db, returned_pvalue=1, returned_adjpvalue = 1, species = "hsa",download_latest=T)[[1]][,c(1,6)]


#Unique_DE_pathways <- unique(c(venn_var$x_only, DE_pathways$Pathway_name))
Unique_DE_pathways <- KEGG_HCV.NEW_diff_pathways
HCV.NEW_df <- KEGG_HCV.NEW_all_pathways[match(Unique_DE_pathways, KEGG_HCV.NEW_all_pathways$Pathway_Name),]
SARS_df <- KEGG_SARS_all_pathways[match(Unique_DE_pathways, KEGG_SARS_all_pathways$Pathway_Name),]
SARS.HCV_df <- KEGG_SARS.HCV_all_pathways[match(Unique_DE_pathways, KEGG_SARS.HCV_all_pathways$Pathway_Name),]
#SARS.HCV_unique_df <- KEGG_SARS.HCV_unique_all_pathways[match(Unique_DE_pathways, KEGG_SARS.HCV_unique_all_pathways$Pathway_Name),]


merged_pathway_df <- data.frame(Pathway_name=Unique_DE_pathways, SARS.HCV=SARS.HCV_df$pvalueAdj, HCV.NEW=HCV.NEW_df$pvalueAdj, SARS=SARS_df$pvalueAdj)
merged_pathway_df[,2:ncol(merged_pathway_df)] <- -log10(merged_pathway_df[,2:ncol(merged_pathway_df)])
merged_pathway_df <- merged_pathway_df %>% arrange(.by_group = T)
merged_pathway_df[is.na(merged_pathway_df)] <- 0

DE_pathways <- merged_pathway_df %>% filter((SARS.HCV-HCV)>1.2 & (SARS.HCV-SARS)>1.2) %>% arrange(rev(SARS.HCV), .by_group = T)
DE_pathways <- merged_pathway_df %>% arrange(SARS.HCV, .by_group = T)  
names(DE_pathways) <- c("Pathway_name", "SARS.HCV","HCV", "SARS")

DE_pathways <- DE_pathways %>% arrange(-(SARS.HCV), .by_group = T)
nrow(DE_pathways)
View(DE_pathways)
melted_DE_pathways <- melt(DE_pathways, id.vars = "Pathway_name")
pdf("SARS.HCV_DE_pathways_non_overlapping.pdf", width = 12)
ggplot(melted_DE_pathways, aes(x = Pathway_name, y = value, fill=variable)) + 
  geom_bar(position=position_dodge(width=0.85), stat = "identity", width=0.7) + coord_flip() + theme_bw()+
  scale_y_continuous(name=bquote(~-Log[10]~ '(adjusted p-value)')) +
  scale_x_discrete(name="", limits=rev(DE_pathways$Pathway_name)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank(), legend.position = c("top"),
        axis.text.x = element_text(face="bold", size=20, angle=0, color = 'black'),
        axis.text.y = element_text(face="bold", size=17, angle=0, color = 'black'),
        axis.title.x = element_text(size=14, face="bold"),
        legend.text = element_text(face="bold", size=17, angle=0, color = 'black'),
        panel.border = element_rect(colour = "black", fill=NA, size=2))
dev.off()
?scale_y_continuous
?geom_bar
vignette("ggplot2-specs")

part1 <- merged_pathway_df[1:18,]
part2 <- merged_pathway_df[19:37,]
names(merged_pathway_df) <- c("Pathway_name", "SARS.HCV","HCV", "SARS")
merged_pathway_df <- merged_pathway_df %>% arrange(SARS.HCV, SARS.HCV-HCV & SARS.HCV-SARS, .by_group = T)
pathway_order <- merged_pathway_df$Pathway_name

library(reshape)
melted_pathway_df <- melt(merged_pathway_df, id.vars = "Pathway_name")
melted_pathway_part1 <- melt(part1, id.vars = "Pathway_name")
melted_pathway_part2 <- melt(part2, id.vars = "Pathway_name")

#pdf("pathway_SARS.HCV_part2_cut_1e-5_filtered.pdf")
pdf("melted_pathway.pdf", width = 14, height = 14)
ggplot(melted_pathway_df, aes(x = Pathway_name, y = value, fill=variable)) + 
  geom_bar(position="dodge", stat = "identity") + coord_flip() + theme_bw()+
  scale_y_continuous(name="") +
  scale_x_discrete(name="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank(), legend.position = c("top"),
        axis.text.x = element_text(face="bold", size=10, angle=0, color = 'black'),
        axis.text.y = element_text(face="bold", size=10, angle=0, color = 'black'))
dev.off()

SARS.HCV_pathways_total <- find_enriched_pathway(SARS.HCV_db, species = "hsa",returned_adjpvalue = 0.05, download_latest=T)
View(SARS.HCV_pathways_total$stastic)
p1 <- SARS.HCV_pathways_total$detail$`04668`
p2 <- SARS.HCV_pathways_total$detail$`05134`
p3 <- SARS.HCV_pathways_total$detail$`05100`
p4 <- SARS.HCV_pathways_total$detail$`00220`
p5 <- SARS.HCV_pathways_total$detail$`04657`
p6 <- SARS.HCV_pathways_total$detail$`04666`
p7 <- SARS.HCV_pathways_total$detail$`04922`
p8 <- SARS.HCV_pathways_total$detail$`01210`
p9 <- SARS.HCV_pathways_total$detail$`04979`
p10 <- SARS.HCV_pathways_total$detail$`00450`
p11 <- SARS.HCV_pathways_total$detail$`05030`
p12 <- SARS.HCV_pathways_total$detail$`03320`
p13 <- SARS.HCV_pathways_total$detail$`00020`


# p1 <- SARS.HCV_pathways_total$detail$`04668`
# p2 <- SARS.HCV_pathways_total$detail$`05134`
# p3 <- SARS.HCV_pathways_total$detail$`05100`
# p4 <- SARS.HCV_pathways_total$detail$`04657`
# p5 <- SARS.HCV_pathways_total$detail$`04666`
# p6 <- SARS.HCV_pathways_total$detail$`03320`
# # p7 <- SARS.HCV_pathways_total$detail$`04666`
# # p8 <- SARS.HCV_pathways_total$detail$`00513`
# # p9 <- SARS.HCV_pathways_total$detail$`03320`
# # p10 <- SARS.HCV_pathways_total$detail$`04137`
# # p11 <- SARS.HCV_pathways_total$detail$`00511`


#Match the genes
p1_genes <- db$Symbol[match(as.numeric(p1), db$GeneID)]
p2_genes <- db$Symbol[match(as.numeric(p2), db$GeneID)]
p3_genes <- db$Symbol[match(as.numeric(p3), db$GeneID)]
p4_genes <- db$Symbol[match(as.numeric(p4), db$GeneID)]
p5_genes <- db$Symbol[match(as.numeric(p5), db$GeneID)]
p6_genes <- db$Symbol[match(as.numeric(p6), db$GeneID)]
p7_genes <- db$Symbol[match(as.numeric(p7), db$GeneID)]
p8_genes <- db$Symbol[match(as.numeric(p8), db$GeneID)]
p9_genes <- db$Symbol[match(as.numeric(p9), db$GeneID)]
p10_genes <- db$Symbol[match(as.numeric(p10), db$GeneID)]
p11_genes <- db$Symbol[match(as.numeric(p11), db$GeneID)]
p12_genes <- db$Symbol[match(as.numeric(p12), db$GeneID)]
p13_genes <- db$Symbol[match(as.numeric(p13), db$GeneID)]

Pathway_names <- DE_pathways$Pathway_name


#Save genes from pathways
write.xlsx(p1_genes, file="p1_genes.xlsx")
write.xlsx(p2_genes, file="p2_genes.xlsx")
write.xlsx(p3_genes, file="p3_genes.xlsx")
write.xlsx(p4_genes, file="p4_genes.xlsx")
write.xlsx(p5_genes, file="p5_genes.xlsx")
write.xlsx(p6_genes, file="p6_genes.xlsx")
write.xlsx(p7_genes, file="p7_genes.xlsx")
write.xlsx(p8_genes, file="p8_genes.xlsx")
write.xlsx(p9_genes, file="p9_genes.xlsx")
write.xlsx(p10_genes, file="p10_genes.xlsx")
write.xlsx(p11_genes, file="p11_genes.xlsx")
write.xlsx(p12_genes, file="p12_genes.xlsx")
write.xlsx(p13_genes, file="p13_genes.xlsx")
Pathway_names <- DE_pathways$Pathway_name
write.xlsx(Pathway_names, file="Pathway_names.xlsx")
View(test_df)


load("NEW_matrix_downsized.Rdata")
#save(NEW_matrix_downsized, file="NEW_matrix_downsized.Rdata")
str(NEW_matrix_downsized)

#Plot SARS.HCV_unique pathways (log ratio vs. control)
p1_matrix <- NEW_matrix_downsized[rownames(NEW_matrix_downsized)%in%p1_genes,]
SARS.HCV <- log(rowMeans(p1_matrix[,1:3]/rowMeans(p1_matrix[,10:12])))
SARS <- log(rowMeans(p1_matrix[,4:6])/rowMeans(p1_matrix[,10:12]))
HCV <- log(rowMeans(p1_matrix[,7:9])/rowMeans(p1_matrix[,10:12]))
p1_df <- cbind(HCV,SARS,SARS.HCV)
p1_df[is.infinite(p1_df)] <- 0
p1_df <- data.frame(p1_df)
p1_df$genes <- rownames(p1_df)
library(reshape)
melted_p1_df <- melt(p1_df, id.vars = "genes")
pdf("pathway_p1.pdf")
boxplot(melted_p1_df$value~melted_p1_df$variable)
dev.off()
#View(p1_df)

p2_matrix <- NEW_matrix_downsized[rownames(NEW_matrix_downsized)%in%p2_genes,]
SARS.HCV <- log(rowMeans(p2_matrix[,1:3]/rowMeans(p2_matrix[,10:12])))
SARS <- log(rowMeans(p2_matrix[,4:6])/rowMeans(p2_matrix[,10:12]))
HCV <- log(rowMeans(p2_matrix[,7:9])/rowMeans(p2_matrix[,10:12]))
p2_df <- cbind(HCV,SARS,SARS.HCV)
p2_df[is.infinite(p2_df)] <- 0
p2_df <- data.frame(p2_df)
p2_df$genes <- rownames(p2_df)
library(reshape)
melted_p2_df <- melt(p2_df, id.vars = "genes")
pdf("pathway_p2.pdf")
boxplot(melted_p2_df$value~melted_p2_df$variable)
dev.off()
#View(p2_df)

p3_matrix <- NEW_matrix_downsized[rownames(NEW_matrix_downsized)%in%p3_genes,]
SARS.HCV <- log(rowMeans(p3_matrix[,1:3]/rowMeans(p3_matrix[,10:12])))
SARS <- log(rowMeans(p3_matrix[,4:6])/rowMeans(p3_matrix[,10:12]))
HCV <- log(rowMeans(p3_matrix[,7:9])/rowMeans(p3_matrix[,10:12]))
p3_df <- cbind(HCV,SARS,SARS.HCV)
p3_df[is.infinite(p3_df)] <- 0
p3_df <- data.frame(p3_df)
p3_df$genes <- rownames(p3_df)
library(reshape)
melted_p3_df <- melt(p3_df, id.vars = "genes")
pdf("pathway_p3.pdf")
boxplot(melted_p3_df$value~melted_p3_df$variable)
dev.off()
#View(p3_df)

p4_matrix <- NEW_matrix_downsized[rownames(NEW_matrix_downsized)%in%p4_genes,]
SARS.HCV <- log(rowMeans(p4_matrix[,1:3]/rowMeans(p4_matrix[,10:12])))
SARS <- log(rowMeans(p4_matrix[,4:6])/rowMeans(p4_matrix[,10:12]))
HCV <- log(rowMeans(p4_matrix[,7:9])/rowMeans(p4_matrix[,10:12]))
p4_df <- cbind(HCV,SARS,SARS.HCV)
p4_df[is.infinite(p4_df)] <- 0
p4_df <- data.frame(p4_df)
p4_df$genes <- rownames(p4_df)
library(reshape)
melted_p4_df <- melt(p4_df, id.vars = "genes")
pdf("pathway_p4.pdf")
boxplot(melted_p4_df$value~melted_p4_df$variable)
dev.off()
#View(p4_df)


p5_matrix <- NEW_matrix_downsized[rownames(NEW_matrix_downsized)%in%p5_genes,]
SARS.HCV <- log(rowMeans(p5_matrix[,1:3]/rowMeans(p5_matrix[,10:12])))
SARS <- log(rowMeans(p5_matrix[,4:6])/rowMeans(p5_matrix[,10:12]))
HCV <- log(rowMeans(p5_matrix[,7:9])/rowMeans(p5_matrix[,10:12]))
p5_df <- cbind(HCV,SARS,SARS.HCV)
p5_df[is.infinite(p5_df)] <- 0
p5_df <- data.frame(p5_df)
p5_df$genes <- rownames(p5_df)
library(reshape)
melted_p5_df <- melt(p5_df, id.vars = "genes")
pdf("pathway_p5.pdf")
boxplot(melted_p5_df$value~melted_p5_df$variable)
dev.off()
#View(p5_df)

p6_matrix <- NEW_matrix_downsized[rownames(NEW_matrix_downsized)%in%p6_genes,]
SARS.HCV <- log(rowMeans(p6_matrix[,1:3]/rowMeans(p6_matrix[,10:12])))
SARS <- log(rowMeans(p6_matrix[,4:6])/rowMeans(p6_matrix[,10:12]))
HCV <- log(rowMeans(p6_matrix[,7:9])/rowMeans(p6_matrix[,10:12]))
p6_df <- cbind(HCV,SARS,SARS.HCV)
p6_df[is.infinite(p6_df)] <- 0
p6_df <- data.frame(p6_df)
p6_df$genes <- rownames(p6_df)
library(reshape)
melted_p6_df <- melt(p6_df, id.vars = "genes")
pdf("pathway_p6.pdf")
boxplot(melted_p6_df$value~melted_p6_df$variable)
dev.off()
#View(p6_df)

p7_matrix <- NEW_matrix_downsized[rownames(NEW_matrix_downsized)%in%p7_genes,]
SARS.HCV <- log(rowMeans(p7_matrix[,1:3]/rowMeans(p7_matrix[,10:12])))
SARS <- log(rowMeans(p7_matrix[,4:6])/rowMeans(p7_matrix[,10:12]))
HCV <- log(rowMeans(p7_matrix[,7:9])/rowMeans(p7_matrix[,10:12]))
p7_df <- cbind(HCV,SARS,SARS.HCV)
p7_df[is.infinite(p7_df)] <- 0
p7_df <- data.frame(p7_df)
p7_df$genes <- rownames(p7_df)
library(reshape)
melted_p7_df <- melt(p7_df, id.vars = "genes")
pdf("pathway_p7.pdf")
boxplot(melted_p7_df$value~melted_p7_df$variable)
dev.off()
#View(p7_df)

p8_matrix <- NEW_matrix_downsized[rownames(NEW_matrix_downsized)%in%p8_genes,]
SARS.HCV <- log(rowMeans(p8_matrix[,1:3]/rowMeans(p8_matrix[,10:12])))
SARS <- log(rowMeans(p8_matrix[,4:6])/rowMeans(p8_matrix[,10:12]))
HCV <- log(rowMeans(p8_matrix[,7:9])/rowMeans(p8_matrix[,10:12]))
p8_df <- cbind(HCV,SARS,SARS.HCV)
p8_df[is.infinite(p8_df)] <- 0
p8_df <- data.frame(p8_df)
p8_df$genes <- rownames(p8_df)
library(reshape)
melted_p8_df <- melt(p8_df, id.vars = "genes")
pdf("pathway_p8.pdf")
boxplot(melted_p8_df$value~melted_p8_df$variable)
dev.off()
#View(p8_df)

p9_matrix <- NEW_matrix_downsized[rownames(NEW_matrix_downsized)%in%p9_genes,]
SARS.HCV <- log(rowMeans(p9_matrix[,1:3]/rowMeans(p9_matrix[,10:12])))
SARS <- log(rowMeans(p9_matrix[,4:6])/rowMeans(p9_matrix[,10:12]))
HCV <- log(rowMeans(p9_matrix[,7:9])/rowMeans(p9_matrix[,10:12]))
p9_df <- cbind(HCV,SARS,SARS.HCV)
p9_df[is.infinite(p9_df)] <- 0
p9_df <- data.frame(p9_df)
p9_df$genes <- rownames(p9_df)
library(reshape)
melted_p9_df <- melt(p9_df, id.vars = "genes")
pdf("pathway_p9.pdf")
boxplot(melted_p9_df$value~melted_p9_df$variable)
dev.off()
#View(p9_df)

p10_matrix <- NEW_matrix_downsized[rownames(NEW_matrix_downsized)%in%p10_genes,]
SARS.HCV <- log(rowMeans(p10_matrix[,1:3]/rowMeans(p10_matrix[,10:12])))
SARS <- log(rowMeans(p10_matrix[,4:6])/rowMeans(p10_matrix[,10:12]))
HCV <- log(rowMeans(p10_matrix[,7:9])/rowMeans(p10_matrix[,10:12]))
p10_df <- cbind(HCV,SARS,SARS.HCV)
p10_df[is.infinite(p10_df)] <- 0
p10_df <- data.frame(p10_df)
p10_df$genes <- rownames(p10_df)
library(reshape)
melted_p10_df <- melt(p10_df, id.vars = "genes")
pdf("pathway_p10.pdf")
boxplot(melted_p10_df$value~melted_p10_df$variable)
dev.off()
#View(p10_df)

p11_matrix <- NEW_matrix_downsized[rownames(NEW_matrix_downsized)%in%p11_genes,]
SARS.HCV <- log(rowMeans(p11_matrix[,1:3]/rowMeans(p11_matrix[,10:12])))
SARS <- log(rowMeans(p11_matrix[,4:6])/rowMeans(p11_matrix[,10:12]))
HCV <- log(rowMeans(p11_matrix[,7:9])/rowMeans(p11_matrix[,10:12]))
p11_df <- cbind(HCV,SARS,SARS.HCV)
p11_df[is.infinite(p11_df)] <- 0
p11_df <- data.frame(p11_df)
p11_df$genes <- rownames(p11_df)
library(reshape)
melted_p11_df <- melt(p11_df, id.vars = "genes")
pdf("pathway_p11.pdf")
boxplot(melted_p11_df$value~melted_p11_df$variable)
dev.off()
#View(p11_df)

p12_matrix <- NEW_matrix_downsized[rownames(NEW_matrix_downsized)%in%p12_genes,]
SARS.HCV <- log(rowMeans(p12_matrix[,1:3]/rowMeans(p12_matrix[,10:12])))
SARS <- log(rowMeans(p12_matrix[,4:6])/rowMeans(p12_matrix[,10:12]))
HCV <- log(rowMeans(p12_matrix[,7:9])/rowMeans(p12_matrix[,10:12]))
p12_df <- cbind(HCV,SARS,SARS.HCV)
p12_df[is.infinite(p12_df)] <- 0
p12_df <- data.frame(p12_df)
p12_df$genes <- rownames(p12_df)
library(reshape)
melted_p12_df <- melt(p12_df, id.vars = "genes")

p13_matrix <- NEW_matrix_downsized[rownames(NEW_matrix_downsized)%in%p13_genes,]
SARS.HCV <- log(rowMeans(p13_matrix[,1:3]/rowMeans(p13_matrix[,10:12])))
SARS <- log(rowMeans(p13_matrix[,4:6])/rowMeans(p13_matrix[,10:12]))
HCV <- log(rowMeans(p13_matrix[,7:9])/rowMeans(p13_matrix[,10:12]))
p13_df <- cbind(HCV,SARS,SARS.HCV)
p13_df[is.infinite(p13_df)] <- 0
p13_df <- data.frame(p13_df)
p13_df$genes <- rownames(p13_df)
library(reshape)
melted_p13_df <- melt(p13_df, id.vars = "genes")



#Plot 6 DE pathways toagether
pdf("SARS.HCV_pathways_1_5.pdf")
par(mfrow=c(3,2), cex.axis=1.5, cex.main=1.5, font.axis=2, las=1, lwd=1.5, mar=c(5, 4, 4, 2))
boxplot(melted_p1_df$value~melted_p1_df$variable, xlab="", ylab="", font=2, main=Pathway_names[1])
boxplot(melted_p2_df$value~melted_p2_df$variable, xlab="", ylab="", font=2, main=Pathway_names[2])
boxplot(melted_p3_df$value~melted_p3_df$variable, xlab="", ylab="", font=2, main=Pathway_names[3])
boxplot(melted_p4_df$value~melted_p4_df$variable, xlab="", ylab="", font=2, main=Pathway_names[4])
boxplot(melted_p5_df$value~melted_p5_df$variable, xlab="", ylab="", font=2, main=Pathway_names[5])
dev.off()

pdf("SARS.HCV_pathways_6_10.pdf")
par(mfrow=c(3,2), cex.axis=1.5, cex.main=1.5, font.axis=2, las=1, lwd=1.5, mar=c(5, 4, 4, 2))
boxplot(melted_p6_df$value~melted_p6_df$variable, xlab="", ylab="", font=2, main=Pathway_names[6])
boxplot(melted_p7_df$value~melted_p7_df$variable, xlab="", ylab="", font=2, main=Pathway_names[7])
boxplot(melted_p8_df$value~melted_p8_df$variable, xlab="", ylab="", font=2, main=Pathway_names[7])
boxplot(melted_p9_df$value~melted_p9_df$variable, xlab="", ylab="", font=2, main=Pathway_names[9])
boxplot(melted_p10_df$value~melted_p10_df$variable, xlab="", ylab="", font=2, main=Pathway_names[10])
dev.off()

pdf("SARS.HCV_pathways_11_13.pdf", width = 3.5)
par(mfrow=c(3,1), cex.axis=1.5, cex.main=1.5, font.axis=2, las=1, lwd=1.5, mar=c(5, 4, 4, 2))
boxplot(melted_p11_df$value~melted_p11_df$variable, xlab="", ylab="", font=2, main=Pathway_names[6])
boxplot(melted_p12_df$value~melted_p12_df$variable, xlab="", ylab="", font=2, main=Pathway_names[7])
boxplot(melted_p13_df$value~melted_p13_df$variable, xlab="", ylab="", font=2, main=Pathway_names[7])
dev.off()







# pdf("SARS.HCV_p7_p11.pdf")
# par(mfrow=c(3,2), cex.axis=1.5, cex.main=1.5, font.axis=2, las=1, lwd=1.5, mar=c(5, 4, 4, 2))
# boxplot(melted_p7_df$value~melted_p7_df$variable, xlab="", ylab="", font=2, main=Pathway_names[7])
# boxplot(melted_p8_df$value~melted_p8_df$variable, xlab="", ylab="", font=2, main=Pathway_names[8])
# boxplot(melted_p9_df$value~melted_p9_df$variable, xlab="", ylab="", font=2, main=Pathway_names[9])
# boxplot(melted_p10_df$value~melted_p10_df$variable, xlab="", ylab="", font=2, main=Pathway_names[10])
# boxplot(melted_p11_df$value~melted_p11_df$variable, xlab="", ylab="", font=2, main=Pathway_names[11])
# dev.off()
# 

?par

# #Plot SARS.HCV_unique pathways
# p1_matrix <- NEW_matrix_downsized[rownames(NEW_matrix_downsized)%in%p1_genes,]
# SARS.HCV <- log(rowMeans(p1_matrix[,1:3]))
# SARS <- log(rowMeans(p1_matrix[,4:6]))
# HCV <- log(rowMeans(p1_matrix[,7:9]))
# Control <- log(rowMeans(p1_matrix[,10:12]))
# p1_df <- cbind(Control,HCV,SARS,SARS.HCV)
# p1_df[is.infinite(p1_df)] <- 0
# p1_df <- data.frame(p1_df)
# p1_df$genes <- rownames(p1_df)
# library(reshape)
# melted_p1_df <- melt(p1_df, id.vars = "genes")
# pdf("pathway_p1.pdf")
# boxplot(melted_p1_df$value~melted_p1_df$variable)
# dev.off()
# View(p1_df)
# 
# p2_matrix <- NEW_matrix_downsized[rownames(NEW_matrix_downsized)%in%p2_genes,]
# SARS.HCV <- log(rowMeans(p2_matrix[,1:3]))
# SARS <- log(rowMeans(p2_matrix[,4:6]))
# HCV <- log(rowMeans(p2_matrix[,7:9]))
# Control <- log(rowMeans(p2_matrix[,10:12]))
# p2_df <- cbind(Control,HCV,SARS,SARS.HCV)
# p2_df[is.infinite(p2_df)] <- 0
# p2_df <- data.frame(p2_df)
# p2_df$genes <- rownames(p2_df)
# library(reshape)
# melted_p2_df <- melt(p2_df, id.vars = "genes")
# pdf("pathway_p2.pdf")
# boxplot(melted_p2_df$value~melted_p2_df$variable)
# dev.off()
# View(p2_df)
# 
# p3_matrix <- NEW_matrix_downsized[rownames(NEW_matrix_downsized)%in%p3_genes,]
# SARS.HCV <- log(rowMeans(p3_matrix[,1:3]))
# SARS <- log(rowMeans(p3_matrix[,4:6]))
# HCV <- log(rowMeans(p3_matrix[,7:9]))
# Control <- log(rowMeans(p3_matrix[,10:12]))
# p3_df <- cbind(Control,HCV,SARS,SARS.HCV)
# p3_df[is.infinite(p3_df)] <- 0
# p3_df <- data.frame(p3_df)
# p3_df$genes <- rownames(p3_df)
# library(reshape)
# melted_p3_df <- melt(p3_df, id.vars = "genes")
# pdf("pathway_p3.pdf")
# boxplot(melted_p3_df$value~melted_p3_df$variable)
# dev.off()
# View(p3_df)
# 
# p4_matrix <- NEW_matrix_downsized[rownames(NEW_matrix_downsized)%in%p4_genes,]
# SARS.HCV <- log(rowMeans(p4_matrix[,1:3]))
# SARS <- log(rowMeans(p4_matrix[,4:6]))
# HCV <- log(rowMeans(p4_matrix[,7:9]))
# Control <- log(rowMeans(p4_matrix[,10:12]))
# p4_df <- cbind(Control,HCV,SARS,SARS.HCV)
# p4_df[is.infinite(p4_df)] <- 0
# p4_df <- data.frame(p4_df)
# p4_df$genes <- rownames(p4_df)
# library(reshape)
# melted_p4_df <- melt(p4_df, id.vars = "genes")
# pdf("pathway_p4.pdf")
# boxplot(melted_p4_df$value~melted_p4_df$variable)
# dev.off()
# View(p4_df)





