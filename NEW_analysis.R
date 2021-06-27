HCV.OLD_matrix_downsized
NEW_matrix_downsized
names(HCV.OLD_matrix_downsized)


#HCV.OLD_res
HCV.OLD_matrix_downsized
cond <- names(HCV.OLD_matrix_downsized)
names(HCV.OLD_matrix_downsized) <- as.character(1:ncol(HCV.OLD_matrix_downsized))
HCV.OLD_res <- DSrun(HCV.OLD_matrix_downsized, cond, 'Control')

#HCV.NEW_res
HCV.NEW_matrix <- NEW_matrix_downsized[,7:12]
names(HCV.NEW_matrix)
cond <- c(rep("HCV",3), rep("Control",3))
load("DSrun.Rdata")
HCV.NEW_res <- DSrun(HCV.NEW_matrix, cond, 'Control')
HCV.NEW_res_df <- data.frame(row=HCV.NEW_res@rownames, col=HCV.NEW_res@listData)
HCV.NEW_up <- (HCV.NEW_res_df %>% select(row,col.log2FoldChange, col.padj) %>% filter(!is.na(col.padj), col.padj<1e-5, col.log2FoldChange>0.5))$row
HCV.NEW_down <- (HCV.NEW_res_df %>% select(row,col.log2FoldChange, col.padj) %>% filter(!is.na(col.padj), col.padj<1e-5, col.log2FoldChange<c(-0.5)))$row
HCV.NEW_DE <- c(HCV.NEW_up, HCV.NEW_down)



#SARS_res
SARS_matrix <- NEW_matrix_downsized[,c(4:6,10:12)]
names(SARS_matrix)
cond <- c(rep("SARS",3), rep("Control",3))
SARS_res <- DSrun(SARS_matrix, cond, 'Control')
SARS_res_df <- data.frame(row=SARS_res@rownames, col=SARS_res@listData)
SARS_up <- (SARS_res_df %>% select(row,col.log2FoldChange, col.padj) %>% filter(!is.na(col.padj), col.padj<1e-5, col.log2FoldChange>0.5))$row
SARS_down <- (SARS_res_df %>% select(row,col.log2FoldChange, col.padj) %>% filter(!is.na(col.padj), col.padj<1e-5, col.log2FoldChange<c(-0.5)))$row
SARS_DE <- c(SARS_up, SARS_down)

#SARS.HCV_res
SARS.HCV_matrix <- NEW_matrix_downsized[,c(1:3,10:12)]
names(SARS.HCV_matrix)
cond <- c(rep("SARS.HCV",3), rep("Control",3))
SARS.HCV_res <- DSrun(SARS.HCV_matrix, cond, 'Control')
SARS.HCV_res_df <- data.frame(row=SARS.HCV_res@rownames, col=SARS.HCV_res@listData)
SARS.HCV_up <- (SARS.HCV_res_df %>% select(row,col.log2FoldChange, col.padj) %>% filter(!is.na(col.padj), col.padj<1e-5, col.log2FoldChange>0.5))$row
SARS.HCV_down <- (SARS.HCV_res_df %>% select(row,col.log2FoldChange, col.padj) %>% filter(!is.na(col.padj), col.padj<1e-5, col.log2FoldChange<c(-0.5)))$row
SARS.HCV_DE <- c(SARS.HCV_up, SARS.HCV_down)

#Volcano plots
keyvals <- ifelse(
  HCV.NEW_res$log2FoldChange < -0.5 & HCV.NEW_res$padj<1e-5, 'darkblue',
  ifelse(HCV.NEW_res$log2FoldChange > 0.5 & HCV.NEW_res$padj<1e-5, 'red3',
         'darkgrey'))
keyvals[is.na(keyvals)] <- 'darkgrey'
names(keyvals)[keyvals == 'red3'] <- 'high'
names(keyvals)[keyvals == 'darkgrey'] <- 'mid'
names(keyvals)[keyvals == 'darkblue'] <- 'low'
colAlpha <- ifelse(names(keyvals)[keyvals == 'darkgrey']=='mid',1,0)
pdf("HCV.NEW_DE_volcano.pdf")
EnhancedVolcano(HCV.NEW_res,
                lab = NA,
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-5, 5),
                ylim= c(0,300),
                title='',
                pCutoff = 1e-5,
                FCcutoff = 0.5,
                pointSize = 3,
                cutoffLineWidth = 0.8,
                colCustom = keyvals,
                legendPosition = 'none',
                gridlines.major = FALSE,
                gridlines.minor = FALSE, subtitle = '')#, colAlpha = colAlpha
dev.off()

# #Volcano plots
# keyvals <- ifelse(
#   HCV.NEW_res$log2FoldChange < -0.5 & HCV.NEW_res$padj<1e-5, 'darkblue',
#   ifelse(HCV.NEW_res$log2FoldChange > 0.5 & HCV.NEW_res$padj<1e-5, 'red3',
#          'grey'))
# keyvals[is.na(keyvals)] <- 'grey'
# names(keyvals)[keyvals == 'red3'] <- 'high'
# names(keyvals)[keyvals == 'grey'] <- 'mid'
# names(keyvals)[keyvals == 'darkblue'] <- 'low'
# colAlpha <- ifelse(rownames(SARS.HCV_res)%in%HCV.NEW_DE,1,0)
# pdf("SARS.HCV_DE_volcano.pdf")
# EnhancedVolcano(SARS.HCV_res,
#                 lab = NA,
#                 x = 'log2FoldChange',
#                 y = 'padj',
#                 xlim = c(-5, 5),
#                 title='',
#                 pCutoff = 1e-5,
#                 FCcutoff = 0.5,
#                 pointSize = 3,
#                 cutoffLineWidth = 0.8,
#                 #colCustom = keyvals,
#                 legendPosition = 'none',
#                 gridlines.major = FALSE,
#                 gridlines.minor = FALSE, subtitle = '')#, colAlpha = colAlpha
# dev.off()

#Genes overlap
HCV.NEW_up, SARS_up, HCV.NEW_up
length(HCV.NEW_down)
length(SARS_down)
length(SARS.HCV_down)

length(Reduce(intersect, list(HCV.NEW_down,SARS.HCV_down)))
SARS.HCV_unique <- setdiff(SARS.HCV_DE, c(HCV.NEW_DE,SARS_DE))
SARS.HCV_uniq_cond <- rownames(SARS.HCV_res)%in%SARS.HCV_unique & SARS.HCV_res$padj<1e-5 & SARS_res$padj>0.05 & HCV.NEW_res$padj>0.05

#Volcano plot of SARS.HCV_unique
keyvals <- ifelse(
  SARS.HCV_uniq_cond & SARS.HCV_res$log2FoldChange < -0.5 & SARS.HCV_res$padj<1e-5, 'darkblue',
  ifelse(SARS.HCV_uniq_cond & SARS.HCV_res$log2FoldChange > 0.5 & SARS.HCV_res$padj<1e-5, 'red3',
         'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'red3'] <- 'high'
names(keyvals)[keyvals == 'grey'] <- 'mid'
names(keyvals)[keyvals == 'darkblue'] <- 'low'
colAlpha <- ifelse(keyvals%in%c('red3','darkblue'),1,0)
SARS.HCV_res_filtered <- rownames(SARS.HCV_res)[keyvals%in%c('red3','darkblue')]
sum(SARS.HCV_res_filtered%in%SARS.HCV_down)

pdf("SARS.HCV_unique_HCV.NEW_DE_volcano.pdf")
EnhancedVolcano(HCV.NEW_res,
                lab = NA,
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-5, 5),
                ylim = c(0,50),
                title='',
                pCutoff = 1e-5,
                FCcutoff = 0.5,
                pointSize = 3,
                cutoffLineWidth = 0.8,
                colCustom = keyvals,
                legendPosition = 'none',
                gridlines.major = FALSE,
                gridlines.minor = FALSE, subtitle = '', colAlpha = colAlpha)#, colAlpha = colAlpha
dev.off()
sum(SARS.HCV_down%in%SARS.HCV_unique)
#
#
#KEGGlibrary pathway analysis
HCV.NEW_DE
SARS_DE
SARS.HCV_DE
SARS.HCV_unique

library(KEGGprofile)
library(openxlsx)
db <- read.csv("human_gene_id.csv")
head(db)

HCV.NEW_db <- as.character(db$GeneID[na.omit(match(HCV.NEW_DE, toupper(db$Symbol)))])
KEGG_HCV.NEW <- find_enriched_pathway(HCV.NEW_db, species = "hsa", download_latest=T)[[1]][,c(1,6)]
KEGG_HCV.NEW_diff_pathways <- as.character(KEGG_HCV.NEW$Pathway_Name)
KEGG_HCV.NEW_all_pathways <- find_enriched_pathway(HCV.NEW_db, returned_pvalue=1, returned_adjpvalue = 1, species = "hsa",download_latest=T)[[1]][,c(1,6)]

SARS_db <- as.character(db$GeneID[na.omit(match(SARS_DE, toupper(db$Symbol)))])
KEGG_SARS <- find_enriched_pathway(SARS_db, species = "hsa",returned_adjpvalue = 0.0001, download_latest=T)[[1]][,c(1,6)]
KEGG_SARS_diff_pathways <- as.character(KEGG_SARS$Pathway_Name)
KEGG_SARS_all_pathways <- find_enriched_pathway(SARS_db, returned_pvalue=1, returned_adjpvalue = 1, species = "hsa",download_latest=T)[[1]][,c(1,6)]

SARS.HCV_db <- as.character(db$GeneID[na.omit(match(SARS.HCV_DE, toupper(db$Symbol)))])
KEGG_SARS.HCV <- find_enriched_pathway(SARS.HCV_db, species = "hsa",returned_adjpvalue = 0.0001, download_latest=T)[[1]][,c(1,6)]
KEGG_SARS.HCV_diff_pathways <- as.character(KEGG_SARS.HCV$Pathway_Name)
KEGG_SARS.HCV_all_pathways <- find_enriched_pathway(SARS.HCV_db, returned_pvalue=1, returned_adjpvalue = 1, species = "hsa",download_latest=T)[[1]][,c(1,6)]

SARS.HCV_unique_db <- as.character(db$GeneID[na.omit(match(SARS.HCV_res_filtered, toupper(db$Symbol)))])
KEGG_SARS.HCV_unique <- find_enriched_pathway(SARS.HCV_unique_db, species = "hsa", download_latest=T)[[1]][,c(1,6)]
KEGG_SARS.HCV_unique_diff_pathways <- as.character(KEGG_SARS.HCV_unique$Pathway_Name)
KEGG_SARS.HCV_unique_all_pathways <- find_enriched_pathway(SARS.HCV_unique_db, returned_pvalue=1, returned_adjpvalue = 1, species = "hsa",download_latest=T)[[1]][,c(1,6)]


Unique_DE_pathways <- unique(c(KEGG_HCV.NEW_diff_pathways,KEGG_SARS_diff_pathways,
                               KEGG_SARS.HCV_diff_pathways, KEGG_SARS.HCV_unique_diff_pathways))
HCV.NEW_df <- KEGG_HCV.NEW_all_pathways[match(Unique_DE_pathways, KEGG_HCV.NEW_all_pathways$Pathway_Name),]
SARS_df <- KEGG_SARS_all_pathways[match(Unique_DE_pathways, KEGG_SARS_all_pathways$Pathway_Name),]
SARS.HCV_df <- KEGG_SARS.HCV_all_pathways[match(Unique_DE_pathways, KEGG_SARS.HCV_all_pathways$Pathway_Name),]
SARS.HCV_unique_df <- KEGG_SARS.HCV_unique_all_pathways[match(Unique_DE_pathways, KEGG_SARS.HCV_unique_all_pathways$Pathway_Name),]


merged_pathway_df <- data.frame(Pathway_name=Unique_DE_pathways, HCV.NEW=HCV.NEW_df$pvalueAdj, SARS=SARS_df$pvalueAdj,
                                SARS.HCV=SARS.HCV_df$pvalueAdj, SARS.HCV_unique=SARS.HCV_unique_df$pvalueAdj)
merged_pathway_df[,2:ncol(merged_pathway_df)] <- -log10(merged_pathway_df[,2:ncol(merged_pathway_df)])
merged_pathway_df[is.na(merged_pathway_df)] <- 0
library(reshape)
melted_pathway_df <- melt(merged_pathway_df, id.vars = "Pathway_name")
pdf("pathway_NEW_cut_1e-5.pdf")
ggplot(melted_pathway_df, aes(x = Pathway_name, y = value, fill=variable)) + 
  geom_bar(position="dodge", stat = "identity") + coord_flip() + theme_bw()+
  scale_y_continuous(name="") +
  scale_x_discrete(name="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(face="bold", size=10, angle=0, color = 'black'),
        axis.text.y = element_text(face="bold", size=9, angle=0, color = 'black'))
dev.off()

