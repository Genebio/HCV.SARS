length(SARS.HCV_up)
length(SARS.HCV_down)
length(SARS_up)
length(SARS_down)
length(HCV.NEW_up)
length(HCV.NEW_down)
length(intersect(SARS.HCV_down, SARS_down))

length(intersect(SARS.HCV_up, HCV.NEW_down))/1154
length(Reduce(intersect, list(SARS.HCV_down, c(SARS_down))))*100/1686
library(openxlsx)
write.xlsx(SARS.HCV_up, file="SARS.HCV_up.xlsx")
write.xlsx(SARS_up, file="SARS_up.xlsx")
write.xlsx(HCV.NEW_up, file="HCV.NEW_up.xlsx")

write.xlsx(SARS.HCV_down, file="SARS.HCV_down.xlsx")
write.xlsx(SARS_down, file="SARS_down.xlsx")
write.xlsx(HCV.NEW_down, file="HCV.NEW_down.xlsx")

install.packages("BioVenn")
??BioVenn
library(BioVenn)
draw.venn(SARS.HCV_up,SARS_up, HCV.NEW_up, 
                     xtitle="SARS+HCV", ytitle="SARS", ztitle="HCV",output = "pdf",
                     filename = "BioVenn_upregulated.pdf",xt_s=1.5, yt_s=1.5, zt_s=1.5, nr_s=1.5,
                     subtitle = "", title="DE upregulated", nrtype = "abs")

draw.venn(SARS.HCV_down,SARS_down, HCV.NEW_down, 
          xtitle="SARS+HCV", ytitle="SARS", ztitle="HCV",output = "pdf",
          filename = "BioVenn_downregulated.pdf",xt_s=1.5, yt_s=1.5, zt_s=1.5, nr_s=1.5,
          subtitle = "", title="DE downregulated", nrtype = "abs")

#Venn overlaps of pathway names
KEGG_SARS.HCV_diff_pathways #37 pathways
KEGG_SARS_diff_pathways <- KEGG_SARS_diff_pathways[-8] #33 pathways
KEGG_HCV.NEW_diff_pathways <- KEGG_HCV.NEW_diff_pathways[-13] #46 pathways

venn_var <- draw.venn(KEGG_SARS.HCV_diff_pathways,KEGG_SARS_diff_pathways, KEGG_HCV.NEW_diff_pathways, 
          xtitle="SARS+HCV", ytitle="SARS", ztitle="HCV",output = "pdf",
          filename = "DE_pathways_perc_overlaps.pdf",xt_s=1.5, yt_s=1.5, zt_s=1.5, nr_s=1.5,
          subtitle = "", title="", nrtype = "abs")
venn_var$x_only
setdiff(c(venn_var$x_only), c(Pathway_names))
setdiff(venn_var$x_only, c(Pathway_names))
