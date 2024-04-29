
homedir = 'D:/CSU/Rare cell detection/Github/scCAD/'
workdir = paste0(homedir, 'Rplots/9. ccRCC datasets analysis/')  

library(pheatmap)

lab = read.table(paste0(workdir, "Fig. 10b_source_data2_Kidney_ccRCC_lab.txt"))[,1]
scCAD_res = read.table(paste0(workdir, "Fig. 10b_source_data3_Kidney_ccRCC_scCAD_res.txt"))[,1]
load(paste0(workdir, "Fig. 10c_source_data1_Kidney_ccRCC_scCAD_degs.Rdata"))
load(paste0(workdir, "Fig. 10c_source_data2_Kidney_ccRCC_subdata.Rdata"))

cn = names(lab) = names(scCAD_res) = colnames(sub_dat)

id = cn[which(lab=="Tcell")]
R4 = names(which(scCAD_res=="R4"))
id = setdiff(id, R4)
tmp_lab = c(lab[id], rep("R4",length(R4)))
comb = c(id, R4)
names(tmp_lab) = comb

annotation = as.data.frame(tmp_lab)
colnames(annotation) = "Types"
annotation$`Types` <- factor(annotation$`Types`, levels = c("Tcell",
                                                            "R4"))
ann_colors = list(
  Types = c("Tcell" = as.character(paletteer_d("ggthemes::Tableau_10")[1]),
            "R4" = as.character(paletteer_d("ggthemes::Tableau_10")[2])))

png(paste0(workdir, "Fig. 10c_R4.png"), width = 8, height = 12, units = 'in', res = 600)
pheatmap(sub_dat[,comb], border_color = NA, cluster_cols = FALSE, 
         annotation_colors = ann_colors,
         annotation_col = annotation, show_colnames = FALSE, annotation_names_col = TRUE)
dev.off()

load(paste0(workdir, "Fig. 10c_source_data3_Kidney_ccRCC_subdata.Rdata"))

id = cn[which(lab=="Macro")]
R5 = names(which(scCAD_res=="R5"))
id = setdiff(id, R5)
tmp_lab = c(lab[id], rep("R5",length(R5)))
comb = c(id, R5)
names(tmp_lab) = comb

annotation = as.data.frame(tmp_lab)
colnames(annotation) = "Types"
annotation$`Types` <- factor(annotation$`Types`, levels = c("Macro",
                                                            "R5"))
ann_colors = list(
  Types = c("Macro" = as.character(paletteer_d("ggthemes::Tableau_10")[1]),
            "R5" = as.character(paletteer_d("ggthemes::Tableau_10")[2])))

png(paste0(workdir, "Fig. 10c_R5.png"), width = 8, height = 12, units = 'in', res = 600)
pheatmap(sub_dat[,comb], border_color = NA, cluster_cols = FALSE, 
         annotation_colors = ann_colors,
         annotation_col = annotation, show_colnames = FALSE, annotation_names_col = TRUE)
dev.off()

load(paste0(workdir, "Fig. 10c_source_data4_Kidney_ccRCC_subdata.Rdata"))

id = cn[which(lab=="ua")]
R7 = names(which(scCAD_res=="R7"))
id = setdiff(id, R7)
tmp_lab = c(lab[id], rep("R7",length(R7)))
comb = c(id, R7)
names(tmp_lab) = comb

annotation = as.data.frame(tmp_lab)
colnames(annotation) = "Types"
annotation$`Types` <- factor(annotation$`Types`, levels = c("ua",
                                                            "R7"))
ann_colors = list(
  Types = c("ua" = as.character(paletteer_d("ggthemes::Tableau_10")[1]),
            "R7" = as.character(paletteer_d("ggthemes::Tableau_10")[2])))

png(paste0(workdir, "Fig. 10c_R7.png"), width = 8, height = 12, units = 'in', res = 600)
pheatmap(sub_dat[,comb], border_color = NA, cluster_cols = FALSE, 
         annotation_colors = ann_colors,
         annotation_col = annotation, show_colnames = FALSE, annotation_names_col = TRUE)
dev.off()




