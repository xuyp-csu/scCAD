
homedir = 'D:/CSU/Rare cell detection/Github/scCAD/'
workdir = paste0(homedir, 'Rplots/10. PBMCs datasets analysis/')  

library(pheatmap)

scCAD_res = read.table(paste0(workdir, "Supplementary Fig. 6a_source_data2_PBMC-bench-1_scCAD_res.txt"))[,1]
lab = read.table(paste0(workdir, "Supplementary Fig. 6a_source_data3_PBMC-bench-1_lab.txt"))[,1]
load(paste0(workdir, "Supplementary Fig. 6c_source_data1_PBMC-bench-1_degs.Rdata"))
load(paste0(workdir, "Supplementary Fig. 6c_source_data2_PBMC-bench-1_subdata.Rdata"))

deg4 = rownames(degs[[4]])[which(degs[[4]]$avg_log2FC>=1.5)]
names(lab) = names(scCAD_res) = colnames(sub_dat)

idx = setdiff(colnames(sub_dat)[which(lab %in% c("CD4+ T activated"))], names(which(scCAD_res=="R4")))
tmp_tl = rep("R4", length(which(scCAD_res=="R4")))
names(tmp_tl) = names(which(scCAD_res=="R4"))
tl = c(lab[idx], tmp_tl)
idx = c(idx, names(which(scCAD_res=="R4")))

annotation = as.data.frame(tl)
colnames(annotation) = "Types"
annotation$`Types` <- factor(annotation$`Types`, levels = c("CD4+ T activated",
                                                            "R4"))
ann_colors = list(
  Types = c("CD4+ T activated"=paletteer_d("ggthemes::Tableau_10")[1],
            "R4"=paletteer_d("ggthemes::Tableau_10")[2]))



png(paste0(workdir, "Supplementary Fig. 6c_R4.png"), width = 8, height = 12, units = 'in', res = 600)
pheatmap(sub_dat[,idx], border_color = NA, cluster_cols = FALSE, 
         annotation_colors = ann_colors,
         annotation_col = annotation, show_colnames = FALSE, annotation_names_col = TRUE)
dev.off()

#### R6

load(paste0(workdir, "Supplementary Fig. 6c_source_data3_PBMC-bench-1_subdata.Rdata"))

deg6 = rownames(degs[[6]])[which(degs[[6]]$avg_log2FC>=1.5)]

idx = setdiff(colnames(sub_dat)[which(lab %in% c("NK"))], names(which(scCAD_res=="R6")))
tmp_tl = rep("R6", length(names(which(scCAD_res=="R6"))))
names(tmp_tl) = names(which(scCAD_res=="R6"))
tl = c(lab[idx], tmp_tl)
idx = c(idx, names(which(scCAD_res=="R6")))

annotation = as.data.frame(tl)
colnames(annotation) = "Types"
annotation$`Types` <- factor(annotation$`Types`, levels = c("NK",
                                                            "R6"))
ann_colors = list(
  Types = c("NK"=paletteer_d("ggthemes::Tableau_10")[1],
            "R6"=paletteer_d("ggthemes::Tableau_10")[2]))

png(paste0(workdir, "Supplementary Fig. 6c_R6.png"), width = 8, height = 12, units = 'in', res = 600)
pheatmap(sub_dat[,idx], border_color = NA, cluster_cols = FALSE, 
         annotation_colors = ann_colors,
         annotation_col = annotation, show_colnames = FALSE, annotation_names_col = TRUE)
dev.off()

#### R10

load(paste0(workdir, "Supplementary Fig. 6c_source_data4_PBMC-bench-1_subdata.Rdata"))

deg10 = rownames(degs[[10]])[which(degs[[10]]$avg_log2FC>=1.5)]

idx = setdiff(colnames(sub_dat)[which(lab %in% c("B1 B"))], names(which(scCAD_res=="R10")))
tmp_tl = rep("R10", length(names(which(scCAD_res=="R10"))))
names(tmp_tl) = names(which(scCAD_res=="R10"))
tl = c(lab[idx], tmp_tl)
idx = c(idx, names(which(scCAD_res=="R10")))

annotation = as.data.frame(tl)
colnames(annotation) = "Types"
annotation$`Types` <- factor(annotation$`Types`, levels = c("B1 B",
                                                            "R10"))
ann_colors = list(
  Types = c("B1 B"=paletteer_d("ggthemes::Tableau_10")[1],
            "R10"=paletteer_d("ggthemes::Tableau_10")[2]))

png(paste0(workdir, "Supplementary Fig. 6c_R10.png"), width = 8, height = 12, units = 'in', res = 600)
pheatmap(sub_dat[,idx], border_color = NA, cluster_cols = FALSE, 
         annotation_colors = ann_colors,
         annotation_col = annotation, show_colnames = FALSE, annotation_names_col = TRUE)
dev.off()






