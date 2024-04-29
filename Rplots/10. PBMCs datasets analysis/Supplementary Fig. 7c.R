
homedir = 'D:/CSU/Rare cell detection/Github/scCAD/'
workdir = paste0(homedir, 'Rplots/10. PBMCs datasets analysis/')  

library(pheatmap)

scCAD_res = read.table(paste0(workdir, "Supplementary Fig. 7a_source_data2_PBMC-bench-2_scCAD_res.txt"))[,1]
lab = read.table(paste0(workdir, "Supplementary Fig. 7a_source_data3_PBMC-bench-2_lab.txt"))[,1]
load(paste0(workdir, "Supplementary Fig. 7c_source_data1_PBMC-bench-2_degs.Rdata"))
load(paste0(workdir, "Supplementary Fig. 7c_source_data2_PBMC-bench-2_subdata.Rdata"))

#### R3

deg3 = rownames(degs[[3]])[which(degs[[3]]$avg_log2FC>=1.5)]
names(lab) = names(scCAD_res) = colnames(sub_dat)

idx = setdiff(colnames(sub_dat)[which(lab %in% c("CD4+ T activated"))], names(which(scCAD_res=="R3")))
tmp_tl = rep("R3", length(names(which(scCAD_res=="R3"))))
names(tmp_tl) = names(which(scCAD_res=="R3"))
tl = c(lab[idx], tmp_tl)
idx = c(idx, names(which(scCAD_res=="R3")))

annotation = as.data.frame(tl)
colnames(annotation) = "Types"
annotation$`Types` <- factor(annotation$`Types`, levels = c("CD4+ T activated",
                                                            "R3"))
ann_colors = list(
  Types = c("CD4+ T activated"=paletteer_d("ggthemes::Tableau_10")[1],
            "R3"=paletteer_d("ggthemes::Tableau_10")[2]))

png(paste0(workdir, "Supplementary Fig. 7c_R3.png"), width = 8, height = 12, units = 'in', res = 600)
pheatmap(sub_dat[,idx], border_color = NA, cluster_cols = FALSE, 
         annotation_colors = ann_colors,
         annotation_col = annotation, show_colnames = FALSE, annotation_names_col = TRUE)
dev.off()

#### R8

load(paste0(workdir, "Supplementary Fig. 7c_source_data3_PBMC-bench-2_subdata.Rdata"))

deg8 = rownames(degs[[8]])[which(degs[[8]]$avg_log2FC>=1.5)]

idx = setdiff(colnames(sub_dat)[which(lab %in% c("Lymph prog"))], names(which(scCAD_res=="R8")))
tmp_tl = rep("R8", length(names(which(scCAD_res=="R8"))))
names(tmp_tl) = names(which(scCAD_res=="R8"))
tl = c(lab[idx], tmp_tl)
idx = c(idx, names(which(scCAD_res=="R8")))

annotation = as.data.frame(tl)
colnames(annotation) = "Types"
annotation$`Types` <- factor(annotation$`Types`, levels = c("Lymph prog",
                                                            "R8"))
ann_colors = list(
  Types = c("Lymph prog"=paletteer_d("ggthemes::Tableau_10")[1],
            "R8"=paletteer_d("ggthemes::Tableau_10")[2]))

png(paste0(workdir, "Supplementary Fig. 7c_R8.png"), width = 8, height = 12, units = 'in', res = 600)
pheatmap(sub_dat[,idx], border_color = NA, cluster_cols = FALSE, 
         annotation_colors = ann_colors,
         annotation_col = annotation, show_colnames = FALSE, annotation_names_col = TRUE)
dev.off()

#### R9

load(paste0(workdir, "Supplementary Fig. 7c_source_data4_PBMC-bench-2_subdata.Rdata"))

deg9 = rownames(degs[[9]])[which(degs[[9]]$avg_log2FC>=1.5)]

idx = setdiff(colnames(sub_dat)[which(lab %in% c("B1 B"))], names(which(scCAD_res=="R9")))
tmp_tl = rep("R9", length(names(which(scCAD_res=="R9"))))
names(tmp_tl) = names(which(scCAD_res=="R9"))
tl = c(lab[idx], tmp_tl)
idx = c(idx, names(which(scCAD_res=="R9")))

annotation = as.data.frame(tl)
colnames(annotation) = "Types"
annotation$`Types` <- factor(annotation$`Types`, levels = c("B1 B",
                                                            "R9"))
ann_colors = list(
  Types = c("B1 B"=paletteer_d("ggthemes::Tableau_10")[1],
            "R9"=paletteer_d("ggthemes::Tableau_10")[2]))


png(paste0(workdir, "Supplementary Fig. 7c_R9.png"), width = 8, height = 12, units = 'in', res = 600)
pheatmap(sub_dat[,idx], border_color = NA, cluster_cols = FALSE, 
         annotation_colors = ann_colors,
         annotation_col = annotation, show_colnames = FALSE, annotation_names_col = TRUE)
dev.off()




