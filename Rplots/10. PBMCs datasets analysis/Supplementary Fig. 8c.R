
homedir = 'D:/CSU/Rare cell detection/Github/scCAD/'
workdir = paste0(homedir, 'Rplots/10. PBMCs datasets analysis/')  

library(pheatmap)

scCAD_res = read.table(paste0(workdir, "Supplementary Fig. 8a_source_data2_PBMC-bench-3_scCAD_res.txt"))[,1]
lab = read.table(paste0(workdir, "Supplementary Fig. 8a_source_data3_PBMC-bench-3_lab.txt"))[,1]
load(paste0(workdir, "Supplementary Fig. 8c_source_data1_PBMC-bench-3_degs.Rdata"))
load(paste0(workdir, "Supplementary Fig. 8c_source_data2_PBMC-bench-3_subdata.Rdata"))

#### R7

deg7 = rownames(degs[[7]])[which(degs[[7]]$avg_log2FC>=1.5)]
names(lab) = names(scCAD_res) = colnames(sub_dat)

idx = setdiff(colnames(sub_dat)[which(lab %in% c("Lymph prog"))], names(which(scCAD_res=="R7")))
tmp_tl = rep("R7", length(names(which(scCAD_res=="R7"))))
names(tmp_tl) = names(which(scCAD_res=="R7"))
tl = c(lab[idx], tmp_tl)
idx = c(idx, names(which(scCAD_res=="R7")))

annotation = as.data.frame(tl)
colnames(annotation) = "Types"
annotation$`Types` <- factor(annotation$`Types`, levels = c("Lymph prog",
                                                            "R7"))
ann_colors = list(
  Types = c("Lymph prog"=paletteer_d("ggthemes::Tableau_10")[1],
            "R7"=paletteer_d("ggthemes::Tableau_10")[2]))

png(paste0(workdir, "Supplementary Fig. 8c_R7.png"), width = 8, height = 12, units = 'in', res = 600)
pheatmap(sub_dat[,idx], border_color = NA, cluster_cols = FALSE, 
         annotation_colors = ann_colors,
         annotation_col = annotation, show_colnames = FALSE, annotation_names_col = TRUE)
dev.off()

#### R8

load(paste0(workdir, "Supplementary Fig. 8c_source_data3_PBMC-bench-3_subdata.Rdata"))

deg8 = rownames(degs[[8]])[which(degs[[8]]$avg_log2FC>=1.5)]

idx = setdiff(colnames(sub_dat)[which(lab %in% c("CD16+ Mono"))], names(which(scCAD_res=="R8")))
tmp_tl = rep("R8", length(names(which(scCAD_res=="R8"))))
names(tmp_tl) = names(which(scCAD_res=="R8"))
tl = c(lab[idx], tmp_tl)
idx = c(idx, names(which(scCAD_res=="R8")))

annotation = as.data.frame(tl)
colnames(annotation) = "Types"
annotation$`Types` <- factor(annotation$`Types`, levels = c("CD16+ Mono",
                                                            "R8"))
ann_colors = list(
  Types = c("CD16+ Mono"=paletteer_d("ggthemes::Tableau_10")[1],
            "R8"=paletteer_d("ggthemes::Tableau_10")[2]))

png(paste0(workdir, "Supplementary Fig. 8c_R8.png"), width = 8, height = 12, units = 'in', res = 600)
pheatmap(sub_dat[,idx], border_color = NA, cluster_cols = FALSE, 
         annotation_colors = ann_colors,
         annotation_col = annotation, show_colnames = FALSE, annotation_names_col = TRUE)
dev.off()




