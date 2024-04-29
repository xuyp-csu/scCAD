
homedir = 'D:/CSU/Rare cell detection/Github/scCAD/'
workdir = paste0(homedir, 'Rplots/10. PBMCs datasets analysis/')  

library(pheatmap)

scCAD_res = read.table(paste0(workdir, "Supplementary Fig. 9a_source_data2_PBMC-test_scCAD_res.txt"))[,1]
lab = read.table(paste0(workdir, "Supplementary Fig. 9a_source_data3_PBMC-test_lab.txt"))[,1]
load(paste0(workdir, "Supplementary Fig. 9c_source_data1_PBMC-test_degs.Rdata"))
load(paste0(workdir, "Supplementary Fig. 9c_source_data2_PBMC-test_subdata.Rdata"))

#### R1

deg1 = rownames(degs[[1]])[which(degs[[1]]$avg_log2FC>=1.5)]
names(lab) = names(scCAD_res) = colnames(sub_dat)

idx = setdiff(colnames(sub_dat)[which(lab %in% c("Lymph prog"))], names(which(scCAD_res=="R1")))
tmp_tl = rep("R1", length(names(which(scCAD_res=="R1"))))
names(tmp_tl) = names(which(scCAD_res=="R1"))
tl = c(lab[idx], tmp_tl)
idx = c(idx, names(which(scCAD_res=="R1")))

annotation = as.data.frame(tl)
colnames(annotation) = "Types"
annotation$`Types` <- factor(annotation$`Types`, levels = c("Lymph prog",
                                                            "R1"))
ann_colors = list(
  Types = c("Lymph prog"=paletteer_d("ggthemes::Tableau_10")[1],
            "R1"=paletteer_d("ggthemes::Tableau_10")[2]))

png(paste0(workdir, "Supplementary Fig. 9c_R1.png"), width = 8, height = 12, units = 'in', res = 600)
pheatmap(sub_dat[,idx], border_color = NA, cluster_cols = FALSE, 
         annotation_colors = ann_colors,
         annotation_col = annotation, show_colnames = FALSE, annotation_names_col = TRUE)
dev.off()

#### R3

load(paste0(workdir, "Supplementary Fig. 9c_source_data3_PBMC-test_subdata.Rdata"))

deg3 = rownames(degs[[3]])[which(degs[[3]]$avg_log2FC>=1.5)]

idx = setdiff(colnames(sub_dat)[which(lab %in% c("B1 B"))], names(which(scCAD_res=="R3")))
tmp_tl = rep("R3", length(names(which(scCAD_res=="R3"))))
names(tmp_tl) = names(which(scCAD_res=="R3"))
tl = c(lab[idx], tmp_tl)
idx = c(idx, names(which(scCAD_res=="R3")))

annotation = as.data.frame(tl)
colnames(annotation) = "Types"
annotation$`Types` <- factor(annotation$`Types`, levels = c("B1 B",
                                                            "R3"))
ann_colors = list(
  Types = c("B1 B"=paletteer_d("ggthemes::Tableau_10")[1],
            "R3"=paletteer_d("ggthemes::Tableau_10")[2]))

png(paste0(workdir, "Supplementary Fig. 9c_R3.png"), width = 8, height = 12, units = 'in', res = 600)
pheatmap(sub_dat[,idx], border_color = NA, cluster_cols = FALSE, 
         annotation_colors = ann_colors,
         annotation_col = annotation, show_colnames = FALSE, annotation_names_col = TRUE)
dev.off()




