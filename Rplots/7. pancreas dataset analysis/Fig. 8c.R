
homedir = 'D:/CSU/Rare cell detection/Github/scCAD/'
workdir = paste0(homedir, 'Rplots/7. pancreas dataset analysis/')  

library(ggplot2)
library(paletteer)
library(reshape2)
library(RColorBrewer)

lab = read.table(paste0(workdir, "Fig. 8a_source_data2_lab.txt"))[,1]
scCAD_res = read.table(paste0(workdir, "Fig. 8a_source_data3_scCAD_res.txt"))[,1]
load(paste0(workdir, "Fig. 8b_source_data1_degs.Rdata"))
load(paste0(workdir, "Fig. 8c_source_data_subdata.Rdata"))

tmp_lab = lab
sub_dat = sub_dat[, order(tmp_lab)]
scCAD_res = scCAD_res[order(tmp_lab)]
tmp_lab = tmp_lab[order(tmp_lab)]
names(scCAD_res) = names(tmp_lab) = colnames(sub_dat)

tmp_lab = tmp_lab[setdiff(names(tmp_lab), names(which(scCAD_res=="R1")))]
rare_list = rep("R1", length(which(scCAD_res=="R1")))
names(rare_list) = names(which(scCAD_res=="R1"))
tmp_lab = c(tmp_lab, rare_list)
sub_dat = sub_dat[, names(tmp_lab)]
id_list = c(which(tmp_lab=="beta"), which(tmp_lab=="R1"))
ann_colors = list(
  Types = c(beta="#BAB0ACFF", 
            R1="#E15759FF"))
annotation = as.data.frame(tmp_lab[id_list])
colnames(annotation) = "Types"
annotation$`Types` <- factor(annotation$`Types`, levels = c("beta", "R1"))

png(paste0(workdir, "Fig. 8c.png"), width = 8, height = 8, units = 'in', res = 600)
pheatmap(sub_dat[,id_list], border_color = NA, cluster_cols = FALSE, 
         annotation_colors = ann_colors,
         annotation_col = annotation, show_colnames = FALSE, annotation_names_col = TRUE)
dev.off()

