
homedir = 'D:/CSU/Rare cell detection/Github/scCAD/'
workdir = paste0(homedir, 'Rplots/4. Airway dataset analysis/')  

library(ggplot2)
library(paletteer)
library(reshape2)
library(pheatmap)

load(paste0(workdir,"Fig. 5b_5c_source_data_degs.Rdata"))
load(paste0(workdir,"Fig. 5c_source_data_subdata.Rdata"))

meta_data  <- read.table(paste0(workdir, "Fig. 5a_source_data1_tsne.txt"), sep = "\t", header = T)
meta_data$mouse[which(meta_data$mouse=="WT_2")] = "WT_M2"
meta_data$mouse[which(meta_data$mouse=="Foxj1_GFP_M1")] = "Foxj1GFP_M1"
meta_data$mouse[which(meta_data$mouse=="Foxj1_GFP_M2")] = "Foxj1GFP_M2"
rownames(meta_data) = paste0(meta_data$mouse,"_",meta_data$barcode)
lab = meta_data$cluster
names(lab) = rownames(meta_data)
lab[which(lab=="Hot")] = "Ionocyte"
scCAD_res = read.table(paste0(workdir, "Fig. 5a_source_data2_scCAD_res.txt"))[,1]
names(scCAD_res) = rownames(meta_data)

tmp_lab = lab
sub_dat = sub_dat[, names(sort(tmp_lab))]
tmp_lab = tmp_lab[colnames(sub_dat)]
tmp_lab = tmp_lab[setdiff(names(tmp_lab),  names(which(scCAD_res=="R1")))]
rare_list = rep("R1", length(names(which(scCAD_res=="R1"))))
names(rare_list) = names(which(scCAD_res=="R1"))
tmp_lab = c(tmp_lab, rare_list)
sub_dat = sub_dat[, names(tmp_lab)]

id_list = c()
n = 300
for (t in unique(tmp_lab)) {
  idx = which(tmp_lab==t)
  if(length(idx)>n){
    set.seed(2023)
    idx = sample(idx, n, replace = FALSE)
  }
  id_list = c(id_list, idx)
}

ann_colors = list(
  Types = c(Basal=paletteer_d("ggthemes::Tableau_10")[1],
            Ciliated=paletteer_d("ggthemes::Tableau_10")[2],
            Club=paletteer_d("ggthemes::Tableau_10")[8],
            Goblet=paletteer_d("ggthemes::Tableau_10")[4],
            Ionocyte=paletteer_d("ggthemes::Tableau_10")[5],
            Neuroendocrine=paletteer_d("ggthemes::Tableau_10")[6],
            Tuft=paletteer_d("ggthemes::Tableau_10")[7],
            R1=paletteer_d("ggthemes::Tableau_10")[3]))

annotation = as.data.frame(tmp_lab[id_list])
colnames(annotation) = "Types"
annotation$`Types` <- factor(annotation$`Types`, levels = c("Basal", "Ciliated", "Club",
                                                            "Goblet", "Ionocyte", "Neuroendocrine",
                                                            "Tuft", "R1"))

png(paste0(workdir, "Fig. 5c"), width = 8, height = 12, units = 'in', res = 600)
pheatmap(sub_dat[,id_list], border_color = NA, cluster_cols = FALSE, 
         annotation_colors = ann_colors,
         annotation_col = annotation, show_colnames = FALSE, annotation_names_col = TRUE)
dev.off()

