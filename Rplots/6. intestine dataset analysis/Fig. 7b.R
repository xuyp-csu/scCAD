
homedir = 'D:/CSU/Rare cell detection/Github/scCAD/'
workdir = paste0(homedir, 'Rplots/6. intestine dataset analysis/')  

library(paletteer)
library(Rtsne)
library(pheatmap)

scCAD_res = read.table(paste0(workdir, "Fig. 7a_source_data2_scCAD_res.txt"))[,1]
load(paste0(workdir,"Fig. 7b_source_data1_scCAD_degs.Rdata"))
ref = read.table(paste0(workdir,"Fig. 7b_source_data2_degs_ref.txt"), fill = TRUE, sep = "\t")

ref_genes = c()
for (i in c(9,10,18,19)) {
  ref_genes = c(ref_genes, ref[1:10,i])
}
ref_genes = unique(ref_genes)
ref_genes = ref_genes[ref_genes!= ""]

load(paste0(workdir,"Fig. 7b_source_data3_subdata.Rdata"))
tmp_lab = scCAD_res[which(scCAD_res!="Abundant")]
sub_dat = sub_dat[, order(tmp_lab)]
tmp_lab = tmp_lab[order(tmp_lab)]
names(tmp_lab) = colnames(sub_dat)

ann_colors = list(
  Types = c(R1=paletteer_d("ggthemes::Tableau_10")[1],
            R2=paletteer_d("ggthemes::Tableau_10")[2],
            R3=paletteer_d("ggthemes::Tableau_10")[3],
            R4=paletteer_d("ggthemes::Tableau_10")[4],
            R5=paletteer_d("ggthemes::Tableau_10")[5],
            R6=paletteer_d("ggthemes::Tableau_10")[6]))

annotation = as.data.frame(tmp_lab)
colnames(annotation) = "Types"
annotation$`Types` <- factor(annotation$`Types`, levels = c("R1","R2","R3","R4","R5","R6"))

png(paste0(workdir, "Fig. 7b.png"), width = 8, height = 8, units = 'in', res = 600)
pheatmap(sub_dat, border_color = NA, cluster_col = FALSE, cluster_rows = FALSE,
         annotation_colors = ann_colors,
         annotation_col = annotation, show_colnames = FALSE, annotation_names_col = TRUE)
dev.off()


