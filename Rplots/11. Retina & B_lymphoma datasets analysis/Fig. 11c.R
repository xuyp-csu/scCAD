
homedir = 'D:/CSU/Rare cell detection/Github/scCAD/'
workdir = paste0(homedir, 'Rplots/11. Retina & B_lymphoma datasets analysis/')  

library(paletteer)
library(readxl)
library(pheatmap)

scCAD_res = read.table(paste0(workdir, "Fig. 11a_source_data2_scCAD_res.txt"))[,1]
load(paste0(workdir, "Fig. 11c_source_data_subdata.Rdata"))

n_cells = ncol(sub_dat)
id = c(2,4,5,6,15,19)
tl = rep("p", n_cells)
names(tl) = colnames(sub_dat)
for (i in id) {
  tl[which(scCAD_res==paste0("R",i))] = paste0("R",i)
}
tl = tl[which(tl!="p")]

annotation = as.data.frame(tl)
colnames(annotation) = "Types"
annotation$`Types` <- factor(annotation$`Types`, levels = c("R2",
                                                            "R4",
                                                            "R5",
                                                            "R6",
                                                            "R15",
                                                            "R19"))
ann_colors = list(
  Types = c("R2"=paletteer_d("ggthemes::Tableau_10")[1],
            "R4"=paletteer_d("ggthemes::Tableau_10")[2],
            "R5"=paletteer_d("ggthemes::Tableau_10")[3],
            "R6"=paletteer_d("ggthemes::Tableau_10")[4],
            "R15"=paletteer_d("ggthemes::Tableau_10")[5],
            "R19"=paletteer_d("ggthemes::Tableau_10")[6]))

png(paste0(workdir, "Fig. 11c.png"), width = 8, height = 12, units = 'in', res = 600)
pheatmap(sub_dat[,names(sort(tl))], border_color = NA, cluster_cols = FALSE, 
         annotation_colors = ann_colors,
         annotation_col = annotation, show_colnames = FALSE, annotation_names_col = TRUE)
dev.off()



