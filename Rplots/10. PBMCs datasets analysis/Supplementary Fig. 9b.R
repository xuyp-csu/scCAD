
homedir = 'D:/CSU/Rare cell detection/Github/scCAD/'
workdir = paste0(homedir, 'Rplots/10. PBMCs datasets analysis/')  

library(ggplot2)
library(RColorBrewer)

load(paste0(workdir, "Supplementary Fig. 9b_source_data1_PBMC-test_umap.Rdata"))
scCAD_res = read.table(paste0(workdir, "Supplementary Fig. 9a_source_data2_PBMC-test_scCAD_res.txt"))[,1]
lab = read.table(paste0(workdir, "Supplementary Fig. 9a_source_data3_PBMC-test_lab.txt"))[,1]

umap = cbind(umap, lab)
colnames(umap) = c("UMAP_1", "UMAP_2", "label")
n <- length(unique(lab))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

set.seed(2023)
pal = sample(col_vector, length(unique(lab)))

ggplot(umap, aes(x=UMAP_1, y=UMAP_2, color=label)) + 
  geom_point(size = 0.8) +
  theme_bw() +
  xlab("Umap-1") +
  ylab("Umap-2") +
  scale_color_manual(values = pal) +
  theme(text=element_text(size=15,  family="serif",face = 'bold'),
        strip.text = element_text(size=15, family="serif",face = 'bold'),
        legend.text = element_text(size=15, family="serif",face = 'bold'),
        legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 4)))
ggsave(path = workdir, filename = "Supplementary Fig. 9b.png",width = 12, height = 6, device='png', dpi=300)

temp_umap = umap
temp_umap$label = factor(scCAD_res, levels = c("Abundant","R1", "R2", "R3"))

ggplot(temp_umap, aes(x=UMAP_1, y=UMAP_2, color=label)) + 
  geom_point(size = 0.8) +
  theme_bw() +
  xlab("tSNE-1") +
  ylab("tSNE-2") +
  scale_color_manual(values = c("#BAB0ACFF", paletteer_dynamic("cartography::multi.pal", 20))) +
  theme(text=element_text(size=15,  family="serif",face = 'bold'),
        strip.text = element_text(size=15, family="serif",face = 'bold'),
        legend.text = element_text(size=15, family="serif",face = 'bold'),
        legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 4)))
ggsave(path = workdir, filename = "Supplementary Fig. 9b_Rlab.png",width = 8, height = 6, device='png', dpi=300)



