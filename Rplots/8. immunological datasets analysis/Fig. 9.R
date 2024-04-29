
homedir = 'D:/CSU/Rare cell detection/Github/scCAD/'
workdir = paste0(homedir, 'Rplots/8. immunological datasets analysis/')  

library(ggplot2)
library(paletteer)
library(reshape2)

#######################  T cells  ##########################

load(paste0(workdir,"Fig. 9_source_data1_Tcells_tsne.Rdata"))
lab = read.table(paste0(workdir, "Fig. 9_source_data2_Tcells_lab.txt"))[,1]

tsne_out = cbind(tsne_out, lab)
colnames(tsne_out) = c("tSNE_1", "tSNE_2", "label")

ggplot(tsne_out, aes(x=tSNE_1, y=tSNE_2, color=label)) + 
  geom_point(size = 0.3) +
  theme_bw() +
  xlab("tSNE-1") +
  ylab("tSNE-2") +
  theme(text=element_text(size=15,  family="serif",face = 'bold'),
        strip.text = element_text(size=15, family="serif",face = 'bold'),
        legend.text = element_text(size=15, family="serif",face = 'bold'),
        legend.title = element_blank(),
        legend.position = "right",
        legend.key.width=unit(1,"cm"),
        legend.key.height=unit(0.5,"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(ncol =1, override.aes = list(size = 4)))
ggsave(path = workdir, filename = "Fig. 9a_Tcells.png",width = 12, height = 8, device='png', dpi=300)

temp_tsne_out = tsne_out
temp_tsne_out$label[which(temp_tsne_out$label!="CD8 Proliferating" &
                            temp_tsne_out$label!="CD4 Proliferating" &
                            temp_tsne_out$label!="dnT")] = "Abundant"
ggplot(temp_tsne_out, aes(x=tSNE_1, y=tSNE_2, color=label)) + 
  geom_point(size = 0.3) +
  theme_bw() +
  xlab("tSNE-1") +
  ylab("tSNE-2") +
  scale_color_manual(values = c("#BAB0ACFF", paletteer_d("ggthemes::Tableau_10"))) +
  theme(text=element_text(size=15,  family="serif",face = 'bold'),
        strip.text = element_text(size=15, family="serif",face = 'bold'),
        legend.text = element_text(size=15, family="serif",face = 'bold'),
        legend.title = element_blank(),
        legend.position = "right",
        legend.key.width=unit(1,"cm"),
        legend.key.height=unit(0.5,"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(ncol =1, override.aes = list(size = 4)))
ggsave(path = workdir, filename = "Fig. 9b_Tcells.png",width = 12, height = 8, device='png', dpi=300)

scCAD_res = read.table(paste0(workdir, "Fig. 9_source_data3_Tcells_scCAD_res.txt"))[,1]
temp_tsne_out = tsne_out
temp_tsne_out$label = scCAD_res

ggplot(temp_tsne_out, aes(x=tSNE_1, y=tSNE_2, color=label)) + 
  geom_point(size = 0.3) +
  theme_bw() +
  xlab("tSNE-1") +
  ylab("tSNE-2") +
  scale_color_manual(values = c("#BAB0ACFF", paletteer_d("ggthemes::Tableau_10"))) +
  theme(text=element_text(size=15,  family="serif",face = 'bold'),
        strip.text = element_text(size=15, family="serif",face = 'bold'),
        legend.text = element_text(size=15, family="serif",face = 'bold'),
        legend.title = element_blank(),
        legend.position = "right",
        legend.key.width=unit(1,"cm"),
        legend.key.height=unit(0.5,"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(ncol =1, override.aes = list(size = 4)))
ggsave(path =  workdir, filename = "Fig. 9c_Tcells.png" ,width = 11, height = 8, device='png', dpi=300)

#######################  Crohn  ##########################

lab = read.table(paste0(workdir, "Fig. 9_source_data5_Crohn_lab.txt"))[,1]
meta = read.table(paste0(workdir, "Fig. 9_source_data4_Crohn_tsne.txt"), sep = "\t")
rownames(meta) = meta[,1]
meta = meta[-c(1,2), c(2,3)]
tsne_out = meta
tsne_out = cbind(tsne_out, lab)
colnames(tsne_out) = c("tSNE_1", "tSNE_2", "label")
tsne_out$tSNE_1 = as.numeric(tsne_out$tSNE_1)
tsne_out$tSNE_2 = as.numeric(tsne_out$tSNE_2)

ggplot(tsne_out, aes(x=tSNE_1, y=tSNE_2, color=label)) + 
  geom_point(size = 0.3) +
  theme_bw() +
  xlab("tSNE-1") +
  ylab("tSNE-2") +
  theme(text=element_text(size=15,  family="serif",face = 'bold'),
        strip.text = element_text(size=15, family="serif",face = 'bold'),
        legend.text = element_text(size=15, family="serif",face = 'bold'),
        legend.title = element_blank(),
        legend.position = "right",
        legend.key.width=unit(1,"cm"),
        legend.key.height=unit(0.5,"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(ncol =1, override.aes = list(size = 4)))
ggsave(path = workdir, filename = "Fig. 9a_Crohn.png" ,width = 12, height = 8, device='png', dpi=300)

temp_tsne_out = tsne_out
temp_tsne_out$label[which(temp_tsne_out$label!="Glial Cells" &
                            temp_tsne_out$label!="Pericytes" &
                            temp_tsne_out$label!="Smooth muscle cells" &
                            temp_tsne_out$label!="Lymphatics" &
                            temp_tsne_out$label!="Mast cells" &
                            temp_tsne_out$label!="CD36+ endothelial cells" &
                            temp_tsne_out$label!="DC2 CD206+" &
                            temp_tsne_out$label!="DC1" &
                            temp_tsne_out$label!="Plasmablasts" &
                            temp_tsne_out$label!="T (gd)" &
                            temp_tsne_out$label!="Fibroblasts" &
                            temp_tsne_out$label!="ACKR1+ endothelial cells")] = "Abundant"
ggplot(temp_tsne_out, aes(x=tSNE_1, y=tSNE_2, color=label)) + 
  geom_point(size = 0.3) +
  theme_bw() +
  xlab("tSNE-1") +
  ylab("tSNE-2") +
  scale_color_manual(values = c("#BAB0ACFF", paletteer_d("ggthemes::Tableau_10"), "#D5E4A2", "#D2AF81", "#78909C")) +
  theme(text=element_text(size=15,  family="serif",face = 'bold'),
        strip.text = element_text(size=15, family="serif",face = 'bold'),
        legend.text = element_text(size=15, family="serif",face = 'bold'),
        legend.title = element_blank(),
        legend.position = "right",
        legend.key.width=unit(1,"cm"),
        legend.key.height=unit(0.5,"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(ncol =1, override.aes = list(size = 4)))
ggsave(path = workdir, filename = "Fig. 9b_Crohn.png" ,width = 12, height = 8, device='png', dpi=300)

scCAD_res = read.table(paste0(workdir, "Fig. 9_source_data6_Crohn_scCAD_res.txt"))[,1]
temp_tsne_out = tsne_out
temp_tsne_out$label = scCAD_res

ggplot(temp_tsne_out, aes(x=tSNE_1, y=tSNE_2, color=label)) + 
  geom_point(size = 0.3) +
  theme_bw() +
  xlab("tSNE-1") +
  ylab("tSNE-2") +
  scale_color_manual(values = c("#BAB0ACFF", paletteer_d("ggthemes::Tableau_10"))) +
  theme(text=element_text(size=15,  family="serif",face = 'bold'),
        strip.text = element_text(size=15, family="serif",face = 'bold'),
        legend.text = element_text(size=15, family="serif",face = 'bold'),
        legend.title = element_blank(),
        legend.position = "right",
        legend.key.width=unit(1,"cm"),
        legend.key.height=unit(0.5,"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(ncol =1, override.aes = list(size = 4)))
ggsave(path = workdir, filename = "Fig. 9c_Crohn.png" ,width = 10, height = 8, device='png', dpi=300)


