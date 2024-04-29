
homedir = 'D:/CSU/Rare cell detection/Github/scCAD/'
workdir = paste0(homedir, 'Rplots/4. Airway dataset analysis/')  

library(ggplot2)
library(paletteer)
library(reshape2)

meta_data  <- read.table(paste0(workdir, "Fig. 5a_source_data1_tsne.txt"), sep = "\t", header = T)
meta_data$mouse[which(meta_data$mouse=="WT_2")] = "WT_M2"
meta_data$mouse[which(meta_data$mouse=="Foxj1_GFP_M1")] = "Foxj1GFP_M1"
meta_data$mouse[which(meta_data$mouse=="Foxj1_GFP_M2")] = "Foxj1GFP_M2"
rownames(meta_data) = paste0(meta_data$mouse,"_",meta_data$barcode)
lab = meta_data$cluster
tsne_out = meta_data[,1:2]
tsne_out = as.data.frame(tsne_out)
tsne_out = cbind(tsne_out, lab)
colnames(tsne_out) = c("tSNE_1", "tSNE_2", "label")

ggplot(tsne_out, aes(x=tSNE_1, y=tSNE_2, color=label)) + 
  geom_point(size = 0.8) +
  theme_bw() +
  xlab("tSNE-1") +
  ylab("tSNE-2") +
  scale_color_manual(values = paletteer_d("ggthemes::Tableau_10")) +
  theme(text=element_text(size=15,  family="serif",face = 'bold'),
        strip.text = element_text(size=15, family="serif",face = 'bold'),
        legend.text = element_text(size=15, family="serif",face = 'bold'),
        legend.title = element_blank(),
        legend.position = c(.99, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.key.width=unit(1,"cm"),
        legend.key.height=unit(0.5,"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 4)))
ggsave(path = workdir, filename = "Fig. 5a.png",width = 8, height = 6, device='png', dpi=300)

##################### scCAD ###########################

scCAD_res = read.table(paste0(workdir, "Fig. 5a_source_data2_scCAD_res.txt"))[,1]
temp_tsne_out = tsne_out
temp_tsne_out$label = scCAD_res
ggplot(temp_tsne_out, aes(x=tSNE_1, y=tSNE_2, color=label)) + 
  geom_point(size = 0.8) +
  theme_bw() +
  xlab("tSNE-1") +
  ylab("tSNE-2") +
  scale_color_manual(values = c("#BAB0ACFF", paletteer_d("ggthemes::Tableau_10"))) +
  theme(text=element_text(size=15,  family="serif",face = 'bold'),
        strip.text = element_text(size=15, family="serif",face = 'bold'),
        legend.text = element_text(size=15, family="serif",face = 'bold'),
        legend.title = element_blank(),
        legend.position = c(.99, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.key.width=unit(1,"cm"),
        legend.key.height=unit(0.5,"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 4)))
ggsave(path = workdir, filename = "Fig. 5a&Supplementary Fig. 2_scCAD_res.png",width = 8, height = 6, device='png', dpi=300)

##################### CellSIUS ###########################

CellSIUS_res = read.table(paste0(workdir, "Supplementary Fig. 2_source_data1_CellSIUS_res.txt"))[,1]
temp_tsne_out = tsne_out
temp_tsne_out$label = CellSIUS_res
ggplot(temp_tsne_out, aes(x=tSNE_1, y=tSNE_2, color=label)) + 
  geom_point(size = 0.8) +
  theme_bw() +
  xlab("tSNE-1") +
  ylab("tSNE-2") +
  scale_color_manual(values = c("#BAB0ACFF", "#E15759FF")) +
  theme(text=element_text(size=15,  family="serif",face = 'bold'),
        strip.text = element_text(size=15, family="serif",face = 'bold'),
        legend.text = element_text(size=15, family="serif",face = 'bold'),
        legend.title = element_blank(),
        legend.position = c(.99, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.key.width=unit(1,"cm"),
        legend.key.height=unit(0.5,"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 4)))
ggsave(path = workdir, filename = "Supplementary Fig. 2_CellSIUS_res.png",width = 8, height = 6, device='png', dpi=300)

##################### CIARA ###########################

CIARA_res = read.table(paste0(workdir, "Supplementary Fig. 2_source_data2_CIARA_res.txt"))[,1]
temp_tsne_out = tsne_out
temp_tsne_out$label = CIARA_res
ggplot(temp_tsne_out, aes(x=tSNE_1, y=tSNE_2, color=label)) + 
  geom_point(size = 0.8) +
  theme_bw() +
  xlab("tSNE-1") +
  ylab("tSNE-2") +
  scale_color_manual(values = c("#BAB0ACFF", "#E15759FF")) +
  theme(text=element_text(size=15,  family="serif",face = 'bold'),
        strip.text = element_text(size=15, family="serif",face = 'bold'),
        legend.text = element_text(size=15, family="serif",face = 'bold'),
        legend.title = element_blank(),
        legend.position = c(.99, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.key.width=unit(1,"cm"),
        legend.key.height=unit(0.5,"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 4)))
ggsave(path = workdir, filename = "Supplementary Fig. 2_CIARA_res.png",width = 8, height = 6, device='png', dpi=300)

##################### EDGE ###########################

EDGE_res = read.table(paste0(workdir, "Supplementary Fig. 2_source_data3_EDGE_res.txt"))[,1]
temp_tsne_out = tsne_out
temp_tsne_out$label = EDGE_res
ggplot(temp_tsne_out, aes(x=tSNE_1, y=tSNE_2, color=label)) + 
  geom_point(size = 0.8) +
  theme_bw() +
  xlab("tSNE-1") +
  ylab("tSNE-2") +
  scale_color_manual(values = c("#BAB0ACFF", "#E15759FF")) +
  theme(text=element_text(size=15,  family="serif",face = 'bold'),
        strip.text = element_text(size=15, family="serif",face = 'bold'),
        legend.text = element_text(size=15, family="serif",face = 'bold'),
        legend.title = element_blank(),
        legend.position = c(.99, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.key.width=unit(1,"cm"),
        legend.key.height=unit(0.5,"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 4)))
ggsave(path = workdir, filename = "Supplementary Fig. 2_EDGE_res.png",width = 8, height = 6, device='png', dpi=300)

##################### FiRE ###########################

FiRE_res = read.table(paste0(workdir, "Supplementary Fig. 2_source_data4_FiRE_res.txt"))[,1]
temp_tsne_out = tsne_out
temp_tsne_out$label = FiRE_res
ggplot(temp_tsne_out, aes(x=tSNE_1, y=tSNE_2, color=label)) + 
  geom_point(size = 0.8) +
  theme_bw() +
  xlab("tSNE-1") +
  ylab("tSNE-2") +
  scale_color_manual(values = c("#BAB0ACFF", "#E15759FF")) +
  theme(text=element_text(size=15,  family="serif",face = 'bold'),
        strip.text = element_text(size=15, family="serif",face = 'bold'),
        legend.text = element_text(size=15, family="serif",face = 'bold'),
        legend.title = element_blank(),
        legend.position = c(.99, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.key.width=unit(1,"cm"),
        legend.key.height=unit(0.5,"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 4)))
ggsave(path = workdir, filename = "Supplementary Fig. 2_FiRE_res.png",width = 8, height = 6, device='png', dpi=300)

##################### GapClust ###########################

GapClust_res = read.table(paste0(workdir, "Supplementary Fig. 2_source_data5_GapClust_res.txt"))[,1]
temp_tsne_out = tsne_out
temp_tsne_out$label = GapClust_res
ggplot(temp_tsne_out, aes(x=tSNE_1, y=tSNE_2, color=label)) + 
  geom_point(size = 0.8) +
  theme_bw() +
  xlab("tSNE-1") +
  ylab("tSNE-2") +
  scale_color_manual(values = c("#BAB0ACFF", "#E15759FF")) +
  theme(text=element_text(size=15,  family="serif",face = 'bold'),
        strip.text = element_text(size=15, family="serif",face = 'bold'),
        legend.text = element_text(size=15, family="serif",face = 'bold'),
        legend.title = element_blank(),
        legend.position = c(.99, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.key.width=unit(1,"cm"),
        legend.key.height=unit(0.5,"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 4)))
ggsave(path = workdir, filename = "Supplementary Fig. 2_GapClust_res.png",width = 8, height = 6, device='png', dpi=300)

##################### GiniClust3 ###########################

GiniClust3_res = read.table(paste0(workdir, "Supplementary Fig. 2_source_data6_GiniClust3_res.txt"))[,1]
temp_tsne_out = tsne_out
temp_tsne_out$label = GiniClust3_res
ggplot(temp_tsne_out, aes(x=tSNE_1, y=tSNE_2, color=label)) + 
  geom_point(size = 0.8) +
  theme_bw() +
  xlab("tSNE-1") +
  ylab("tSNE-2") +
  scale_color_manual(values = c("#BAB0ACFF", "#E15759FF")) +
  theme(text=element_text(size=15,  family="serif",face = 'bold'),
        strip.text = element_text(size=15, family="serif",face = 'bold'),
        legend.text = element_text(size=15, family="serif",face = 'bold'),
        legend.title = element_blank(),
        legend.position = c(.99, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.key.width=unit(1,"cm"),
        legend.key.height=unit(0.5,"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 4)))
ggsave(path = workdir, filename = "Supplementary Fig. 2_GiniClust3_res.png",width = 8, height = 6, device='png', dpi=300)

##################### SCA ###########################

SCA_res = read.table(paste0(workdir, "Supplementary Fig. 2_source_data7_SCA_res.txt"))[,1]
temp_tsne_out = tsne_out
temp_tsne_out$label = SCA_res
ggplot(temp_tsne_out, aes(x=tSNE_1, y=tSNE_2, color=label)) + 
  geom_point(size = 0.8) +
  theme_bw() +
  xlab("tSNE-1") +
  ylab("tSNE-2") +
  scale_color_manual(values = c("#BAB0ACFF", "#E15759FF")) +
  theme(text=element_text(size=15,  family="serif",face = 'bold'),
        strip.text = element_text(size=15, family="serif",face = 'bold'),
        legend.text = element_text(size=15, family="serif",face = 'bold'),
        legend.title = element_blank(),
        legend.position = c(.99, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.key.width=unit(1,"cm"),
        legend.key.height=unit(0.5,"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 4)))
ggsave(path = workdir, filename = "Supplementary Fig. 2_SCA_res.png",width = 8, height = 6, device='png', dpi=300)


