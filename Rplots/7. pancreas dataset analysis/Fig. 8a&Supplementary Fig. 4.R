
homedir = 'D:/CSU/Rare cell detection/Github/scCAD/'
workdir = paste0(homedir, 'Rplots/7. pancreas dataset analysis/')  

library(ggplot2)
library(paletteer)
library(reshape2)
library(RColorBrewer)

load(paste0(workdir, "Fig. 8a_source_data1_tsne.Rdata"))
lab = read.table(paste0(workdir, "Fig. 8a_source_data2_lab.txt"))[,1]
n <- length(unique(lab))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
set.seed(2023)
pal = sample(col_vector, length(unique(lab)))

ggplot(tsne_out, aes(x=tSNE_1, y=tSNE_2, color=label)) + 
  geom_point(size = 1) +
  theme_bw() +
  xlab("tSNE-1") +
  ylab("tSNE-2") +
  scale_color_manual(values = pal) +
  theme(text=element_text(size=20,  family="serif",face = 'bold'),
        strip.text = element_text(size=20, family="serif",face = 'bold'),
        legend.text = element_text(size=20, family="serif",face = 'bold'),
        legend.title = element_blank(),
        legend.position = "right",
        legend.key.width=unit(1,"cm"),
        legend.key.height=unit(0.5,"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(ncol =1, override.aes = list(size = 4)))
ggsave(path = workdir, filename = "Fig. 8a",width = 12, height = 8, device='png', dpi=300)

##################### rare cell type ###########################

rare_h = 0.01*length(lab)
rare_cell_lab = rep("Abundant", length(lab))
n_p = c()
for (l in names(which(table(lab)<rare_h))) {
  n_p = c(n_p, pal[which(sort(as.character(unique(lab)))==l)])
  rare_cell_lab[which(lab==l)] = l
}

temp_tsne_out = tsne_out
temp_tsne_out$label = rare_cell_lab
ggplot(temp_tsne_out, aes(x=tSNE_1, y=tSNE_2, color=label)) + 
  geom_point(size = 1) +
  theme_bw() +
  xlab("tSNE-1") +
  ylab("tSNE-2") +
  scale_color_manual(values = c("#BAB0ACFF", n_p)) +
  theme(text=element_text(size=20,  family="serif",face = 'bold'),
        strip.text = element_text(size=20, family="serif",face = 'bold'),
        legend.text = element_text(size=20, family="serif",face = 'bold'),
        legend.title = element_blank(),
        legend.position = "right",
        legend.key.width=unit(1,"cm"),
        legend.key.height=unit(0.5,"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(ncol =1, override.aes = list(size = 4)))
ggsave(path = workdir, filename = "Fig. 8a_Rlab.png",width = 12, height = 8, device='png', dpi=300)

##################### scCAD ###########################

scCAD_res = read.table(paste0(workdir, "Fig. 8a_source_data3_scCAD_res.txt"))[,1]
temp_tsne_out = tsne_out
temp_tsne_out$label = scCAD_res

ggplot(temp_tsne_out, aes(x=tSNE_1, y=tSNE_2, color=label)) + 
  geom_point(size = 1) +
  theme_bw() +
  xlab("tSNE-1") +
  ylab("tSNE-2") +
  scale_color_manual(values = c("#BAB0ACFF", paletteer_d("ggthemes::Tableau_10"))) +
  theme(text=element_text(size=20,  family="serif",face = 'bold'),
        strip.text = element_text(size=20, family="serif",face = 'bold'),
        legend.text = element_text(size=20, family="serif",face = 'bold'),
        legend.title = element_blank(),
        legend.position = "right",
        legend.key.width=unit(1,"cm"),
        legend.key.height=unit(0.5,"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(ncol =1, override.aes = list(size = 4)))
ggsave(path = workdir, filename = "Fig. 8a_scCAD_res.png",width = 12, height = 8, device='png', dpi=300)

##################### CellSIUS ###########################

CellSIUS_res = read.table(paste0(workdir, "Supplementary Fig. 4_source_data1_CellSIUS_res.txt"))[,1]
temp_tsne_out = tsne_out
temp_tsne_out$label = CellSIUS_res

ggplot(temp_tsne_out, aes(x=tSNE_1, y=tSNE_2, color=label)) + 
  geom_point(size = 1) +
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
  guides(color = guide_legend(ncol =1, override.aes = list(size = 4)))
ggsave(path = workdir, filename = "Supplementary Fig. 4_CellSIUS_res.png",width = 8, height = 6, device='png', dpi=300)

##################### EDGE ###########################

EDGE_res = read.table(paste0(workdir, "Supplementary Fig. 4_source_data2_EDGE_res.txt"))[,1]
temp_tsne_out = tsne_out
temp_tsne_out$label = EDGE_res

ggplot(temp_tsne_out, aes(x=tSNE_1, y=tSNE_2, color=label)) + 
  geom_point(size = 1) +
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
  guides(color = guide_legend(ncol =1, override.aes = list(size = 4)))
ggsave(path = workdir, filename = "Supplementary Fig. 4_EDGE_res.png",width = 8, height = 6, device='png', dpi=300)

##################### FiRE ###########################

FiRE_res = read.table(paste0(workdir, "Supplementary Fig. 4_source_data3_FiRE_res.txt"))[,1]
temp_tsne_out = tsne_out
temp_tsne_out$label = FiRE_res

ggplot(temp_tsne_out, aes(x=tSNE_1, y=tSNE_2, color=label)) + 
  geom_point(size = 1) +
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
  guides(color = guide_legend(ncol =1, override.aes = list(size = 4)))
ggsave(path = workdir, filename = "Supplementary Fig. 4_FiRE_res.png",width = 8, height = 6, device='png', dpi=300)

##################### GapClust ###########################

GapClust_res = read.table(paste0(workdir, "Supplementary Fig. 4_source_data4_GapClust_res.txt"))[,1]
temp_tsne_out = tsne_out
temp_tsne_out$label = GapClust_res

ggplot(temp_tsne_out, aes(x=tSNE_1, y=tSNE_2, color=label)) + 
  geom_point(size = 1) +
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
  guides(color = guide_legend(ncol =1, override.aes = list(size = 4)))
ggsave(path = workdir, filename = "Supplementary Fig. 4_GapClust_res.png",width = 8, height = 6, device='png', dpi=300)

##################### GiniClust3 ###########################

GiniClust3_res = read.table(paste0(workdir, "Supplementary Fig. 4_source_data5_GiniClust3_res.txt"))[,1]
temp_tsne_out = tsne_out
temp_tsne_out$label = GiniClust3_res

ggplot(temp_tsne_out, aes(x=tSNE_1, y=tSNE_2, color=label)) + 
  geom_point(size = 1) +
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
  guides(color = guide_legend(ncol =1, override.aes = list(size = 4)))
ggsave(path = workdir, filename = "Supplementary Fig. 4_GiniClust3_res.png",width = 8, height = 6, device='png', dpi=300)

##################### RaceID3 ###########################

RaceID3_res = read.table(paste0(workdir, "Supplementary Fig. 4_source_data6_RaceID3_res.txt"))[,1]
temp_tsne_out = tsne_out
temp_tsne_out$label = RaceID3_res

ggplot(temp_tsne_out, aes(x=tSNE_1, y=tSNE_2, color=label)) + 
  geom_point(size = 1) +
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
  guides(color = guide_legend(ncol =1, override.aes = list(size = 4)))
ggsave(path = workdir, filename = "Supplementary Fig. 4_RaceID3_res.png",width = 8, height = 6, device='png', dpi=300)

##################### SCA ###########################

SCA_res = read.table(paste0(workdir, "Supplementary Fig. 4_source_data7_SCA_res.txt"))[,1]
temp_tsne_out = tsne_out
temp_tsne_out$label = SCA_res

ggplot(temp_tsne_out, aes(x=tSNE_1, y=tSNE_2, color=label)) + 
  geom_point(size = 1) +
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
  guides(color = guide_legend(ncol =1, override.aes = list(size = 4)))
ggsave(path = workdir, filename = "Supplementary Fig. 4_SCA_res.png",width = 8, height = 6, device='png', dpi=300)

##################### SCISSORS ###########################

SCISSORS_res = read.table(paste0(workdir, "Supplementary Fig. 4_source_data8_SCISSORS_res.txt"))[,1]
temp_tsne_out = tsne_out
temp_tsne_out$label = SCISSORS_res

ggplot(temp_tsne_out, aes(x=tSNE_1, y=tSNE_2, color=label)) + 
  geom_point(size = 1) +
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
  guides(color = guide_legend(ncol =1, override.aes = list(size = 4)))
ggsave(path = workdir, filename = "Supplementary Fig. 4_SCISSORS_res.png",width = 8, height = 6, device='png', dpi=300)


