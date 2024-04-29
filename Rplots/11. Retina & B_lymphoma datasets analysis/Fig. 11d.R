
homedir = 'D:/CSU/Rare cell detection/Github/scCAD/'
workdir = paste0(homedir, 'Rplots/11. Retina & B_lymphoma datasets analysis/')  

library(paletteer)
library(ggplot2)

lab = read.table(paste0(workdir, "Fig. 11d_source_data1_lab.txt"))[,1]
scCAD_res = read.table(paste0(workdir, "Fig. 11d_source_data2_scCAD_res.txt"))[,1]
load(paste0(workdir, "Fig. 11d_source_data3_tsne.Rdata"))

temp_tsne_out = cbind(tsne_out$Y, lab)
colnames(temp_tsne_out) = c("V1","V2","label")
temp_tsne_out = as.data.frame(temp_tsne_out)
temp_tsne_out$V1 = as.numeric(temp_tsne_out$V1)
temp_tsne_out$V2 = as.numeric(temp_tsne_out$V2)

ggplot(temp_tsne_out, aes(x=V1, y=V2, color=label)) + 
  geom_point(size = 0.8) +
  theme_bw() +
  xlab("tSNE-1") +
  ylab("tSNE-2") +
  scale_color_manual(values = c(paletteer_d("ggthemes::Tableau_10"), paletteer_dynamic("cartography::multi.pal", 20))) +
  theme(text=element_text(size=15,  family="serif",face = 'bold'),
        strip.text = element_text(size=15, family="serif",face = 'bold'),
        legend.text = element_text(size=15, family="serif",face = 'bold'),
        legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 4)))
ggsave(path = workdir, filename = "Fig. 11d.png",width = 10, height = 6, device='png', dpi=300)

temp_tsne_out = cbind(tsne_out$Y, scCAD_res)
colnames(temp_tsne_out) = c("V1","V2","label")
temp_tsne_out = as.data.frame(temp_tsne_out)
temp_tsne_out$label = factor(scCAD_res, levels = c("Abundant","R1", "R2", "R3", "R4", "R5"))
temp_tsne_out$V1 = as.numeric(temp_tsne_out$V1)
temp_tsne_out$V2 = as.numeric(temp_tsne_out$V2)

ggplot(temp_tsne_out, aes(x=V1, y=V2, color=label)) + 
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
ggsave(path = workdir, filename = "Fig. 11d_Rlab_tsne.png",width = 8, height = 6, device='png', dpi=300)





