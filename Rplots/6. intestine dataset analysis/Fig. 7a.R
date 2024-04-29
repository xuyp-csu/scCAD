
homedir = 'D:/CSU/Rare cell detection/Github/scCAD/'
workdir = paste0(homedir, 'Rplots/6. intestine dataset analysis/')  

library(ggplot2)
library(paletteer)
library(reshape2)
library(RColorBrewer)

load(paste0(workdir, "Fig. 7a_source_data1_tsne.Rdata"))
scCAD_res = read.table(paste0(workdir, "Fig. 7a_source_data2_scCAD_res.txt"))[,1]
tsne_out = cbind(tsne_out, scCAD_res)
colnames(tsne_out) = c("tSNE_1", "tSNE_2", "label")

ggplot(tsne_out, aes(x=tSNE_1, y=tSNE_2, color=label)) + 
  geom_point(size = 1) +
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
ggsave(path = workdir, filename = "Fig. 7a.png",width = 8, height = 6, device='png', dpi=300)


