
homedir = 'D:/CSU/Rare cell detection/Github/scCAD/'
workdir = paste0(homedir, 'Rplots/12. Parameters test/')  

library(paletteer)
library(ggplot2)

load(paste0(workdir, "Supplementary Fig. 10_source_data.Rdata"))

ggplot(df, aes(x = M, y = proportion, fill=type)) +
  geom_col(colour="black",width=0.5,    
           position=position_dodge(0.5)) +
  geom_line(aes(x=M, y=Num/(max(Num)/0.75), group=type),stat="identity",size=1)+
  theme_bw() +
  xlab("M") +
  ylab("Proportion") +
  scale_fill_manual(labels = c("All type", "Rare type"), values=paletteer_d("ggthemes::Tableau_10")) +
  scale_y_continuous(sec.axis=sec_axis(~.*(max(df$Num)/0.75),name="The number of clusters"))+
  theme(text=element_text(size=20,  family="serif",face = 'bold'),
        axis.title.x = element_text(size=12, family="serif",face = 'bold.italic'),
        axis.text.x = element_text(size=12, family="serif",face = 'bold'),
        strip.text = element_text(size=20, family="serif",face = 'bold'),
        legend.text = element_text(size=12, family="serif",face = 'bold'),
        legend.title = element_blank(),
        legend.position = c(.8, .95),
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.margin = margin(0, 0, 0, 0),
        legend.key.width=unit(0.3,"cm"),
        legend.key.height=unit(0.3,"cm"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(path = workdir, filename = "Supplementary Fig. 10.png",width = 10, height = 6, device='png', dpi=300)


