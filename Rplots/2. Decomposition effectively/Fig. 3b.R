
homedir = 'D:/CSU/Rare cell detection/Github/scCAD/'
workdir = paste0(homedir, 'Rplots/2. Decomposition effectively/')  

library(ggplot2)
library(paletteer)
library(readxl)
library(reshape2)
library(scales)

df = read_excel(paste0(workdir, "Fig. 3b_source_data.xlsx"))
df = as.data.frame(df)
df$...1 = c("10X_PBMC", "Airway", "Arc-ME",  "Cao",  "Chung", "Colliculus",  "Darmanis", 
                 "Deng", "Goolam",  "Hippocampus", "iLNs", "Koh", "Li", "Livers", "Macosko",
                 "MacParland", "Mammary", "Pancreas", "Pollen", "Puram", "Shekhar", "Tonsil",
                 "UUOkidney", "Yang", "Zelsel")
df = melt(df)
colnames(df) = c("ds", "med", "proportion")
df[which(df$proportion==0),3] = 0.001
ggplot(df, aes(x = ds, y = proportion, fill = med)) +
  geom_col(colour="black",width=0.5,    
           position=position_dodge(0.5)) +
  theme_bw() +
  xlab("") +
  ylab("Proportion") +
  scale_fill_manual(labels = c("I-clusters", "M-clusters"), values=paletteer_d("ggthemes::Tableau_10")) +
  theme(text=element_text(size=20,  family="serif",face = 'bold'),
        axis.text.x = element_text(size=12, family="serif",face = 'bold', angle = 90, hjust = 1),
        strip.text = element_text(size=20, family="serif",face = 'bold'),
        legend.text = element_text(size=12, family="serif",face = 'bold'),
        legend.title = element_blank(),
        legend.position = c(.01, .99),
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.margin = margin(0, 0, 0, 0),
        legend.key.width=unit(0.3,"cm"),
        legend.key.height=unit(0.3,"cm"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(path = workdir, filename = "Fig. 3b.png", width = 10, height = 6, device='png', dpi=300)

