
homedir = 'D:/CSU/Rare cell detection/Github/scCAD/'
workdir = paste0(homedir, 'Rplots/12. Parameters test/')  

library(paletteer)
library(ggplot2)

load(paste0(workdir, "Supplementary Fig. 11_source_data.Rdata"))

ggplot(df, aes(x=Dataset, y=Independence_score)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  theme_bw() +
  xlab("Dataset") +
  ylab("Independence score") +
  geom_hline(yintercept=0.7, linetype="dashed", color = "red") +
  theme(text=element_text(size=20,  family="serif",face = 'bold'),
        axis.title.x = element_text(size=12, family="serif",face = 'bold'),
        axis.text.x = element_text(size=12, family="serif",face = 'bold'),
        strip.text = element_text(size=20, family="serif",face = 'bold'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(path = workdir, filename = "Supplementary Fig. 11.png",width = 10, height = 6, device='png', dpi=300)

