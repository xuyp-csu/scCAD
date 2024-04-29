
homedir = 'D:/CSU/Rare cell detection/Github/scCAD/'
workdir = paste0(homedir, 'Rplots/1. Benchmarking/')  

library(ggplot2)
library(paletteer)
library(readxl)
library(reshape2)

df = read_excel(paste0(workdir, "Fig. 2a_source_data.xlsx"))
df = as.data.frame(df)
rownames(df) = df[,1]
df = df[,-1]
df = melt(df, measure.vars = colnames(df))
colnames(df) = c("med","F1")
df$med = factor(df$med, levels = c("scCAD", setdiff(sort(unique(as.character(df$med))),"scCAD")))
df[which(df=="NA", arr.ind = T)]=0
df$F1 = as.numeric(df$F1)

ggplot(df, aes(x = med, y = F1, fill = med)) +
  geom_boxplot(width = 0.5, outlier.shape = NA)+
  theme_bw() +
  scale_y_continuous(name = expression(F[1]~"score")) +
  scale_fill_manual(values = c(paletteer_d("ggthemes::Tableau_10"), "#CAB2D6FF")) +
  # scale_fill_manual(values = paletteer_d("tidyquant::tq_dark")) +
  theme(text=element_text(size=20,  family="serif",face = 'bold'),
        strip.text = element_text(size=20, family="serif",face = 'bold'),
        legend.position = "none",
        axis.title.x=element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(path = workdir, filename = "Fig. 2a.png",width = 14, height = 6, device='png', dpi=300)



