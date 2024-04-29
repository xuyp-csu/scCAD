
homedir = 'D:/CSU/Rare cell detection/Github/scCAD/'
workdir = paste0(homedir, 'Rplots/3. Robustness test/')  

library(ggplot2)
library(paletteer)

load(paste0(workdir, "Fig. 4b_source_data.Rdata"))

df$n = factor(df$n, levels = unique(df$n))
df_ = matrix(nrow = 0,ncol = 3)
for (i in 1:length(unique(df$n))) {
  tmp = df[df$n==i,]
  for (m in c("scCAD", "GapClust", "FiRE", "GiniClust3")) {
    df_ = rbind(df_, c(m,i,mean(tmp[tmp$med==m,4])))
  }
}
colnames(df_) = c("med", "n", "f1score")
df_ = as.data.frame(df_)
df_$n= factor(df_$n, levels = unique(df_$n))
df_$f1score = as.numeric(df_$f1score)

ggplot(df_, aes(x=n, y=f1score, group=med, color=med)) +
  geom_line(size=1.5)+
  theme_bw()+
  scale_x_discrete(name = "# of differential genes", breaks = seq(0,length(unique(df_$n)),10)) +
  ylab(expression(F[1]~"score")) +
  ylim(c(0,1))+
  scale_color_manual(values = paletteer_d("ggthemes::Tableau_10")) +
  labs(colour="Methods") +
  theme(text=element_text(size=20,  family="serif",face = 'bold'),
        strip.text = element_text(size=20, family="serif",face = 'bold'),
        legend.text = element_text(size=15, family="serif",face = 'bold'),
        legend.title = element_text(size=15, family="serif",face = 'bold'),
        strip.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(path = homedir, filename = "Fig. 4b.png",width = 12, height = 6, device='png', dpi=300)


