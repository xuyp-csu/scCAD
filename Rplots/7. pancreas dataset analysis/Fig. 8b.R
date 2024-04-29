
homedir = 'D:/CSU/Rare cell detection/Github/scCAD/'
workdir = paste0(homedir, 'Rplots/7. pancreas dataset analysis/')  

library(ggplot2)
library(paletteer)
library(reshape2)
library(RColorBrewer)

lab = read.table(paste0(workdir, "Fig. 8a_source_data2_lab.txt"))[,1]
scCAD_res = read.table(paste0(workdir, "Fig. 8a_source_data3_scCAD_res.txt"))[,1]
load(paste0(workdir, "Fig. 8b_source_data1_degs.Rdata"))
load(paste0(workdir, "Fig. 8b_source_data2_subdata.Rdata"))

deg1 = rownames(degs[[1]])[which(degs[[1]]$avg_log2FC>=1.5)][1:10]
deg2 = rownames(degs[[2]])[which(degs[[2]]$avg_log2FC>=1.5)][1:10]
deg3 = rownames(degs[[3]])[which(degs[[3]]$avg_log2FC>=1.5)][1:10]
deg4 = rownames(degs[[4]])[which(degs[[4]]$avg_log2FC>=1.5)][1:10]

df = data.frame()
for (g in c(deg1,deg2,deg3,deg4)) {
  tmp_df = data.frame()
  clu = c()
  for (t in setdiff(unique(scCAD_res), "Abundant")) {
    idx = which(scCAD_res==t)
    tmp = as.data.frame(sub_dat[g,idx])
    rownames(tmp) = NULL
    colnames(tmp) = "Expression"
    tmp_df = rbind(tmp_df,tmp)
    clu = c(clu, rep(t, length(idx)))
  }
  for (t in unique(lab)) {
    idx = which(lab==t)
    tmp = as.data.frame(sub_dat[g,idx])
    rownames(tmp) = NULL
    colnames(tmp) = "Expression"
    tmp_df = rbind(tmp_df, tmp)
    clu = c(clu, rep(t, length(idx)))
  }
  
  tmp_df = cbind(tmp_df, clu, rep(g,nrow(tmp_df)))
  rownames(tmp_df) = NULL
  colnames(tmp_df) = c("Expressions","Clusters","Genes")
  
  df = rbind(df, tmp_df)
}

df$Genes <- factor(df$Genes, levels = unique(df$Genes))
tdf = rep("A", nrow(df))
tdf[which(df$Genes %in% deg2)] = "B"
tdf[which(df$Genes %in% deg3)] = "C"
tdf[which(df$Genes %in% deg4)] = "D"
df = cbind(df, tdf)
colnames(df)[4] = "gClusters"

ggplot(data = df, aes(x = Expressions, y = Clusters, fill = gClusters)) +
  geom_violin(scale = 'width',
              trim = TRUE,
              width = 0.5) +
  theme_bw() +
  scale_fill_manual(values = paletteer_d("ggthemes::Tableau_10")) +
  scale_y_discrete(limits = rev(c("R1", "R2", "R3","R4", names(table(lab)))),
                   name = "") +
  facet_grid(cols = vars(Genes), scales = 'free') +
  xlab("")+
  theme(text=element_text(size=20,  family="serif",face = 'bold'),
        strip.text = element_text(size=20, family="serif",face = 'bold'),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x.bottom = element_blank(),
        axis.ticks = element_blank(),
        strip.text.x = element_text(size=18,  family="serif",face = 'bold', angle = 90, hjust = 0),
        strip.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(path = workdir, filename = "Fig. 8b.png",width = 12, height = 8, device='png', dpi=300)


