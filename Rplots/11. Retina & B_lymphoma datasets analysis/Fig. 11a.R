
homedir = 'D:/CSU/Rare cell detection/Github/scCAD/'
workdir = paste0(homedir, 'Rplots/11. Retina & B_lymphoma datasets analysis/')  

library(ggplot2)
library(paletteer)
library(reshape2)

lab = read.table(paste0(workdir, "Fig. 11a_source_data1_lab.txt"))[,1]
scCAD_res = read.table(paste0(workdir, "Fig. 11a_source_data2_scCAD_res.txt"))[,1]
load(paste0(workdir, "Fig. 11a_source_data3_subdata.Rdata"))

ref = c("Cabp5","Cacna1s","Car10","Car8","Cck","Cntn4","Erbb4","Gabrr1","Glra1","Gng13","Gpr179","Grik1",
        "Grm6","Kcng4","Lhx4","Neto1","Neurod4","Scgn","St18","Syt2","Trpm1","Vsx1","Vsx2")

df = data.frame()
for (g in ref) {
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
  for (t in paste0("BC",1:10)) {
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
df = df[which(df$Clusters %in% c(paste0("R",1:(length(unique(scCAD_res))-1)), paste0("BC",1:10))),]
df = df[which(df$Clusters %in% c("R1","R3","R11","R12","R14","R16", paste0("BC",1:10))),]
ggplot(data = df, aes(x = Expressions, y = Clusters, fill = Genes)) +
  geom_violin(scale = 'width',
              trim = TRUE,
              width = 0.5) +
  theme_bw() +
  scale_fill_manual(values = c(paletteer_d("ggthemes::Tableau_10"), paletteer_dynamic("cartography::multi.pal", 20))) +
  scale_y_discrete(limits = rev(c("R1","R3","R11","R12","R14","R16", paste0("BC",1:10))),
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
ggsave(path = workdir, filename = "Fig. 11a.png",width = 12, height = 10, device='png', dpi=300)



