
homedir = 'D:/CSU/Rare cell detection/Github/scCAD/'
workdir = paste0(homedir, 'Rplots/11. Retina & B_lymphoma datasets analysis/')  

library(ggplot2)
library(paletteer)
library(reshape2)

lab = read.table(paste0(workdir, "Fig. 11a_source_data1_lab.txt"))[,1]
scCAD_res = read.table(paste0(workdir, "Fig. 11a_source_data2_scCAD_res.txt"))[,1]
load(paste0(workdir, "Fig. 11a_source_data3_subdata.Rdata"))

tmp = sub_dat[, which(scCAD_res==paste0("R", 1))]
tmp = apply(tmp, 1, mean)
for (i in c(3,11,12,14,16)) {
  tmp = rbind(tmp, apply(sub_dat[, which(scCAD_res==paste0("R", i))], 1, mean))
}
for (t in paste0("BC",1:10)) {
  idx = which(lab==t)
  tmp = rbind(tmp, apply(sub_dat[, idx], 1, mean))
}
rownames(tmp) = c("R1","R3","R11","R12","R14","R16", paste0("BC",1:10))
png(paste0(workdir, "Fig. 11b.png"), width = 8, height = 6, units = 'in', res = 600)
pheatmap::pheatmap(cor(t(tmp))[1:6,7:16],fontsize = 20)
dev.off()


