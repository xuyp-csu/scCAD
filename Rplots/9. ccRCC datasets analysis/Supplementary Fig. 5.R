
homedir = 'D:/CSU/Rare cell detection/Github/scCAD/'
workdir = paste0(homedir, 'Rplots/9. ccRCC datasets analysis/')  

library(Seurat)

load(paste0(workdir, "Supplementary Fig. 5_source_data1_Kidney_ccRCC_subdata.Rdata"))
load(paste0(workdir, "Supplementary Fig. 5_source_data2_Kidney_ccRCC_umap_seurat.Rdata"))

obj <- CreateSeuratObject(counts = sub_dat, min.cells = 3)
obj@reductions = reductions
rownames(obj@reductions$umap@cell.embeddings) = colnames(sub_dat)

FeaturePlot(obj, features = c("GZMB","GZMH","GNLY", "CD1C", "CD207", "FCER1A", "AHSP","HBD","HEMGN"))
ggsave(path = workdir, filename = "Supplementary Fig. 5.png",width = 12, height = 8, device='png', dpi=300)


