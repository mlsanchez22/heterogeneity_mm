library(hipathia)

gtex_v8_norm <- readRDS("/mnt/lustre/scratch/CBRA/research/projects/heterogeneity_mm/results/00_mm_gtex/gtex_v8_norm.rds")

pathways <- load_pathways(species = "hsa")

trans_data <- translate_data(gtex_v8_norm, "hsa")
exp_data <- normalize_data(trans_data)
results <- hipathia(exp_data, pathways, decompose = FALSE, verbose=FALSE)
path_vals <- get_paths_data(results)
path_normalized <- normalize_paths(path_vals, pathways)

saveRDS(path_normalized, file="/mnt/lustre/scratch/CBRA/research/projects/heterogeneity_mm/results/00_mm_gtex/pathvals_gtexv8.rds")
