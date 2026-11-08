# ==============================================================================
# PaaSc 复现工程总调度入口 (One-Click Reproduction Runner)
#
# 使用方式:
#   Rscript repro_scripts/run_all.R              # 默认 quick 模式 (约 2-3 分钟)
#   Rscript repro_scripts/run_all.R --mode quick  # 快速模式: 运行 01, 02, 04, 05, 06
#   Rscript repro_scripts/run_all.R --mode full   # 全量模式: 包含 28,871 细胞 03 重算 (约 15 分钟)
# ==============================================================================

suppressPackageStartupMessages({
  library(utils)
})

args <- commandArgs(trailingOnly = TRUE)
mode <- "quick"
if ("--mode" %in% args) {
  idx <- which(args == "--mode")
  if (idx < length(args)) mode <- args[idx + 1]
} else if (any(grepl("^--mode=", args))) {
  mode <- sub("^--mode=", "", args[grep("^--mode=", args)][1])
}

cat("===================================================================\n")
cat("  PaaSc (PLOS Comput Biol, 2025-11) 全流程自动化复现总调度器\n")
cat("  运行模式:", toupper(mode), if (mode == "quick") "(跳过12分钟的Fig6全量重算)" else "(包含全部2.8万细胞重算)", "\n")
cat("===================================================================\n\n")

base_dir <- "d:/project/gemini/reproduction/PaaSc"
setwd(base_dir)

scripts_to_run <- c(
  "01_pbmc3k_repro.R",
  "02_spatial_stxBrain_repro.R",
  if (mode == "full") "03_fig6_ifn_repro.R" else character(0),
  "04_fig8_st_cortex_repro.R",
  "05_fig3_fig4_fig5_fig7_fig9_benchmarks.R",
  "06_fig2_and_supfig2_repro.R"
)

results_status <- list()

for (s in scripts_to_run) {
  s_path <- file.path("repro_scripts", s)
  cat(sprintf("\n>>> 正在执行: %s ...\n", s))
  t_start <- proc.time()
  
  res <- tryCatch({
    system2("Rscript", args = c(s_path), stdout = "", stderr = "")
  }, error = function(e) {
    cat("执行出错:", conditionMessage(e), "\n")
    return(1)
  })
  
  elapsed <- round((proc.time() - t_start)[3], 1)
  status_str <- if (res == 0) "SUCCESS" else "FAILED"
  cat(sprintf("<<< 完成: %s | 状态: %s | 耗时: %s 秒\n", s, status_str, elapsed))
  results_status[[s]] <- data.frame(Script = s, Status = status_str, Elapsed_sec = elapsed)
}

if (mode == "quick") {
  cat("\n[提示] Quick 模式已跳过 03_fig6_ifn_repro.R (如需全量重算 28,871 细胞，请运行: Rscript repro_scripts/run_all.R --mode full)\n")
}

# ------------------------------------------------------------------------------
# 导出全环境依赖软件包版本清单 (repro_results/environment_versions.tsv)
# ------------------------------------------------------------------------------
cat("\n=== 收集并导出项目完整环境版本快照 ===\n")
pkgs <- c(
  "R", "PaaSc", "Seurat", "SeuratObject", "SeuratData", "stxBrain.SeuratData",
  "CelliD", "msigdbr", "ROCit", "pROC", "ComplexHeatmap", "circlize",
  "pheatmap", "ggh4x", "ggplot2", "patchwork", "tidyverse", "data.table",
  "mixtools", "Matrix", "sctransform", "uwot"
)

ver_list <- list()
for (p in pkgs) {
  v <- if (p == "R") {
    paste0(R.version$major, ".", R.version$minor)
  } else {
    tryCatch(as.character(packageVersion(p)), error = function(e) "Not Installed")
  }
  
  src <- if (p %in% c("R", "stats", "utils")) {
    "Base R"
  } else if (p %in% c("PaaSc")) {
    "GitHub (yoyoong/PaaSc v1.0.1)"
  } else if (p %in% c("CelliD", "ComplexHeatmap")) {
    "Bioconductor 3.23"
  } else if (p %in% c("SeuratData", "stxBrain.SeuratData")) {
    "GitHub / SatijaLab"
  } else {
    "CRAN"
  }
  ver_list[[p]] <- data.frame(Package = p, Version = v, Repository = src)
}
ver_df <- do.call(rbind, ver_list)
env_path <- file.path("repro_results", "environment_versions.tsv")
write.table(ver_df, env_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat("环境版本表已写入:", env_path, "\n")
print(ver_df, row.names = FALSE)

# 打印执行汇总
cat("\n===================================================================\n")
cat("  PaaSc 全部指定复现脚本执行完毕！\n")
cat("===================================================================\n")
if (length(results_status) > 0) {
  print(do.call(rbind, results_status), row.names = FALSE)
}
