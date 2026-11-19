# ==============================================================================
# PaaSc 复现工程总调度入口 (One-Click Reproduction Runner)
#
# 设计说明:
#   - 通过 Rscript --file= 自动解析 PaaSc/ 项目根目录，无硬编码绝对路径。
#   - 显式校验 --mode 参数 (quick / full)，前置校验 Zenodo 大文件存在性。
#   - 使用当前 R 运行时对应的 Rscript 可执行文件，任一子脚本退出码非 0 即熔断并退出非 0 状态。
#   - 自动收集全环境软件包版本写入 repro_results/environment_versions.tsv。
#
# 使用方式:
#   Rscript repro_scripts/run_all.R              # 默认 quick 模式 (约 2-3 分钟)
#   Rscript repro_scripts/run_all.R --mode quick # 快速模式: 运行 01, 02, 04, 05, 06
#   Rscript repro_scripts/run_all.R --mode full  # 全量模式: 包含 28,871 细胞 03 重算 (约 12-15 分钟)
# ==============================================================================

suppressPackageStartupMessages({
  library(utils)
})

# 解析并严格校验 --mode 参数 (支持 "--mode quick" 与 "--mode=quick")
args <- commandArgs(trailingOnly = TRUE)
mode <- "quick"
if ("--mode" %in% args) {
  idx <- which(args == "--mode")[1]
  if (idx < length(args)) {
    mode <- args[idx + 1]
  } else {
    stop("参数错误: --mode 后必须指定 quick 或 full。")
  }
} else if (any(grepl("^--mode=", args))) {
  mode <- sub("^--mode=", "", args[grep("^--mode=", args)][1])
}
mode <- tolower(trimws(mode))
if (!(mode %in% c("quick", "full"))) {
  stop(sprintf("非法运行模式: '%s'。仅支持 --mode quick 或 --mode full。", mode))
}

# 动态定位 PaaSc/ 项目根目录，彻底消除本机硬编码路径
get_project_root <- function() {
  c_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", c_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
    return(normalizePath(file.path(script_dir, ".."), mustWork = TRUE))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

base_dir <- get_project_root()
setwd(base_dir)

rscript_bin <- file.path(
  R.home("bin"),
  if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
)
if (!file.exists(rscript_bin)) rscript_bin <- "Rscript"

cat("===================================================================\n")
cat("  PaaSc (PLOS Comput Biol, 2025-11) 全流程自动化复现总调度器\n")
cat("  工作目录:", base_dir, "\n")
cat("  运行模式:", toupper(mode),
    if (mode == "quick") "(跳过约12分钟的Fig6全量重算)" else "(包含全部28,871细胞Fig6端到端重算)", "\n")
cat("===================================================================\n\n")

# 前置校验 .gitignore 屏蔽的大体积 Zenodo RDS 文件是否存在
required_rds <- c("zenodo_data/PaaSc_data/Fig8/ST_braincortex.RDS")
if (mode == "full") {
  required_rds <- c("zenodo_data/PaaSc_data/Fig6/IFN.RDS", required_rds)
}
missing_rds <- required_rds[!file.exists(required_rds)]
if (length(missing_rds) > 0) {
  stop(
    "缺少未纳入 Git 的大体积 Zenodo 数据文件:\n  - ",
    paste(missing_rds, collapse = "\n  - "),
    "\n请从 Zenodo (10.5281/zenodo.17138447 / 17136774) 下载 PaaSc_data.zip 并解压至 zenodo_data/PaaSc_data/。"
  )
}

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
  if (!file.exists(s_path)) {
    stop("找不到目标脚本: ", s_path)
  }
  cat(sprintf("\n>>> 正在执行: %s ...\n", s))
  t_start <- proc.time()

  res <- tryCatch({
    system2(rscript_bin, args = c(s_path), stdout = "", stderr = "")
  }, error = function(e) {
    cat("执行出错:", conditionMessage(e), "\n")
    return(1L)
  })

  elapsed <- round((proc.time() - t_start)[3], 1)
  status_str <- if (identical(as.integer(res), 0L)) "SUCCESS" else "FAILED"
  cat(sprintf("<<< 完成: %s | 状态: %s | 耗时: %s 秒\n", s, status_str, elapsed))
  results_status[[s]] <- data.frame(Script = s, Status = status_str, Elapsed_sec = elapsed, stringsAsFactors = FALSE)

  if (!identical(as.integer(res), 0L)) {
    cat("\n[FAIL-FAST] 流水线中止于脚本:", s, "(退出码:", res, ")\n")
    print(do.call(rbind, results_status), row.names = FALSE)
    quit(save = "no", status = 1L)
  }
}

if (mode == "quick") {
  cat("\n[提示] Quick 模式已跳过 03_fig6_ifn_repro.R (如需全量重算 28,871 细胞，请运行: Rscript repro_scripts/run_all.R --mode full)\n")
}

# ------------------------------------------------------------------------------
# 导出全环境依赖软件包版本清单 (repro_results/environment_versions.tsv)
# ------------------------------------------------------------------------------
cat("\n=== 收集并导出项目完整环境版本快照 ===\n")
dir.create("repro_results", recursive = TRUE, showWarnings = FALSE)

bioc_label <- tryCatch(
  paste0("Bioconductor ", as.character(BiocManager::version())),
  error = function(e) "Bioconductor"
)

pkgs <- c(
  "R", "PaaSc", "Seurat", "SeuratObject", "SeuratData", "stxBrain.SeuratData",
  "CelliD", "msigdbr", "ROCit", "pROC", "ComplexHeatmap", "circlize",
  "pheatmap", "ggh4x", "ggplot2", "patchwork", "tidyverse", "data.table",
  "mixtools", "Matrix", "sctransform", "uwot"
)

ver_list <- lapply(pkgs, function(p) {
  v <- if (p == "R") {
    paste0(R.version$major, ".", R.version$minor)
  } else {
    tryCatch(as.character(packageVersion(p)), error = function(e) "Not Installed")
  }
  src <- if (p %in% c("R", "stats", "utils")) {
    "Base R"
  } else if (p == "PaaSc") {
    paste0("GitHub (yoyoong/PaaSc v", v, ")")
  } else if (p %in% c("CelliD", "ComplexHeatmap")) {
    bioc_label
  } else if (p %in% c("SeuratData", "stxBrain.SeuratData")) {
    "GitHub / SatijaLab"
  } else {
    "CRAN"
  }
  data.frame(Package = p, Version = v, Repository = src, stringsAsFactors = FALSE)
})

ver_df <- do.call(rbind, ver_list)
env_path <- file.path("repro_results", "environment_versions.tsv")
write.table(ver_df, env_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat("环境版本表已写入:", env_path, "\n")
print(ver_df, row.names = FALSE)

cat("\n===================================================================\n")
cat("  PaaSc 全部指定复现脚本执行完毕！\n")
cat("===================================================================\n")
if (length(results_status) > 0) {
  print(do.call(rbind, results_status), row.names = FALSE)
}
