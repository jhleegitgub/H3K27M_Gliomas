############################################################
# 08_program_modules_analysis.R
#
#  Aim3: cell-cycle & PDGFRA module in the K27M networks
#
#  (1) 전체 malignant 네트워크에서
#      - cell-cycle gene / PDGFRA module gene 의
#        * 네트워크 포함 개수
#        * 평균 degree
#        * driver 비율
#        을 계산
#
#  (2) 그룹별(OPC / AC / OC) 네트워크에서도
#      동일한 summary를 계산
#
# 사용 파일:
#   data/raw/cellcycle_genes.txt
#   data/raw/PDGFRA_module_genes.txt
#
#   results/tables/K27M_netctrl_gene_summary_thr0p7.tsv
#
#   data/processed/K27M_group-<grp>_gene_degrees_thr0p7.tsv
#   data/processed/K27M_group-<grp>_netctrl_nodes_thr0p7.tsv
#   data/processed/K27M_group-<grp>_netctrl_driver_nodes_thr0p7.txt
#
# 산출물:
#   results/tables/K27M_module_summary_thr0p7.tsv
#   results/tables/K27M_module_group_summary_thr0p7.tsv
############################################################

## 0. 경로 설정 --------------------------------------------------------------

get_script_path <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  fileArg <- "--file="
  match <- grep(fileArg, cmdArgs)
  if (length(match) > 0) {
    return(normalizePath(sub(fileArg, "", cmdArgs[match])))
  }
  if (!is.null(sys.frames()[[1]]) &&
      "ofile" %in% names(sys.frames()[[1]])) {
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
  normalizePath(getwd())
}

script_path  <- get_script_path()
script_dir   <- dirname(script_path)
scripts_dir  <- dirname(script_dir)
project_root <- dirname(scripts_dir)

raw_dir            <- file.path(project_root, "data", "raw")
processed_dir      <- file.path(project_root, "data", "processed")
results_tables_dir <- file.path(project_root, "results", "tables")

cat("Project root       :", project_root,        "\n")
cat("Raw data dir       :", raw_dir,            "\n")
cat("Processed data dir :", processed_dir,      "\n")
cat("Results tables dir :", results_tables_dir, "\n\n")

threshold <- 0.7
thr_label <- gsub("\\.", "p", as.character(threshold))

if (!dir.exists(results_tables_dir)) {
  dir.create(results_tables_dir, recursive = TRUE, showWarnings = FALSE)
}

## 1. gene set 파일 읽기 ------------------------------------------------------

read_gene_set <- function(path) {
  if (!file.exists(path)) {
    stop("gene set 파일이 없습니다: ", path)
  }
  dat <- read.table(path, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE)
  if (ncol(dat) == 1) {
    genes <- dat[[1]]
  } else if ("gene" %in% names(dat)) {
    genes <- dat$gene
  } else if ("GeneSymbol" %in% names(dat)) {
    genes <- dat$GeneSymbol
  } else {
    stop("gene set 형식을 알 수 없습니다: ", path,
         " (한 컬럼 또는 gene / GeneSymbol 컬럼 필요)")
  }
  unique(toupper(genes))
}

cellcycle_path <- file.path(raw_dir, "cellcycle_genes.txt")
pdgfra_path    <- file.path(raw_dir, "PDGFRA_module_genes.txt")

cellcycle_genes <- read_gene_set(cellcycle_path)
pdgfra_genes    <- read_gene_set(pdgfra_path)

cat("cell-cycle gene 수 (리스트):", length(cellcycle_genes), "\n")
cat("PDGFRA module gene 수      :", length(pdgfra_genes),    "\n\n")

## 2. 전체 malignant 네트워크에서 모듈 summary -------------------------------

gene_sum_path <- file.path(
  results_tables_dir,
  paste0("K27M_netctrl_gene_summary_thr", thr_label, ".tsv")
)

if (!file.exists(gene_sum_path)) {
  stop("gene_summary 파일을 찾을 수 없습니다: ", gene_sum_path)
}

gene_sum <- read.table(gene_sum_path, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)

# gene 이름 대문자로 정규화
genes_upper <- toupper(gene_sum$gene)

gene_sum$is_cellcycle <- as.integer(genes_upper %in% cellcycle_genes)
gene_sum$is_pdgfra    <- as.integer(genes_upper %in% pdgfra_genes)

module_names <- c("cellcycle", "PDGFRA_module")
module_flags <- c("is_cellcycle", "is_pdgfra")

summaries_all <- list()

for (i in seq_along(module_names)) {
  mname <- module_names[i]
  flag  <- module_flags[i]

  in_set <- gene_sum[[flag]] == 1

  n_in_list  <- if (mname == "cellcycle") length(cellcycle_genes) else length(pdgfra_genes)
  n_in_net   <- sum(in_set)
  n_total    <- nrow(gene_sum)

  deg_set    <- gene_sum$deg_total[in_set]
  deg_other  <- gene_sum$deg_total[!in_set]

  drv_set    <- gene_sum$is_driver[in_set]      # 0/1
  drv_other  <- gene_sum$is_driver[!in_set]

  mean_deg_set   <- if (length(deg_set)   > 0) mean(deg_set)   else NA_real_
  mean_deg_other <- if (length(deg_other) > 0) mean(deg_other) else NA_real_

  frac_drv_set   <- if (length(drv_set)   > 0) mean(drv_set > 0)   else NA_real_
  frac_drv_other <- if (length(drv_other) > 0) mean(drv_other > 0) else NA_real_

  summaries_all[[i]] <- data.frame(
    level          = "All_malignant",
    module         = mname,
    threshold      = threshold,
    n_genes_total  = n_total,
    n_genes_in_set_list = n_in_list,
    n_genes_in_set_net  = n_in_net,
    mean_deg_in_set     = mean_deg_set,
    mean_deg_other      = mean_deg_other,
    frac_driver_in_set  = frac_drv_set,
    frac_driver_other   = frac_drv_other,
    stringsAsFactors = FALSE
  )
}

summary_all_df <- do.call(rbind, summaries_all)

## 3. 그룹별(OPC / AC / OC) 모듈 summary ------------------------------------

groups <- c("OPC", "AC", "OC")
group_summaries <- list()

for (g in groups) {

  cat("=== 그룹:", g, "===\n")

  # degree 파일
  deg_path <- file.path(
    processed_dir,
    paste0("K27M_group-", g, "_gene_degrees_thr", thr_label, ".tsv")
  )
  nodes_path <- file.path(
    processed_dir,
    paste0("K27M_group-", g, "_netctrl_nodes_thr", thr_label, ".tsv")
  )
  drivers_path <- file.path(
    processed_dir,
    paste0("K27M_group-", g, "_netctrl_driver_nodes_thr", thr_label, ".txt")
  )

  if (!file.exists(deg_path) ||
      !file.exists(nodes_path) ||
      !file.exists(drivers_path)) {
    warning("그룹 ", g, "의 입력 파일이 부족합니다. 건너뜀.")
    next
  }

  deg_df <- read.table(deg_path, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)
  nodes_df <- read.table(nodes_path, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)
  driver_ids <- scan(drivers_path, what = integer(), quiet = TRUE)

  # gene 이름 upper
  deg_df$gene_upper   <- toupper(deg_df$gene)
  nodes_df$gene_upper <- toupper(nodes_df$gene)

  # driver id -> gene
  driver_genes <- nodes_df$gene_upper[match(driver_ids, nodes_df$id)]

  # 그룹 네트워크 유전자 universe
  genes_g <- deg_df$gene_upper

  is_driver_g <- as.integer(genes_g %in% driver_genes)

  deg_df$is_driver   <- is_driver_g
  deg_df$is_cellcycle <- as.integer(genes_g %in% cellcycle_genes)
  deg_df$is_pdgfra    <- as.integer(genes_g %in% pdgfra_genes)

  n_total_g <- nrow(deg_df)

  for (i in seq_along(module_names)) {
    mname <- module_names[i]
    flag  <- module_flags[i]

    in_set <- deg_df[[flag]] == 1

    n_in_list  <- if (mname == "cellcycle") length(cellcycle_genes) else length(pdgfra_genes)
    n_in_net   <- sum(in_set)

    deg_set   <- deg_df$deg_total[in_set]
    deg_other <- deg_df$deg_total[!in_set]

    drv_set   <- deg_df$is_driver[in_set]
    drv_other <- deg_df$is_driver[!in_set]

    mean_deg_set   <- if (length(deg_set)   > 0) mean(deg_set)   else NA_real_
    mean_deg_other <- if (length(deg_other) > 0) mean(deg_other) else NA_real_

    frac_drv_set   <- if (length(drv_set)   > 0) mean(drv_set > 0)   else NA_real_
    frac_drv_other <- if (length(drv_other) > 0) mean(drv_other > 0) else NA_real_

    group_summaries[[length(group_summaries) + 1]] <- data.frame(
      level          = g,
      module         = mname,
      threshold      = threshold,
      n_genes_total  = n_total_g,
      n_genes_in_set_list = n_in_list,
      n_genes_in_set_net  = n_in_net,
      mean_deg_in_set     = mean_deg_set,
      mean_deg_other      = mean_deg_other,
      frac_driver_in_set  = frac_drv_set,
      frac_driver_other   = frac_drv_other,
      stringsAsFactors = FALSE
    )
  }
}

if (length(group_summaries) > 0) {
  summary_group_df <- do.call(rbind, group_summaries)
} else {
  summary_group_df <- NULL
}

## 4. 파일 저장 ---------------------------------------------------------------

module_summary_path <- file.path(
  results_tables_dir,
  paste0("K27M_module_summary_thr", thr_label, ".tsv")
)
write.table(
  summary_all_df,
  file      = module_summary_path,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

cat("전체 malignant module summary 저장:", module_summary_path, "\n")

if (!is.null(summary_group_df)) {
  module_group_summary_path <- file.path(
    results_tables_dir,
    paste0("K27M_module_group_summary_thr", thr_label, ".tsv")
  )
  write.table(
    summary_group_df,
    file      = module_group_summary_path,
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )
  cat("그룹별 module summary 저장:", module_group_summary_path, "\n")
}

cat("\n=== Aim3 module 분석 완료 ===\n")
############################################################
