############################################################
# 05_build_group_networks.R
#
#  - malignant HVG 2000 데이터를 이용해
#    AC-like / OPC-like / OC-like 그룹별 co-expression 네트워크를 만들고
#    netctrl용 노드/엣지 파일을 생성한다.
#
# 사용 파일:
#   data/processed/K27M_malignant_HVG2000_preprocessed.RData
#     - expr_hvg (genes x cells, malignant only)
#     - meta_mal (cells x metadata, Cellcycle / AC-like / OPC-like / OC-like 등)
#
# 산출물 (data/processed):
#   K27M_group-<grp>_netctrl_nodes_thr0p7.tsv
#   K27M_group-<grp>_netctrl_edges_thr0p7_ids.txt
#   K27M_group-<grp>_gene_degrees_thr0p7.tsv
#   K27M_group-<grp>_gene_network_edges_undirected_thr0p7.tsv
#
# 산출물 (results/tables):
#   K27M_gene_network_group_summary_thr0p7.tsv
#
# 여기서 <grp> ∈ {OPC, AC, OC}
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

processed_dir      <- file.path(project_root, "data", "processed")
results_tables_dir <- file.path(project_root, "results", "tables")

if (!dir.exists(processed_dir)) {
  stop("data/processed 폴더가 없습니다: ", processed_dir)
}
if (!dir.exists(results_tables_dir)) {
  dir.create(results_tables_dir, recursive = TRUE, showWarnings = FALSE)
}

cat("Project root       :", project_root,        "\n")
cat("Processed data dir :", processed_dir,      "\n")
cat("Results tables dir :", results_tables_dir, "\n\n")

threshold <- 0.7
thr_label <- gsub("\\.", "p", as.character(threshold))

## 1. 전처리된 malignant RData 로드 ------------------------------------------

rdata_path <- file.path(
  processed_dir,
  "K27M_malignant_HVG2000_preprocessed.RData"
)

if (!file.exists(rdata_path)) {
  stop("전처리 RData 파일이 없습니다: ", rdata_path)
}

load(rdata_path)
# 필요 객체: expr_hvg (genes x cells), meta_mal (cells x metadata)

if (!exists("expr_hvg") || !exists("meta_mal")) {
  stop("expr_hvg 혹은 meta_mal 객체를 찾을 수 없습니다.")
}

cat("expr_hvg dim (genes x cells):",
    paste(dim(expr_hvg), collapse = " x "), "\n")
cat("meta_mal dim (cells x vars):",
    paste(dim(meta_mal), collapse = " x "), "\n\n")

## 2. AC-like / OPC-like / OC-like 그룹 라벨링 ------------------------------

# Portal 메타데이터의 프로그램 점수 컬럼 이름:
prog_cols <- c("OPC-like", "AC-like", "OC-like")
missing <- setdiff(prog_cols, colnames(meta_mal))
if (length(missing) > 0) {
  stop("meta_mal 에서 다음 컬럼을 찾을 수 없습니다: ",
       paste(missing, collapse = ", "))
}

prog_mat <- as.matrix(meta_mal[, prog_cols])

# 각 셀에 대해 가장 높은 프로그램을 그룹으로 지정
grp_idx <- max.col(prog_mat, ties.method = "first")
group3  <- c("OPC", "AC", "OC")[grp_idx]

meta_mal$Group3 <- group3

table_group <- table(meta_mal$Group3)
cat("Group3 분포:\n")
print(table_group)
cat("\n")

## 3. 그룹별 네트워크 생성 함수 ----------------------------------------------

build_group_network <- function(group_label) {

  cat("=== 그룹:", group_label, "===\n")

  # 3-1. 해당 그룹 셀 "번호" 선택 (meta_mal 행 인덱스) ----
  idx_in_group <- which(meta_mal$Group3 == group_label)

  if (length(idx_in_group) < 10) {
    warning("그룹 ", group_label, " 에 셀 수가 너무 적습니다: ",
            length(idx_in_group))
  }

  # expr_hvg 열 순서와 meta_mal 행 순서가 일치한다고 가정하고
  # 인덱스로 컬럼 선택
  expr_g <- expr_hvg[, idx_in_group, drop = FALSE]

  cat("  셀 수:", ncol(expr_g), "\n")

  # 3-2. 변동이 없는 유전자 제거 ------------------------
  vars <- apply(expr_g, 1, var)
  keep_genes <- vars > 0
  expr_g <- expr_g[keep_genes, , drop = FALSE]

  genes <- rownames(expr_g)
  n_genes <- length(genes)

  cat("  유전자 수 (var>0):", n_genes, "\n")

  # 3-3. 상관계수 계산 -------------------------------
  cat("  상관계수 계산 중...\n")
  cor_mat <- cor(t(expr_g), method = "pearson")
  diag(cor_mat) <- 0

  # 3-4. threshold 적용하여 undirected edge 선정 ------
  sel <- which(abs(cor_mat) >= threshold, arr.ind = TRUE)
  sel <- sel[sel[, "row"] < sel[, "col"], , drop = FALSE]

  if (nrow(sel) == 0) {
    warning("그룹 ", group_label, " 에서 threshold 를 넘는 edge 가 없습니다.")
  }

  edges_undir <- data.frame(
    source = genes[sel[, "row"]],
    target = genes[sel[, "col"]],
    weight = cor_mat[sel],
    stringsAsFactors = FALSE
  )

  n_edges_undir <- nrow(edges_undir)

  cat("  선택된 undirected edge 수:", n_edges_undir, "\n")

  # 3-5. degree 계산 --------------------------------
  deg_tab <- table(c(edges_undir$source, edges_undir$target))
  deg_tot <- rep(0L, n_genes)
  names(deg_tot) <- genes
  deg_tot[names(deg_tab)] <- as.integer(deg_tab)

  degrees_df <- data.frame(
    gene      = genes,
    deg_total = as.integer(deg_tot),
    group     = group_label,
    stringsAsFactors = FALSE
  )

  # 3-6. netctrl 노드/엣지 파일 생성 -----------------
  nodes_df <- data.frame(
    id   = seq_len(n_genes),
    gene = genes,
    stringsAsFactors = FALSE
  )

  u <- match(edges_undir$source, genes)
  v <- match(edges_undir$target, genes)

  # 양방향 edge 추가 (directed)
  edges_dir <- rbind(
    data.frame(u = u, v = v, stringsAsFactors = FALSE),
    data.frame(u = v, v = u, stringsAsFactors = FALSE)
  )

  # 파일 이름 prefix
  prefix <- paste0("K27M_group-", group_label)

  nodes_path <- file.path(
    processed_dir,
    paste0(prefix, "_netctrl_nodes_thr", thr_label, ".tsv")
  )
  edges_ids_path <- file.path(
    processed_dir,
    paste0(prefix, "_netctrl_edges_thr", thr_label, "_ids.txt")
  )
  degrees_path <- file.path(
    processed_dir,
    paste0(prefix, "_gene_degrees_thr", thr_label, ".tsv")
  )
  edges_undir_path <- file.path(
    processed_dir,
    paste0(prefix, "_gene_network_edges_undirected_thr", thr_label, ".tsv")
  )

  write.table(
    nodes_df,
    file      = nodes_path,
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )
  write.table(
    edges_dir,
    file      = edges_ids_path,
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  write.table(
    degrees_df,
    file      = degrees_path,
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )
  write.table(
    edges_undir,
    file      = edges_undir_path,
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )

  cat("  파일 저장 완료:\n")
  cat("    nodes   :", nodes_path,      "\n")
  cat("    edges   :", edges_ids_path,  "\n")
  cat("    degrees :", degrees_path,    "\n")
  cat("    undir   :", edges_undir_path,"\n\n")

  # 3-7. 네트워크 요약 값 리턴 ----------------------
  avg_deg <- mean(deg_tot)
  density <- if (n_genes > 1) {
    n_edges_undir / (n_genes * (n_genes - 1) / 2)
  } else {
    NA_real_
  }

  data.frame(
    group      = group_label,
    threshold  = threshold,
    n_nodes    = n_genes,
    n_edges    = n_edges_undir,
    avg_degree = avg_deg,
    density    = density,
    stringsAsFactors = FALSE
  )
}

## 4. 그룹별 네트워크 생성 루프 ----------------------------------------------

groups <- c("OPC", "AC", "OC")
summary_list <- vector("list", length(groups))

for (i in seq_along(groups)) {
  summary_list[[i]] <- build_group_network(groups[[i]])
}

summary_df <- do.call(rbind, summary_list)

summary_path <- file.path(
  results_tables_dir,
  paste0("K27M_gene_network_group_summary_thr", thr_label, ".tsv")
)

write.table(
  summary_df,
  file      = summary_path,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

cat("=== 그룹별 네트워크 생성 완료 ===\n")
cat("요약 테이블:", summary_path, "\n")
############################################################
