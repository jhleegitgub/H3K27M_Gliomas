############################################################
# 02_build_gene_network.R
#
#  - 01_preprocess_HVG.R에서 만든 HVG 행렬(expr_hvg)을 읽어옴
#  - gene-gene 피어슨 상관계수 계산
#  - |r| >= threshold 인 쌍만 edge로 사용
#  - undirected edge + 임의 방향(directed) edge 리스트 저장
#
# 디렉터리 구조 가정:
#   프로젝트 루트/
#     data/processed/K27M_malignant_HVG2000_preprocessed.RData
#     results/tables/
############################################################

## 0. 경로 설정 함수 ---------------------------------------------------------

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

processed_dir     <- file.path(project_root, "data", "processed")
results_tables_dir <- file.path(project_root, "results", "tables")

if (!dir.exists(processed_dir)) {
  stop("data/processed 폴더가 없습니다: ", processed_dir)
}
if (!dir.exists(results_tables_dir)) {
  dir.create(results_tables_dir, recursive = TRUE, showWarnings = FALSE)
}

cat("Project root        :", project_root,        "\n")
cat("Processed data dir  :", processed_dir,      "\n")
cat("Results tables dir  :", results_tables_dir, "\n\n")

## 1. 전처리 RData 로드 -------------------------------------------------------

rdata_path <- file.path(processed_dir,
                        "K27M_malignant_HVG2000_preprocessed.RData")
if (!file.exists(rdata_path)) {
  stop("전처리 RData 파일을 찾을 수 없습니다: ", rdata_path)
}

load(rdata_path)
# 여기서 expr_hvg, expr_mal_log, expr_mal, meta_mal, tsne_mal 이 로드됨

if (!exists("expr_hvg")) {
  stop("expr_hvg 객체가 RData 안에 없습니다.")
}

cat("expr_hvg dim (genes x cells):",
    paste(dim(expr_hvg), collapse = " x "), "\n\n")

## 2. gene-gene 상관계수 계산 -------------------------------------------------

# expr_hvg: genes x cells -> t() 해서 cells x genes
cat("피어슨 상관계수 계산 중...\n")
cor_mat <- cor(t(expr_hvg), method = "pearson", use = "pairwise.complete.obs")

# 자기 자신 상관은 0으로
diag(cor_mat) <- 0

cat("cor_mat dim:", paste(dim(cor_mat), collapse = " x "), "\n\n")

## 3. threshold 적용해서 edge 리스트 만들기 -----------------------------------

# 필요하면 여기서 threshold 값만 바꿔서 여러 버전 만들 수 있음
threshold <- 0.7
cat("Threshold |r| >=", threshold, "로 edge 선택.\n")

# NA는 0 취급
cor_mat[is.na(cor_mat)] <- 0

# 상삼각 행렬에서 조건 만족하는 인덱스만 사용 (중복 방지)
sel <- which(abs(cor_mat) >= threshold & upper.tri(cor_mat), arr.ind = TRUE)

cat("선택된 edge 수 (undirected):", nrow(sel), "\n\n")

if (nrow(sel) == 0) {
  stop("선택된 edge가 0개입니다. threshold를 낮춰야 할 수 있습니다.")
}

genes <- rownames(cor_mat)

edges_undirected <- data.frame(
  gene1  = genes[sel[, 1]],
  gene2  = genes[sel[, 2]],
  weight = cor_mat[sel],
  stringsAsFactors = FALSE
)

## 4. netctrl용 방향 그래프(임의 방향) 만들기 --------------------------------

# 간단하게 알파벳 순으로 더 앞선 유전자를 source, 뒤를 target으로 둠
source <- ifelse(edges_undirected$gene1 < edges_undirected$gene2,
                 edges_undirected$gene1, edges_undirected$gene2)
target <- ifelse(edges_undirected$gene1 < edges_undirected$gene2,
                 edges_undirected$gene2, edges_undirected$gene1)

edges_directed <- data.frame(
  source = source,
  target = target,
  weight = edges_undirected$weight,
  stringsAsFactors = FALSE
)

## 5. 간단한 네트워크 요약 ----------------------------------------------------

n_nodes <- nrow(cor_mat)
n_edges <- nrow(edges_directed)

avg_degree <- 2 * n_edges / n_nodes
density    <- 2 * n_edges / (n_nodes * (n_nodes - 1))

summary_tbl <- data.frame(
  threshold   = threshold,
  n_nodes     = n_nodes,
  n_edges     = n_edges,
  avg_degree  = avg_degree,
  density     = density,
  stringsAsFactors = FALSE
)

cat("노드 수:", n_nodes, "\n")
cat("엣지 수:", n_edges, "\n")
cat("평균 degree:", round(avg_degree, 3), "\n")
cat("density:", signif(density, 3), "\n\n")

## 6. 파일로 저장 -------------------------------------------------------------

thr_label <- gsub("\\.", "p", as.character(threshold))

edges_undir_path <- file.path(
  processed_dir,
  paste0("K27M_gene_network_edges_undirected_thr", thr_label, ".tsv")
)
edges_dir_path <- file.path(
  processed_dir,
  paste0("K27M_gene_network_edges_directed_thr", thr_label, ".tsv")
)
summary_path <- file.path(
  results_tables_dir,
  paste0("K27M_gene_network_summary_thr", thr_label, ".tsv")
)

write.table(
  edges_undirected,
  file      = edges_undir_path,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

write.table(
  edges_directed,
  file      = edges_dir_path,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

write.table(
  summary_tbl,
  file      = summary_path,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

cat("=== 네트워크 생성 완료 ===\n")
cat("Undirected edges :", edges_undir_path, "\n")
cat("Directed edges   :", edges_dir_path,   "\n")
cat("Summary table    :", summary_path,     "\n")
############################################################
