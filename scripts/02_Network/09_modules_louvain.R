#!/usr/bin/env Rscript

## ------------------------------------------------------------
## 03_modules_louvain.R
## - malignant 전체 co-expression network에서 Louvain 모듈 찾기
## - netctrl gene summary와 합쳐서 모듈별 요약 만들기
## - 어디 디렉토리에서 실행해도 동작하게 경로 자동 설정
## ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
})

## ------------------------------------------------------------
## 0. 스크립트 위치 기준으로 레포 루트 경로 찾기
## ------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])

if (length(script_path) == 0) {
  # RStudio/source() 등일 때는 현재 작업 디렉토리 사용
  script_dir <- getwd()
} else {
  script_dir <- dirname(normalizePath(script_path))
}

# script_dir = .../H3K27M_Gliomas/scripts/02_Network
# → 루트 = 그 기준 두 단계 위
root_dir <- normalizePath(file.path(script_dir, "..", ".."))

cat("script_dir :", script_dir,  "\n")
cat("root_dir   :", root_dir,    "\n")

## 파일 경로 (루트 기준)
edge_path <- file.path(
  root_dir,
  "data", "processed",
  "K27M_gene_network_edges_undirected_thr0p7.tsv"
)

gene_sum_path <- file.path(
  root_dir,
  "results", "tables",
  "K27M_netctrl_gene_summary_thr0p7.tsv"
)

cat("edge_path  :", edge_path,    "\n")
cat("gene_sum   :", gene_sum_path, "\n\n")

## ------------------------------------------------------------
## 1. Edge 리스트 읽어서 igraph 객체 만들기
## ------------------------------------------------------------

edges <- fread(edge_path)

# edge 파일에서 어떤 컬럼이 gene1/gene2 역할인지 자동 추론
possible_pairs <- list(
  c("gene1", "gene2"),
  c("source", "target"),
  c("from", "to")
)

pair_found <- FALSE
for (p in possible_pairs) {
  if (all(p %in% names(edges))) {
    setnames(edges, p, c("source", "target"), skip_absent = TRUE)
    pair_found <- TRUE
    break
  }
}

if (!pair_found) {
  stop("edge 파일에서 gene 쌍 컬럼을 찾지 못했습니다. (gene1/gene2 또는 source/target 등 확인 필요)")
}

g <- graph_from_data_frame(
  d = edges[, .(source, target)],
  directed = FALSE
)

g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

cat("노드 수 :", gorder(g), " / 엣지 수 :", gsize(g), "\n")

## ------------------------------------------------------------
## 2. Louvain community detection
## ------------------------------------------------------------

set.seed(42)
lc <- cluster_louvain(g)

cat("모듈 개수:", length(lc), "\n")
cat("상위 10개 모듈 사이즈:\n")
print(sort(sizes(lc), decreasing = TRUE)[1:min(10, length(lc))])

membership_vec <- membership(lc)

modules_dt <- data.table(
  gene   = names(membership_vec),
  module = as.integer(membership_vec)
)

## ------------------------------------------------------------
## 3. netctrl gene summary와 merge
## ------------------------------------------------------------

gene_sum <- fread(gene_sum_path)

if (!"gene" %in% names(gene_sum)) {
  stop("gene summary 파일에 'gene' 컬럼이 없습니다. 컬럼명을 확인해주세요.")
}

annot <- merge(
  gene_sum,
  modules_dt,
  by = "gene",
  all.y = TRUE   # 모듈에 포함된 gene은 모두 유지
)

cat("annot 테이블 차원:", dim(annot)[1], "x", dim(annot)[2], "\n")

## ------------------------------------------------------------
## 4. TF / cell-cycle / PDGFRA / driver 컬럼 자동 인식 + 0/1 변환
##    (실제 컬럼 이름에 맞춰서 후보만 추가해주면 됨)
## ------------------------------------------------------------

# 후보 이름들 정의
candidate_TF        <- c("is_TF", "is_tf", "TF", "isTF")
candidate_cellcycle <- c("is_cellcycle", "is_cell_cycle", "cell_cycle", "isCellCycle")
candidate_PDGFRA    <- c("in_PDGFRA_module", "in_PDGFRA", "PDGFRA_module")
candidate_driver    <- c("is_driver", "driver", "isDriver")

get_col <- function(cands, nms) {
  x <- intersect(cands, nms)
  if (length(x) == 0) NA_character_ else x[1]
}

tf_col        <- get_col(candidate_TF,        names(annot))
cellcycle_col <- get_col(candidate_cellcycle, names(annot))
pdgfra_col    <- get_col(candidate_PDGFRA,    names(annot))
driver_col    <- get_col(candidate_driver,    names(annot))

logical_cols <- c(tf_col, cellcycle_col, pdgfra_col, driver_col)
logical_cols <- logical_cols[!is.na(logical_cols)]

for (col in logical_cols) {
  annot[[col]] <- as.integer(as.logical(annot[[col]]))
}

## ------------------------------------------------------------
## 5. 모듈 단위 summary 테이블
## ------------------------------------------------------------

module_summary <- annot[, .(
  n_genes = .N,
  n_TF        = if (!is.na(tf_col))        sum(get(tf_col),        na.rm = TRUE) else NA_integer_,
  n_cellcycle = if (!is.na(cellcycle_col)) sum(get(cellcycle_col), na.rm = TRUE) else NA_integer_,
  n_PDGFRA    = if (!is.na(pdgfra_col))    sum(get(pdgfra_col),    na.rm = TRUE) else NA_integer_,
  n_driver    = if (!is.na(driver_col))    sum(get(driver_col),    na.rm = TRUE) else NA_integer_
), by = module][order(-n_genes)]

cat("모듈 summary 상위 10개:\n")
print(head(module_summary, 10))

## ------------------------------------------------------------
## 6. 파일 저장 (루트 기준 results/tables)
## ------------------------------------------------------------

out_gene_annot_path <- file.path(
  root_dir,
  "results", "tables",
  "K27M_modules_louvain_gene_table.tsv"
)

out_module_summary_path <- file.path(
  root_dir,
  "results", "tables",
  "K27M_modules_louvain_summary.tsv"
)

fwrite(annot,         out_gene_annot_path,    sep = "\t")
fwrite(module_summary, out_module_summary_path, sep = "\t")

cat("gene-level 모듈 annot 저장:", out_gene_annot_path, "\n")
cat("module summary 저장      :", out_module_summary_path, "\n")

## ------------------------------------------------------------
## 7. 상위 모듈 몇 개 간단히 프린트
## ------------------------------------------------------------

top_k <- 5
top_modules <- module_summary[1:min(top_k, .N), module]

cat("상위", top_k, "개 모듈 번호:", paste(top_modules, collapse = " "), "\n")

for (m in top_modules) {
  cat("\n--- module", m, "---\n")
  sub <- annot[module == m]
  cat("gene 수 :", nrow(sub), "\n")
  
  if (!is.na(driver_col)) {
    cat("  driver 수:", sum(sub[[driver_col]], na.rm = TRUE), "\n")
  }
}

## (옵션) 인터랙티브 세션에서만 히스토그램
if (interactive()) {
  hist(
    module_summary$n_genes,
    breaks = 30,
    main = "Module size distribution (Louvain, malignant network)",
    xlab  = "Number of genes in module"
  )
}

cat("\n완료.\n")
