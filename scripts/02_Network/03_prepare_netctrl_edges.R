############################################################
# 03_prepare_netctrl_edges.R
#
# - 02_build_gene_network.R 에서 만든 directed edge 파일을 읽음
# - netctrl 에서 쓰기 좋게
#     (1) 노드 리스트(id, gene)
#     (2) id 기반 edge 리스트 (u v)
#   를 생성해서 data/processed/ 에 저장
#
# 디렉터리 구조:
#   프로젝트 루트/
#     data/processed/
#       K27M_gene_network_edges_directed_thr0p7.tsv
#     scripts/02_Network/03_prepare_netctrl_edges.R
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

processed_dir <- file.path(project_root, "data", "processed")

if (!dir.exists(processed_dir)) {
  stop("data/processed 폴더가 없습니다: ", processed_dir)
}

cat("Project root   :", project_root,   "\n")
cat("Processed dir  :", processed_dir, "\n\n")

## 1. directed edge 파일 읽기 -------------------------------------------------

# 02_build_gene_network.R 에서 사용한 threshold 값과 라벨
threshold  <- 0.7
thr_label  <- gsub("\\.", "p", as.character(threshold))

edges_dir_path <- file.path(
  processed_dir,
  paste0("K27M_gene_network_edges_directed_thr", thr_label, ".tsv")
)

if (!file.exists(edges_dir_path)) {
  stop("directed edge 파일을 찾을 수 없습니다: ", edges_dir_path)
}

edges <- read.table(
  edges_dir_path,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Edge 파일 로드 완료. dim:",
    paste(dim(edges), collapse = " x "), "\n")
head(edges)

if (!all(c("source", "target") %in% colnames(edges))) {
  stop("edge 파일에 'source', 'target' 컬럼이 없습니다. 컬럼명을 확인하세요.")
}

## 2. 노드 리스트(id, gene) 생성 ---------------------------------------------

nodes <- sort(unique(c(edges$source, edges$target)))
n_nodes <- length(nodes)

node_df <- data.frame(
  id   = seq_len(n_nodes),
  gene = nodes,
  stringsAsFactors = FALSE
)

cat("노드 수:", n_nodes, "\n\n")

## 3. id 기반 edge 리스트(u, v) 생성 -----------------------------------------

# gene 이름을 id로 매핑하는 벡터
gene_to_id <- setNames(node_df$id, node_df$gene)

u <- gene_to_id[edges$source]
v <- gene_to_id[edges$target]

if (any(is.na(u)) || any(is.na(v))) {
  stop("일부 source/target gene 이 node 리스트에 없습니다. 매핑을 확인하세요.")
}

edges_id <- data.frame(
  u = as.integer(u),
  v = as.integer(v)
)

cat("엣지 수:", nrow(edges_id), "\n\n")

## 4. 간단 요약 ---------------------------------------------------------------

avg_outdeg <- nrow(edges_id) / n_nodes
cat("평균 out-degree (approx):", round(avg_outdeg, 3), "\n\n")

## 5. 파일 저장 ---------------------------------------------------------------

nodes_out_path <- file.path(
  processed_dir,
  paste0("K27M_netctrl_nodes_thr", thr_label, ".tsv")
)
edges_ids_out_path <- file.path(
  processed_dir,
  paste0("K27M_netctrl_edges_thr", thr_label, "_ids.txt")
)
edges_names_out_path <- file.path(
  processed_dir,
  paste0("K27M_netctrl_edges_thr", thr_label, "_withNames.tsv")
)

# (1) id, gene 매핑
write.table(
  node_df,
  file      = nodes_out_path,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

# (2) netctrl 입력용: u, v (헤더 없이 두 정수 컬럼)
write.table(
  edges_id,
  file      = edges_ids_out_path,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# (3) 사람이 보기 좋은 이름 버전 (optional)
edges_with_names <- cbind(
  edges[ , c("source", "target")],
  u = edges_id$u,
  v = edges_id$v,
  weight = if ("weight" %in% colnames(edges)) edges$weight else NA
)

write.table(
  edges_with_names,
  file      = edges_names_out_path,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

cat("=== netctrl 입력 파일 준비 완료 ===\n")
cat("노드 매핑   :", nodes_out_path,      "\n")
cat("엣지 (id)   :", edges_ids_out_path,  "\n")
cat("엣지 (이름) :", edges_names_out_path,"\n")
############################################################
