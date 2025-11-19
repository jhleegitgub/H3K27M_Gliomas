############################################################
# 05_analyze_netctrl_results.R
#
#  - 04_compute_controllability_python.py 가 만든 결과를 읽어서
#    * n_D (driver node fraction)
#    * edge type 비율
#    * 각 gene의 degree / driver 여부 / critical-source 여부
#    를 요약 테이블로 저장하는 스크립트
#
# 사용 파일 (data/processed/):
#   K27M_netctrl_nodes_thr0p7.tsv           (id, gene)
#   K27M_netctrl_edges_thr0p7_ids.txt       (u, v)
#   K27M_netctrl_driver_nodes_thr0p7.txt    (id)
#   K27M_netctrl_critical_edges_thr0p7.txt  (u, v)
#   K27M_netctrl_redundant_edges_thr0p7.txt (u, v)
#   K27M_netctrl_ordinary_edges_thr0p7.txt  (u, v)
#
# 산출물 (results/tables/):
#   K27M_netctrl_summary_thr0p7.tsv
#   K27M_netctrl_gene_summary_thr0p7.tsv
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

## 1. 노드/엣지 기본 정보 로드 -----------------------------------------------

nodes_path <- file.path(
  processed_dir,
  paste0("K27M_netctrl_nodes_thr", thr_label, ".tsv")
)
edges_ids_path <- file.path(
  processed_dir,
  paste0("K27M_netctrl_edges_thr", thr_label, "_ids.txt")
)

if (!file.exists(nodes_path)) stop("노드 파일이 없습니다: ", nodes_path)
if (!file.exists(edges_ids_path)) stop("엣지 파일이 없습니다: ", edges_ids_path)

nodes <- read.table(
  nodes_path,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)
edges_all <- read.table(
  edges_ids_path,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)
colnames(edges_all) <- c("u", "v")

n_nodes <- nrow(nodes)
n_edges <- nrow(edges_all)

cat("노드 수:", n_nodes, "\n")
cat("엣지 수:", n_edges, "\n\n")

## 2. controllability 결과 파일 로드 -----------------------------------------

driver_path    <- file.path(processed_dir,
                            paste0("K27M_netctrl_driver_nodes_thr", thr_label, ".txt"))
critical_path  <- file.path(processed_dir,
                            paste0("K27M_netctrl_critical_edges_thr", thr_label, ".txt"))
redundant_path <- file.path(processed_dir,
                            paste0("K27M_netctrl_redundant_edges_thr", thr_label, ".txt"))
ordinary_path  <- file.path(processed_dir,
                            paste0("K27M_netctrl_ordinary_edges_thr", thr_label, ".txt"))

for (p in c(driver_path, critical_path, redundant_path, ordinary_path)) {
  if (!file.exists(p)) {
    stop("controllability 결과 파일이 없습니다 (이름/위치 확인): ", p)
  }
}

driver_ids <- scan(driver_path, what = integer(), quiet = TRUE)

critical_edges  <- read.table(critical_path,  header = FALSE, sep = "\t",
                              stringsAsFactors = FALSE)
redundant_edges <- read.table(redundant_path, header = FALSE, sep = "\t",
                              stringsAsFactors = FALSE)
ordinary_edges  <- read.table(ordinary_path,  header = FALSE, sep = "\t",
                              stringsAsFactors = FALSE)

colnames(critical_edges)  <- c("u", "v")
colnames(redundant_edges) <- c("u", "v")
colnames(ordinary_edges)  <- c("u", "v")

cat("driver 노드 수 :", length(driver_ids), "\n")
cat("critical edge 수 :",  nrow(critical_edges),  "\n")
cat("redundant edge 수:", nrow(redundant_edges), "\n")
cat("ordinary edge 수 :", nrow(ordinary_edges),  "\n\n")

## 3. 요약 통계 (n_D, edge type 비율 등) -------------------------------------

nD <- length(driver_ids) / n_nodes

edge_counts <- c(
  critical  = nrow(critical_edges),
  redundant = nrow(redundant_edges),
  ordinary  = nrow(ordinary_edges)
)
edge_frac <- edge_counts / sum(edge_counts)

summary_netctrl <- data.frame(
  threshold     = threshold,
  n_nodes       = n_nodes,
  n_edges       = n_edges,
  n_driver      = length(driver_ids),
  nD            = nD,
  n_critical    = edge_counts["critical"],
  n_redundant   = edge_counts["redundant"],
  n_ordinary    = edge_counts["ordinary"],
  frac_critical  = edge_frac["critical"],
  frac_redundant = edge_frac["redundant"],
  frac_ordinary  = edge_frac["ordinary"],
  stringsAsFactors = FALSE
)

print(summary_netctrl)
cat("\n")

## 4. gene 단위 요약 (degree, driver, critical-source) -----------------------

# out/in/total degree
deg_out <- tabulate(edges_all$u, nbins = n_nodes)
deg_in  <- tabulate(edges_all$v, nbins = n_nodes)
deg_tot <- deg_out + deg_in

# driver 플래그
is_driver <- integer(n_nodes)
is_driver[driver_ids] <- 1

# critical edge의 source id
crit_source_ids <- unique(critical_edges$u)
is_crit_source  <- integer(n_nodes)
is_crit_source[crit_source_ids] <- 1

gene_summary <- data.frame(
  id             = nodes$id,
  gene           = nodes$gene,
  deg_out        = deg_out,
  deg_in         = deg_in,
  deg_total      = deg_tot,
  is_driver      = is_driver,
  is_crit_source = is_crit_source,
  stringsAsFactors = FALSE
)

## 5. 파일 저장 ---------------------------------------------------------------

summary_path <- file.path(
  results_tables_dir,
  paste0("K27M_netctrl_summary_thr", thr_label, ".tsv")
)
gene_summary_path <- file.path(
  results_tables_dir,
  paste0("K27M_netctrl_gene_summary_thr", thr_label, ".tsv")
)

write.table(
  summary_netctrl,
  file      = summary_path,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

write.table(
  gene_summary,
  file      = gene_summary_path,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

cat("=== controllability 결과 분석 완료 ===\n")
cat("요약 테이블:", summary_path,      "\n")
cat("유전자 요약:", gene_summary_path, "\n")
############################################################
