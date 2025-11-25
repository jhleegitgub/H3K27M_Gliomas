############################################################
# 07_TF_annotation_and_summary.R
#
#  - human TF 리스트를 gene_summary 에 붙여서
#    * is_TF
#    * driver TF / non-driver TF
#    * critical-source TF
#    * high-degree TF (상위 x%)
#    를 요약하는 스크립트
#
# 사용 파일:
#   data/raw/human_TF_list.txt      (한 줄에 하나 TF 심볼)
#   results/tables/K27M_netctrl_gene_summary_thr0p7.tsv
#
# 산출물:
#   results/tables/K27M_TF_summary_thr0p7.tsv
#   results/tables/K27M_TF_lists_thr0p7.tsv   (각 카테고리별 gene 리스트)
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
results_tables_dir <- file.path(project_root, "results", "tables")

if (!dir.exists(raw_dir)) {
  stop("data/raw 폴더가 없습니다: ", raw_dir)
}

cat("Project root   :", project_root,        "\n")
cat("Raw data dir   :", raw_dir,            "\n")
cat("Results tables :", results_tables_dir, "\n\n")

threshold <- 0.7
thr_label <- gsub("\\.", "p", as.character(threshold))

## 1. gene_summary 로드 -------------------------------------------------------

gene_summary_path <- file.path(
  results_tables_dir,
  paste0("K27M_netctrl_gene_summary_thr", thr_label, ".tsv")
)
if (!file.exists(gene_summary_path)) {
  stop("gene_summary 파일을 찾을 수 없습니다: ", gene_summary_path)
}

gene_summary <- read.table(
  gene_summary_path,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

cat("gene_summary dim:", paste(dim(gene_summary), collapse = " x "), "\n\n")

## 2. human TF 리스트 로드 ----------------------------------------------------

tf_path <- file.path(raw_dir, "human_TF_list.txt")
if (!file.exists(tf_path)) {
  stop("TF 리스트 파일(human_TF_list.txt)을 data/raw 에 넣어주세요: ", tf_path)
}

# 헤더가 있을 수도, 없을 수도 있으니 둘 다 처리
tf_raw <- read.table(tf_path, header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE)
if (ncol(tf_raw) == 1) {
  # 한 컬럼이면 그게 바로 심볼이라고 가정
  tf_symbols <- unique(tf_raw[[1]])
} else if ("GeneSymbol" %in% colnames(tf_raw)) {
  tf_symbols <- unique(tf_raw$GeneSymbol)
} else {
  stop("TF 파일 형식을 파악할 수 없습니다. 한 컬럼 또는 'GeneSymbol' 컬럼이 필요.")
}

tf_symbols <- toupper(tf_symbols)  # 대문자로 통일

cat("TF 개수 (리스트):", length(tf_symbols), "\n\n")

## 3. gene_summary에 TF 여부 붙이기 -----------------------------------------

genes_upper <- toupper(gene_summary$gene)
is_TF <- genes_upper %in% tf_symbols

gene_summary$is_TF <- as.integer(is_TF)

cat("네트워크 노드 중 TF 개수:",
    sum(gene_summary$is_TF), "/", nrow(gene_summary), "\n\n")

## 4. high-degree TF, driver TF, critical-source TF 정의 ----------------------

# high-degree TF: degree 상위 10% 안에 들어가는 TF
deg <- gene_summary$deg_total
cutoff_q <- 0.9
deg_cutoff <- quantile(deg, probs = cutoff_q)

gene_summary$is_highdeg <- as.integer(deg >= deg_cutoff)

# 이미 is_driver, is_crit_source 컬럼 있음
# (05_analyze_netctrl_results.R 에서 생성)

## 5. 카테고리별 개수 요약 ----------------------------------------------------

total_TF         <- sum(gene_summary$is_TF)
driver_TF        <- sum(gene_summary$is_TF == 1 & gene_summary$is_driver == 1)
crit_source_TF   <- sum(gene_summary$is_TF == 1 & gene_summary$is_crit_source == 1)
highdeg_TF       <- sum(gene_summary$is_TF == 1 & gene_summary$is_highdeg == 1)

summary_TF <- data.frame(
  threshold       = threshold,
  n_genes         = nrow(gene_summary),
  n_TF_in_list    = length(tf_symbols),
  n_TF_in_network = total_TF,
  driver_TF       = driver_TF,
  crit_source_TF  = crit_source_TF,
  highdeg_TF      = highdeg_TF,
  frac_TF_in_net  = total_TF / nrow(gene_summary),
  frac_driver_TF  = ifelse(total_TF > 0, driver_TF / total_TF, NA),
  frac_crit_TF    = ifelse(total_TF > 0, crit_source_TF / total_TF, NA),
  frac_highdeg_TF = ifelse(total_TF > 0, highdeg_TF / total_TF, NA),
  deg_cutoff_top10pct = as.numeric(deg_cutoff),
  stringsAsFactors = FALSE
)

print(summary_TF)
cat("\n")

## 6. 카테고리별 TF 리스트 테이블 생성 ---------------------------------------

gene_summary$category <- "other"

gene_summary$category[gene_summary$is_TF == 1] <- "TF_only"
gene_summary$category[gene_summary$is_TF == 1 & gene_summary$is_driver == 1] <- "driver_TF"
gene_summary$category[gene_summary$is_TF == 1 & gene_summary$is_crit_source == 1] <- "critical_source_TF"
gene_summary$category[gene_summary$is_TF == 1 & gene_summary$is_highdeg == 1] <- "highdeg_TF"

# (중복 카테고리 가능하니, 엄밀히 disjoint 필요하면 나중에 조건 더 세밀히 나눠도 됨)

TF_lists <- subset(gene_summary, is_TF == 1,
                   select = c(gene, deg_total, is_driver, is_crit_source,
                              is_highdeg, category))

## 7. 파일 저장 ---------------------------------------------------------------

TF_summary_path <- file.path(
  results_tables_dir,
  paste0("K27M_TF_summary_thr", thr_label, ".tsv")
)
TF_lists_path <- file.path(
  results_tables_dir,
  paste0("K27M_TF_lists_thr", thr_label, ".tsv")
)

write.table(
  summary_TF,
  file      = TF_summary_path,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

write.table(
  TF_lists,
  file      = TF_lists_path,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

cat("=== TF annotation & summary 완료 ===\n")
cat("TF summary :", TF_summary_path, "\n")
cat("TF lists  :", TF_lists_path,   "\n")
############################################################
