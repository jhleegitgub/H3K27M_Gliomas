############################################################
# 01_preprocess_HVG.R
#
# H3K27M_Gliomas 프로젝트용 전처리 스크립트
# 레포 구조 가정:
#   프로젝트 루트/
#     data/
#       raw/        (입력 파일)
#       processed/  (전처리 산출물)
#     scripts/
#       01_Preprocess/01_preprocess_HVG.R
#
# 1) Metadata / tSNE / RSEM 읽기
# 2) 공통 샘플 맞추기
# 3) Malignant cell만 선택
# 4) log2 변환 + HVG 상위 2000개 추출
# 5) data/processed/ 에 결과 저장
############################################################

## 0. 스크립트/프로젝트 경로 설정 ---------------------------------------------

get_script_path <- function() {
  # Rscript로 실행할 때
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  fileArg <- "--file="
  match <- grep(fileArg, cmdArgs)
  if (length(match) > 0) {
    return(normalizePath(sub(fileArg, "", cmdArgs[match])))
  }
  # RStudio 등에서 source() 할 때
  if (!is.null(sys.frames()[[1]]) &&
      "ofile" %in% names(sys.frames()[[1]])) {
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
  # 실패하면 현재 작업 디렉터리 사용
  normalizePath(getwd())
}

script_path  <- get_script_path()
script_dir   <- dirname(script_path)
scripts_dir  <- dirname(script_dir)
project_root <- dirname(scripts_dir)

raw_dir       <- file.path(project_root, "data", "raw")
processed_dir <- file.path(project_root, "data", "processed")

if (!dir.exists(raw_dir)) {
  stop("data/raw 폴더를 찾을 수 없습니다: ", raw_dir)
}
if (!dir.exists(processed_dir)) {
  dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
}

cat("Project root : ", project_root,  "\n")
cat("Raw data dir : ", raw_dir,       "\n")
cat("Proc data dir: ", processed_dir, "\n\n")

## 1. 패키지 로드 -------------------------------------------------------------

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(dplyr)

## 2. Metadata 읽기 -----------------------------------------------------------

meta_path <- file.path(raw_dir, "PortalK27M_Metadata.vh20180223.txt")
if (!file.exists(meta_path)) {
  stop("Metadata 파일을 찾을 수 없습니다: ", meta_path)
}

meta <- read.table(
  meta_path,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 첫 데이터 행이 TYPE 설명행이면 제거
if (nrow(meta) > 0 && grepl("TYPE", meta[1, 1], ignore.case = TRUE)) {
  meta <- meta[-1, ]
}

# 숫자 칼럼 numeric 변환 (존재하는 것만)
num_cols <- c(
  "GenesExpressed",
  "HousekeepingGeneExpression",
  "Cellcycle",
  "OPC-variable",
  "OC-like",
  "AC-like",
  "OPC-like"
)
for (cn in num_cols) {
  if (cn %in% colnames(meta)) {
    meta[[cn]] <- as.numeric(meta[[cn]])
  }
}

cat("Metadata dim:", paste(dim(meta), collapse = " x "), "\n")
cat("Metadata Type table:\n")
print(table(meta$Type))
cat("\n")

## 3. t-SNE 좌표 읽기 ---------------------------------------------------------

tsne_path <- file.path(raw_dir, "PortalK27M_tSNE.vh20180223.txt")
if (!file.exists(tsne_path)) {
  stop("tSNE 파일을 찾을 수 없습니다: ", tsne_path)
}

tsne <- read.table(
  tsne_path,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 첫 행이 TYPE 설명행이면 제거
if (nrow(tsne) > 0 && (grepl("TYPE", tsne[1, 1], ignore.case = TRUE) ||
                       tsne[1, "NAME"] == "TYPE")) {
  tsne <- tsne[-1, ]
}

if ("X" %in% colnames(tsne)) tsne$X <- as.numeric(tsne$X)
if ("Y" %in% colnames(tsne)) tsne$Y <- as.numeric(tsne$Y)

cat("tSNE dim:", paste(dim(tsne), collapse = " x "), "\n\n")

## 4. RSEM 발현 데이터 읽기 ---------------------------------------------------

expr_path <- file.path(raw_dir, "K27Mproject.RSEM.vh20170621.txt")
if (!file.exists(expr_path)) {
  stop("Expression 파일을 찾을 수 없습니다: ", expr_path)
}

expr <- read.table(
  expr_path,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 첫 열(유전자 이름)을 rownames로 옮기고 제거
rownames(expr) <- expr[[1]]
expr <- expr[, -1]

cat("Expression dim (genes x cells):",
    paste(dim(expr), collapse = " x "), "\n\n")

## 5. 공통 샘플 맞추기 --------------------------------------------------------

if (!("NAME" %in% colnames(meta))) {
  stop("Metadata에 NAME 칼럼이 없습니다. 실제 샘플 ID 칼럼 이름을 확인하세요.")
}
if (!("NAME" %in% colnames(tsne))) {
  stop("tSNE에 NAME 칼럼이 없습니다. 실제 샘플 ID 칼럼 이름을 확인하세요.")
}

samples_meta <- meta$NAME
samples_tsne <- tsne$NAME
samples_expr <- colnames(expr)

common <- Reduce(intersect, list(samples_expr, samples_meta, samples_tsne))
cat("공통 샘플 수:", length(common), "\n\n")

if (length(common) == 0) {
  stop("세 파일 간 공통 샘플이 0개입니다. 샘플 이름을 다시 확인하세요.")
}

expr <- expr[, common]
meta <- meta[match(common, meta$NAME), ]
tsne <- tsne[match(common, tsne$NAME), ]

stopifnot(all(colnames(expr) == meta$NAME))
stopifnot(all(meta$NAME == tsne$NAME))

cat("정렬 완료: expr / meta / tsne 샘플 순서 일치.\n\n")

## 6. Malignant cell만 선택 ---------------------------------------------------

if (!("Type" %in% colnames(meta))) {
  stop("Metadata에 Type 칼럼이 없습니다. malignant 정보를 가진 칼럼 이름을 확인하세요.")
}

cat("Type 분포:\n")
print(table(meta$Type))
cat("\n")

malignant_idx <- meta$Type == "Malignant"
if (sum(malignant_idx) == 0) {
  stop("Malignant cell이 0개입니다. Type 값이 정확한지 확인하세요.")
}

expr_mal <- expr[, malignant_idx]
meta_mal <- meta[malignant_idx, ]
tsne_mal <- tsne[malignant_idx, ]

cat("Malignant expression dim (genes x cells):",
    paste(dim(expr_mal), collapse = " x "), "\n")
cat("Malignant cell 수:", ncol(expr_mal), "\n\n")

## 7. log2 변환 ---------------------------------------------------------------

expr_mal_log <- log2(as.matrix(expr_mal) + 1)

## 8. HVG 상위 2000개 선택 ----------------------------------------------------

gene_mean <- rowMeans(expr_mal_log)
gene_sd   <- apply(expr_mal_log, 1, sd)
cv        <- gene_sd / gene_mean

keep    <- is.finite(cv) & gene_mean > 0
cv_keep <- cv[keep]

top_n <- 2000
if (length(cv_keep) < top_n) {
  warning("CV를 가진 유전자가 2000개보다 적습니다. 가능한 전부를 사용합니다.")
  top_n <- length(cv_keep)
}
top_genes <- names(sort(cv_keep, decreasing = TRUE))[1:top_n]

expr_hvg <- expr_mal_log[top_genes, ]

cat("HVG matrix dim (genes x malignant cells):",
    paste(dim(expr_hvg), collapse = " x "), "\n\n")

## 9. 결과 저장 ---------------------------------------------------------------

rdata_out <- file.path(processed_dir,
                       "K27M_malignant_HVG2000_preprocessed.RData")
txt_out   <- file.path(processed_dir,
                       "K27M_malignant_HVG2000_expr.txt")

save(
  expr_hvg,
  expr_mal_log,
  expr_mal,
  meta_mal,
  tsne_mal,
  file = rdata_out
)

write.table(
  expr_hvg,
  file = txt_out,
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

cat("=== 전처리 완료 ===\n")
cat("저장 파일:\n")
cat(" - ", rdata_out, "\n")
cat(" - ", txt_out,   "\n")
############################################################
