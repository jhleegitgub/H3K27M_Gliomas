############################################################
# 06_make_figures.R
#
#  - 전처리(HVG), 네트워크, controllability 요약을 읽어서
#    네트워크바이올로지 프로젝트용 그림/표 생성
#
#   1) degree 분포 (histogram + log-log plot)
#   2) 상위 hub/driver gene barplot
#   3) degree tail로부터 대략적인 gamma 추정
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
results_fig_dir    <- file.path(project_root, "results", "figures")

if (!dir.exists(results_fig_dir)) {
  dir.create(results_fig_dir, recursive = TRUE, showWarnings = FALSE)
}

cat("Project root   :", project_root,        "\n")
cat("Processed dir  :", processed_dir,      "\n")
cat("Figures dir    :", results_fig_dir,    "\n\n")

threshold <- 0.7
thr_label <- gsub("\\.", "p", as.character(threshold))

## 1. 데이터 로드 -------------------------------------------------------------

# 1-1) 전처리 RData (expr_hvg, meta_mal, tsne_mal 등)
rdata_path <- file.path(
  processed_dir,
  "K27M_malignant_HVG2000_preprocessed.RData"
)
if (!file.exists(rdata_path)) {
  stop("전처리 RData를 찾을 수 없습니다: ", rdata_path)
}
load(rdata_path)
# expr_hvg, expr_mal_log, expr_mal, meta_mal, tsne_mal

# 1-2) gene-level controllability 요약
gene_summary_path <- file.path(
  results_tables_dir,
  paste0("K27M_netctrl_gene_summary_thr", thr_label, ".tsv")
)
if (!file.exists(gene_summary_path)) {
  stop("gene_summary 파일이 없습니다: ", gene_summary_path)
}
gene_summary <- read.table(
  gene_summary_path,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

summary_path <- file.path(
  results_tables_dir,
  paste0("K27M_netctrl_summary_thr", thr_label, ".tsv")
)
net_summary <- read.table(
  summary_path,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

print(net_summary)
cat("\n")

## 2. degree 분포 그림 --------------------------------------------------------

deg <- gene_summary$deg_total

# 2-1) 일반 histogram
pdf(file.path(results_fig_dir, "K27M_degree_histogram_thr0p7.pdf"),
    width = 6, height = 5)
hist(deg,
     breaks = 50,
     main = "Degree distribution (total degree)",
     xlab = "Degree k",
     ylab = "Frequency")
dev.off()

# 2-2) log-log degree 분포 (P(k) vs k)
deg_tab <- table(deg)
k_vals  <- as.numeric(names(deg_tab))
pk      <- as.numeric(deg_tab) / sum(deg_tab)

# k>=1 만 사용
keep    <- k_vals >= 1
k_tail  <- k_vals[keep]
pk_tail <- pk[keep]

pdf(file.path(results_fig_dir, "K27M_degree_loglog_thr0p7.pdf"),
    width = 6, height = 5)
plot(k_tail, pk_tail,
     log = "xy",
     xlab = "k (log scale)",
     ylab = "P(k) (log scale)",
     main = "Degree distribution (log-log)",
     pch = 16)
dev.off()

## 3. power-law tail에서 gamma 대충 추정 ------------------------------------
#   log P(k) = -gamma * log k + c  → gamma ~ -slope
#   너무 작은 degree는 제외: k >= k_min (예: 3)

k_min <- 3
keep_tail <- k_tail >= k_min
x <- log10(k_tail[keep_tail])
y <- log10(pk_tail[keep_tail])

fit <- lm(y ~ x)

gamma_hat <- -coef(fit)[2]
cat("Estimated gamma (tail, k >=", k_min, "):", as.numeric(gamma_hat), "\n\n")

# 회귀선까지 같이 그림 하나 더
pdf(file.path(results_fig_dir, "K27M_degree_loglog_tail_fit_thr0p7.pdf"),
    width = 6, height = 5)
plot(k_tail, pk_tail,
     log = "xy",
     xlab = "k (log scale)",
     ylab = "P(k) (log scale)",
     main = paste0("Degree tail fit (k>=", k_min, ")\n gamma ~ ",
                   round(gamma_hat, 2)),
     pch = 16)
# 회귀선: y = a + b x
a <- coef(fit)[1]
b <- coef(fit)[2]
xx <- seq(min(x), max(x), length.out = 100)
yy <- a + b * xx
lines(10^xx, 10^yy, col = "red", lwd = 2)
dev.off()

## 4. 상위 hub / driver gene barplot -----------------------------------------

# 상위 30개 high-degree gene
top_n <- 30
ord   <- order(gene_summary$deg_total, decreasing = TRUE)
top_idx <- ord[1:top_n]
top_genes <- gene_summary[top_idx, ]

# driver 여부를 색으로 구분할 수 있게 factor로 변환
top_genes$driver_flag <- ifelse(top_genes$is_driver == 1, "driver", "non-driver")

pdf(file.path(results_fig_dir, "K27M_top30_degree_genes_thr0p7.pdf"),
    width = 8, height = 5)
par(mar = c(8, 4, 4, 2))
barplot(height = top_genes$deg_total,
        names.arg = top_genes$gene,
        las = 2,
        main = "Top 30 genes by total degree",
        ylab = "Total degree",
        col = ifelse(top_genes$is_driver == 1, "red", "grey"))
legend("topright", legend = c("driver", "non-driver"),
       fill = c("red", "grey"), cex = 0.8)
dev.off()

# Critical-edge source 여부까지 보고 싶은 경우: scatter plot
pdf(file.path(results_fig_dir, "K27M_degree_vs_flags_thr0p7.pdf"),
    width = 6, height = 5)
plot(gene_summary$deg_total,
     gene_summary$is_crit_source,
     xlab = "Total degree",
     ylab = "Is critical-source (0/1)",
     main = "Critical-source vs degree",
     pch = 16)
dev.off()

## 5. gamma / 네트워크 요약을 테이블에 추가 저장 ----------------------------

net_summary$gamma_tail <- as.numeric(gamma_hat)
net_summary$k_min      <- k_min

summary_with_gamma_path <- file.path(
  results_tables_dir,
  paste0("K27M_netctrl_summary_thr", thr_label, "_withGamma.tsv")
)

write.table(
  net_summary,
  file      = summary_with_gamma_path,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

cat("=== 그림 및 요약 생성 완료 ===\n")
cat("Figures saved in:", results_fig_dir, "\n")
cat("Summary+gamma   :", summary_with_gamma_path, "\n")
############################################################
