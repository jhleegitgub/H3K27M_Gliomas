#!/usr/bin/env Rscript

# 09_make_summary_figures.R
# - Fig 1: 네트워크 구조 + overlap(Jaccard heatmap)
# - Fig 2: nD 비교 (malignant vs OPC/AC/OC)
# - (옵션) Fig 2B: module driver fraction barplot

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(purrr)

theme_set(theme_bw())

# ------------------------------------------------------------------
# 프로젝트 경로 설정 (다른 스크립트랑 동일하게)
# 현재 위치: H3K27M_Gliomas/scripts/03_Analysis 라고 가정
# ------------------------------------------------------------------
project_root       <- "/home/jihwanlee/research/H3K27M_Gliomas"
data_processed_dir <- file.path(project_root, "data", "processed")
results_fig_dir    <- file.path(project_root, "results", "figures")
results_tbl_dir    <- file.path(project_root, "results", "tables")

cat("Project root      :", project_root, "\n")
cat("Processed data dir:", data_processed_dir, "\n")
cat("Results fig dir   :", results_fig_dir, "\n")
cat("Results tbl dir   :", results_tbl_dir, "\n\n")

if (!dir.exists(results_fig_dir)) dir.create(results_fig_dir, recursive = TRUE)

# ------------------------------------------------------------------
# 1. 네트워크 구조 요약 (n_nodes, n_edges, density, mean degree)
#    - malignant + OPC/AC/OC 비교 (Fig 1A)
# ------------------------------------------------------------------

net_all <- read.table(
  file.path(results_tbl_dir, "K27M_gene_network_summary_thr0p7.tsv"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

net_grp <- read.table(
  file.path(results_tbl_dir, "K27M_gene_network_group_summary_thr0p7.tsv"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

net_all$network <- "Malignant"

# group 컬럼 이름 맞추기 (필요시)
if (!"network" %in% colnames(net_grp) && "group" %in% colnames(net_grp)) {
  net_grp <- net_grp %>% rename(network = group)
}

net_summary <- bind_rows(net_all, net_grp) %>%
  mutate(network = factor(network, levels = c("Malignant", "OPC", "AC", "OC")))

# ---- Fig 1A-1: n_nodes / n_edges barplot ----
p_nodes <- ggplot(net_summary, aes(x = network, y = n_nodes)) +
  geom_col() +
  labs(x = "", y = "Number of nodes", title = "Network size (nodes)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave(file.path(results_fig_dir, "Fig1A_nodes_by_network.pdf"),
       p_nodes, width = 4, height = 4)

p_edges <- ggplot(net_summary, aes(x = network, y = n_edges)) +
  geom_col() +
  labs(x = "", y = "Number of edges", title = "Network size (edges)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave(file.path(results_fig_dir, "Fig1A_edges_by_network.pdf"),
       p_edges, width = 4, height = 4)

# ---- Fig 1A-2: density / mean degree barplot ----
if ("density" %in% colnames(net_summary)) {
  p_density <- ggplot(net_summary, aes(x = network, y = density)) +
    geom_col() +
    labs(x = "", y = "Density", title = "Network density") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

  ggsave(file.path(results_fig_dir, "Fig1A_density_by_network.pdf"),
         p_density, width = 4, height = 4)
}

if ("mean_degree" %in% colnames(net_summary)) {
  p_mdeg <- ggplot(net_summary, aes(x = network, y = mean_degree)) +
    geom_col() +
    labs(x = "", y = "Mean degree", title = "Average degree") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

  ggsave(file.path(results_fig_dir, "Fig1A_mean_degree_by_network.pdf"),
         p_mdeg, width = 4, height = 4)
}

# ------------------------------------------------------------------
# 2. Jaccard overlap (nodes / edges) – Fig 1B, 1C
# ------------------------------------------------------------------

edge_files <- list(
  Malignant = file.path(data_processed_dir, "K27M_gene_network_edges_undirected_thr0p7.tsv"),
  OPC       = file.path(data_processed_dir, "K27M_group-OPC_gene_network_edges_undirected_thr0p7.tsv"),
  AC        = file.path(data_processed_dir, "K27M_group-AC_gene_network_edges_undirected_thr0p7.tsv"),
  OC        = file.path(data_processed_dir, "K27M_group-OC_gene_network_edges_undirected_thr0p7.tsv")
)

read_edge_list <- function(path) {
  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # source / target 컬럼 이름이 없으면 1,2번째 컬럼을 자동으로 source/target으로 바꿔줌
  if (!all(c("source", "target") %in% colnames(df))) {
    if (ncol(df) < 2) {
      stop("Edge file ", path, " 에서 최소 두 개의 컬럼이 필요합니다.")
    }
    old_names <- colnames(df)
    old_names[1:2] <- c("source", "target")
    colnames(df) <- old_names
  }

  df
}


edge_lists <- lapply(edge_files, read_edge_list)

node_sets <- lapply(edge_lists, function(df) unique(c(df$source, df$target)))

edge_sets <- lapply(edge_lists, function(df) {
  apply(df[, c("source", "target")], 1, function(x) paste(sort(x), collapse = "_"))
})

net_names <- names(edge_files)

jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  inter <- length(intersect(a, b))
  uni   <- length(union(a, b))
  if (uni == 0) return(NA_real_)
  inter / uni
}

jaccard_matrix_nodes <- matrix(NA, nrow = length(net_names), ncol = length(net_names),
                               dimnames = list(net_names, net_names))
jaccard_matrix_edges <- jaccard_matrix_nodes

for (i in seq_along(net_names)) {
  for (j in seq_along(net_names)) {
    jaccard_matrix_nodes[i, j] <- jaccard(node_sets[[i]], node_sets[[j]])
    jaccard_matrix_edges[i, j] <- jaccard(edge_sets[[i]], edge_sets[[j]])
  }
}

df_node_jacc <- as.data.frame(as.table(jaccard_matrix_nodes))
colnames(df_node_jacc) <- c("net1", "net2", "jaccard")

df_edge_jacc <- as.data.frame(as.table(jaccard_matrix_edges))
colnames(df_edge_jacc) <- c("net1", "net2", "jaccard")

p_node_jacc <- ggplot(df_node_jacc, aes(x = net1, y = net2, fill = jaccard)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", jaccard)), size = 3) +
  scale_fill_gradient(low = "white", high = "black", na.value = "grey90") +
  labs(x = "", y = "", fill = "Jaccard", title = "Node overlap (Jaccard index)") +
  theme_bw()

ggsave(file.path(results_fig_dir, "Fig1B_node_jaccard_heatmap.pdf"),
       p_node_jacc, width = 4, height = 4)

p_edge_jacc <- ggplot(df_edge_jacc, aes(x = net1, y = net2, fill = jaccard)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", jaccard)), size = 3) +
  scale_fill_gradient(low = "white", high = "black", na.value = "grey90") +
  labs(x = "", y = "", fill = "Jaccard", title = "Edge overlap (Jaccard index)") +
  theme_bw()

ggsave(file.path(results_fig_dir, "Fig1C_edge_jaccard_heatmap.pdf"),
       p_edge_jacc, width = 4, height = 4)

# ------------------------------------------------------------------
# 3. nD (driver fraction) 비교 – Fig 2A
# ------------------------------------------------------------------

nd_all <- read.table(
  file.path(results_tbl_dir, "K27M_netctrl_summary_thr0p7_withGamma.tsv"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)
nd_grp <- read.table(
  file.path(results_tbl_dir, "K27M_group_netctrl_summary_thr0p7.tsv"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

nd_all$network <- "Malignant"
if (!"network" %in% colnames(nd_grp) && "group" %in% colnames(nd_grp)) {
  nd_grp <- nd_grp %>% rename(network = group)
}

nd_summary <- bind_rows(
  nd_all %>% select(network, nD),
  nd_grp  %>% select(network, nD)
) %>%
  mutate(network = factor(network, levels = c("Malignant", "OPC", "AC", "OC")))

p_nd <- ggplot(nd_summary, aes(x = network, y = nD)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.2f", nD)), vjust = -0.3, size = 3) +
  ylim(0, max(nd_summary$nD, na.rm = TRUE) * 1.2) +
  labs(x = "", y = "Driver fraction (nD)",
       title = "Controllability (driver fraction)") +
  theme_bw()

ggsave(file.path(results_fig_dir, "Fig2A_nD_by_network.pdf"),
       p_nd, width = 4, height = 4)

# ------------------------------------------------------------------
# 4. (옵션) module gene driver fraction – Fig 2B
# ------------------------------------------------------------------

module_file <- file.path(results_tbl_dir, "K27M_module_summary_thr0p7.tsv")
if (file.exists(module_file)) {
  mod <- read.table(module_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # 1열이 module 이름이라고 가정
  if (!"module" %in% colnames(mod)) {
    colnames(mod)[1] <- "module"
  }

  # ---- 컬럼 이름 자동 탐색 ----
  # driver 비율 관련 컬럼 중에서 "in_set" / "other" 가 들어간 걸 찾는다.
  frac_cols <- grep("frac", colnames(mod), value = TRUE)

  col_in_set <- frac_cols[grepl("in_set", frac_cols)][1]
  col_other  <- frac_cols[grepl("other",  frac_cols)][1]

  if (is.na(col_in_set) || is.na(col_other)) {
    message("모듈 요약에서 driver fraction 컬럼을 찾지 못했습니다. Fig2B는 건너뜁니다.")
  } else {
    message("모듈 driver fraction 컬럼 사용: in_set = ", col_in_set,
            ", other = ", col_other)

    # 그림용 데이터프레임 구성
    df_plot <- data.frame(
      module         = mod$module,
      Module_genes   = mod[[col_in_set]],
      Rest_of_network = mod[[col_other]]
    )

    mod_long <- df_plot %>%
      tidyr::pivot_longer(
        cols      = c("Module_genes", "Rest_of_network"),
        names_to  = "set",
        values_to = "frac_driver"
      )

    p_mod <- ggplot(mod_long, aes(x = module, y = frac_driver, fill = set)) +
      geom_col(position = "dodge") +
      geom_text(aes(label = sprintf("%.2f", frac_driver)),
                position = position_dodge(width = 0.9),
                vjust = -0.3, size = 3) +
      ylim(0, max(mod_long$frac_driver, na.rm = TRUE) * 1.2) +
      labs(x = "Module", y = "Driver fraction",
           title = "Driver fraction in cell-cycle / PDGFRA modules") +
      theme_bw()

    ggsave(file.path(results_fig_dir, "Fig2B_module_driver_fraction.pdf"),
           p_mod, width = 5, height = 4)
  }
} else {
  message("module summary 파일이 없어 Fig2B는 건너뜁니다.")
}
