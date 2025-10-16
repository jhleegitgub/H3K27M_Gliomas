# ---- CellChat with group-mean expression ------------------------------------
suppressPackageStartupMessages({
  library(CellChat)
  library(data.table)
  library(ComplexHeatmap)
})

# repo root 기준 경로
expr_path   <- "data/processed/mean_expr_by_group.tsv.gz"
counts_path <- "results/tables/cell_counts_by_group.tsv"
out_dir_f   <- "results/figures"
out_dir_t   <- "results/tables"
dir.create(out_dir_f, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_t, recursive = TRUE, showWarnings = FALSE)

# 1) 입력 불러오기 --------------------------------------------------------------
expr <- fread(expr_path) |> as.data.frame()
rownames(expr) <- expr[[1]]; expr[[1]] <- NULL
# expr: genes x (SampleType|State) 그룹 열들

counts <- fread(counts_path) |> as.data.frame()
# counts: Sample × State × n_cells
# 그룹 키를 CellChat의 "group"으로 사용
if (!"group" %in% colnames(counts)) {
  key_cols <- colnames(counts)[1:2]                     # sample, state
  counts$group <- apply(counts[, key_cols], 1, paste, collapse="|")
}

# expr의 열 이름도 동일 키 형식으로 맞추기
# (전처리 스크립트에서 group_by(c("SampleType(or Type)", "StateLabel"))로 만들었음)
# 열 이름이 "Type|StateLabel" 꼴이어야 아래 병합이 잘 맞음
grp_cols <- colnames(expr)
stopifnot(!anyDuplicated(grp_cols))

# 2) CellChat 오브젝트 구성 -----------------------------------------------------
# (중요) CellChat은 gene x cell 행렬을 기대하지만, 여기서는
# '각 그룹을 하나의 대표 세포'로 간주한다. 그룹 크기는 population.size로 가중치.
data.input <- as.matrix(expr)               # genes x groups
meta <- data.frame(labels = colnames(data.input), row.names = colnames(data.input))
colnames(meta) <- "group"

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

CellChatDB <- CellChatDB.human
# 분비 신호만 먼저 (원하면 'ECM-Receptor','Cell-Cell Contact'도 가능)
CellChatDB.use <- CellChatDB[CellChatDB$annotation %in% c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"), ]
cellchat@DB <- CellChatDB.use

# CellChat 기본 전처리
cellchat <- subsetData(cellchat)                 # DB에 없는 유전자는 자동 제거
cellchat <- identifyOverExpressedGenes(cellchat) # (그룹 평균이라도 동작함)
cellchat <- identifyOverExpressedInteractions(cellchat)

# 그룹 크기 가중치(= 실제 세포 수)
pop_size <- tapply(counts$n_cells, counts$group, sum)
pop_size <- pop_size[colnames(cellchat@data)]    # 정렬
pop_size[is.na(pop_size)] <- 1

cellchat <- computeCommunProb(cellchat, population.size = pop_size, raw.use = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 1)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)               # 신호 경로 단위로 집계

# 3) 결과 저장 ------------------------------------------------------------------
# (i) LR 쌍별 엣지 테이블
df.net  <- subsetCommunication(cellchat)                      # 모든 LR 엣지
df.path <- subsetCommunication(cellchat, slot.name = "netP")  # 경로 단위

write.table(df.net,  file.path(out_dir_t, "cellchat_edges_by_LR.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)
write.table(df.path, file.path(out_dir_t, "cellchat_edges_by_pathway.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

# (ii) 네트워크 강도 행렬(발신/수신)
group.names <- levels(cellchat@idents)
mat.send <- netAnalysis_computeCentrality(cellchat, slot.name = "net", type = "outgoing")$centrality
mat.recv <- netAnalysis_computeCentrality(cellchat, slot.name = "net", type = "incoming")$centrality

# (iii) 그림 몇 개
pdf(file.path(out_dir_f, "cellchat_overview.pdf"), width=10, height=8)
netVisual_circle(cellchat@net$count, vertex.weight=pop_size, weight.scale=T, label.edge=F, title.name="Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight=pop_size, weight.scale=T, label.edge=F, title.name="Interaction weights")
dev.off()

pdf(file.path(out_dir_f, "cellchat_role_heatmap.pdf"), width=9, height=6)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
dev.off()

# 특정 경로(예: PDGF)만 네트워크 그림
if ("PDGF" %in% cellchat@netP$pathways) {
  pdf(file.path(out_dir_f, "cellchat_PDGF_network.pdf"), width=7, height=6)
  netVisual_aggregate(cellchat, signaling = "PDGF", layout = "circle")
  dev.off()
}

message("[OK] CellChat finished. Tables in '", out_dir_t, "', figures in '", out_dir_f, "'.")

