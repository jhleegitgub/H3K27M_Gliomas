# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # GUI 없이 파일로 저장
import matplotlib.pyplot as plt
import networkx as nx
from pathlib import Path

OUTDIR = Path("./_out"); OUTDIR.mkdir(exist_ok=True)

# ---------------- Load expression data ----------------
expr_path = "K27Mproject.RSEM.vh20170621.txt"   # 또는 GSE102130_...txt
data_table = pd.read_table(expr_path, sep="\t", index_col=0)
data_matrix = data_table.to_numpy(dtype=float)
cell_id = data_table.columns.to_numpy().astype(str)
gene_name = data_table.index.to_numpy().astype(str)
log_data = np.log2(data_matrix + 1.0)

# ---------------- Load annotation data (cell_type.txt) ----------------
# cell_type.txt는 "cell_id \t type" 두 열이라고 가정
ct = pd.read_table("cell_type.txt", sep="\t", header=None, names=["cell_id","type"])
ct = ct.set_index("cell_id").reindex(cell_id)       # 표현행렬의 컬럼 순서에 맞춤
cell_type = ct["type"].fillna("NA").to_numpy()      # ← 이게 최종 라벨 벡터

# ---------------- Gene filtering ----------------
existenceCutoff = 0.1
geneFilter = (np.sum(log_data > 1, axis=1) / log_data.shape[1]) > existenceCutoff
log_data = log_data[geneFilter, :]
gene_name = gene_name[geneFilter]

eps = 1e-12

# ---------------- Correlation network (genes) ----------------
mean_g = np.mean(log_data, axis=1); std_g = np.std(log_data, axis=1)
cv_g = std_g / (mean_g + eps)
order_g = np.argsort(cv_g)[::-1]
highCutoff_g = min(2000, log_data.shape[0])
sel_g = order_g[:highCutoff_g]

corr_g = np.corrcoef(log_data[sel_g, :])
np.fill_diagonal(corr_g, 0.0)
corrCutoff_g = 0.6
corr_g[np.abs(corr_g) < corrCutoff_g] = 0.0

Gg = nx.to_networkx_graph(corr_g, create_using=nx.Graph)
Gg = nx.relabel_nodes(Gg, {i: gene_name[sel_g[i]] for i in range(highCutoff_g)})

if Gg.number_of_nodes() and Gg.number_of_edges():
    lcc_nodes = max(nx.connected_components(Gg), key=len)
    Gg = Gg.subgraph(lcc_nodes).copy()

nx.write_weighted_edgelist(
    Gg, OUTDIR / f"corrNetwork4geneCVhigh{highCutoff_g}pcc{corrCutoff_g}.txt", delimiter="\t"
)
plt.figure(figsize=(10,10)); nx.draw_networkx(Gg, with_labels=False, node_size=10)
plt.tight_layout(); plt.savefig(OUTDIR / f"corrNetwork4geneCVhigh{highCutoff_g}_pcc{corrCutoff_g}.png", dpi=300); plt.close()

# ---------------- Correlation network (cells) ----------------
mean_c = np.mean(log_data, axis=0); std_c = np.std(log_data, axis=0)
cv_c = std_c / (mean_c + eps)
order_c = np.argsort(cv_c)[::-1]
highCutoff_c = min(200, log_data.shape[1])
sel_c = order_c[:highCutoff_c]

corr_c = np.corrcoef(log_data[:, sel_c].T)
np.fill_diagonal(corr_c, 0.0)
corrCutoff_c = 0.35
corr_c[np.abs(corr_c) < corrCutoff_c] = 0.0

Gc = nx.to_networkx_graph(corr_c, create_using=nx.Graph)
Gc = nx.relabel_nodes(Gc, {i: cell_id[sel_c[i]] for i in range(highCutoff_c)})

if Gc.number_of_nodes() and Gc.number_of_edges():
    lcc_nodes = max(nx.connected_components(Gc), key=len)
    Gc = Gc.subgraph(lcc_nodes).copy()

nx.write_weighted_edgelist(
    Gc, OUTDIR / f"corrNetwork4cellCVhigh{highCutoff_c}pcc{corrCutoff_c}.txt", delimiter="\t"
)
plt.figure(figsize=(10,10)); nx.draw_networkx(Gc, with_labels=False, node_size=10)
plt.tight_layout(); plt.savefig(OUTDIR / f"corrNetwork4cellCVhigh{highCutoff_c}_pcc{corrCutoff_c}.png", dpi=300); plt.close()

print("✅ 완료: 결과는 _out/ 폴더에 저장되었어요.")
