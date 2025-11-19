#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
04_compute_controllability_python.py

netctrl이 없을 때, Python + networkx 로
Liu et al. 구조적 controllability를 계산해서

  - driver node 리스트
  - critical / redundant / ordinary edge 리스트

를 생성하고, netctrl이 뱉는 것과 같은 파일 이름으로 저장한다.

사용 파일 (data/processed/):
  K27M_netctrl_nodes_thr0p7.tsv           (id, gene)
  K27M_netctrl_edges_thr0p7_ids.txt       (u, v)

생성 파일 (data/processed/):
  K27M_netctrl_driver_nodes_thr0p7.txt
  K27M_netctrl_critical_edges_thr0p7.txt
  K27M_netctrl_redundant_edges_thr0p7.txt
  K27M_netctrl_ordinary_edges_thr0p7.txt
"""

import os
import copy
import pandas as pd
import networkx as nx
from networkx.algorithms import bipartite as nx_bip

# ---------------- 0. 경로 설정 ---------------- #

def get_script_dir():
    return os.path.dirname(os.path.abspath(__file__))

script_dir   = get_script_dir()
scripts_dir  = os.path.dirname(script_dir)
project_root = os.path.dirname(scripts_dir)

processed_dir = os.path.join(project_root, "data", "processed")

print("Project root   :", project_root)
print("Processed dir  :", processed_dir)
print()

threshold  = 0.7
thr_label  = str(threshold).replace(".", "p")

nodes_path = os.path.join(
    processed_dir, f"K27M_netctrl_nodes_thr{thr_label}.tsv"
)
edges_path = os.path.join(
    processed_dir, f"K27M_netctrl_edges_thr{thr_label}_ids.txt"
)

if not os.path.exists(nodes_path):
    raise FileNotFoundError(f"노드 파일이 없습니다: {nodes_path}")
if not os.path.exists(edges_path):
    raise FileNotFoundError(f"엣지 파일이 없습니다: {edges_path}")

# ---------------- 1. 노드/엣지 로드 ---------------- #

nodes = pd.read_csv(nodes_path, sep="\t")
edges = pd.read_csv(edges_path, sep="\t", header=None, names=["u", "v"])

n_nodes = nodes.shape[0]
n_edges = edges.shape[0]

print(f"노드 수: {n_nodes}")
print(f"엣지 수: {n_edges}")
print()

# ---------------- 2. bipartite 그래프 구성 ---------------- #
# 왼쪽 U: u1..un, 오른쪽 V: v1..vn
# 디렉티드 그래프의 edge i->j 를 bipartite edge (ui, vj) 로 넣는다.

B = nx.Graph()
U = [f"u{i}" for i in range(1, n_nodes + 1)]
V = [f"v{i}" for i in range(1, n_nodes + 1)]

B.add_nodes_from(U, bipartite=0)
B.add_nodes_from(V, bipartite=1)

for _, row in edges.iterrows():
    ui = f"u{int(row.u)}"
    vj = f"v{int(row.v)}"
    B.add_edge(ui, vj)

print("bipartite graph nodes:", B.number_of_nodes())
print("bipartite graph edges:", B.number_of_edges())
print()

# ---------------- 3. 최대 매칭 & driver nodes ---------------- #

# networkx 버전에 따라 matching.maximum_matching 이 없을 수 있으므로
# bipartite.maximum_matching 사용
M0 = nx_bip.maximum_matching(B, top_nodes=U)
max_matching_edges = len(M0) // 2  # 매칭된 edge 개수

matched_V = {n for n in M0 if n in V}
num_matched_right = len(matched_V)

# unmatched right nodes = driver nodes
drivers = [i for i in range(1, n_nodes + 1) if f"v{i}" not in matched_V]
num_drivers = len(drivers)

nD = num_drivers / n_nodes

print(f"최대 매칭 edge 수: {max_matching_edges}")
print(f"driver node 수: {num_drivers}")
print(f"n_D = {nD:.4f}")
print()

# driver node 파일 저장
driver_path = os.path.join(
    processed_dir, f"K27M_netctrl_driver_nodes_thr{thr_label}.txt"
)
with open(driver_path, "w") as f:
    for i in drivers:
        f.write(f"{i}\n")

# ---------------- 4. edge 분류 ---------------- #
# critical: edge 제거 시 최대 매칭 크기가 줄어드는 edge
# 그 외 edge에 대해:
#   - maxcardinality + max_weight_matching으로
#     최대 매칭 중 이 edge를 포함하는 매칭이 있으면 ordinary
#     없으면 redundant

critical_edges = []
ordinary_edges = []
redundant_edges = []

B0 = B  # 원본 공유

print("edge classification 시작... (조금 시간 걸릴 수 있음)\n")

for idx, row in edges.iterrows():
    u_id = int(row.u)
    v_id = int(row.v)
    ui = f"u{u_id}"
    vj = f"v{v_id}"

    # --- 4-1. critical 여부 체크: edge 제거 후 최대 매칭 크기 비교 --- #
    B_minus = copy.deepcopy(B0)
    if B_minus.has_edge(ui, vj):
        B_minus.remove_edge(ui, vj)

    M_minus = nx_bip.maximum_matching(B_minus, top_nodes=U)
    max_edges_minus = len(M_minus) // 2

    if max_edges_minus < max_matching_edges:
        critical_edges.append((u_id, v_id))
        continue  # critical이면 나머지 테스트는 패스

    # --- 4-2. ordinary vs redundant --- #
    # e 를 최대 매칭에 포함하려고 weight = 2, 나머지 edge는 1
    B_weight = copy.deepcopy(B0)
    for (a, b) in B_weight.edges():
        B_weight[a][b]["weight"] = 1.0
    if B_weight.has_edge(ui, vj):
        B_weight[ui][vj]["weight"] = 2.0

    Mw = nx.algorithms.matching.max_weight_matching(
        B_weight, maxcardinality=True, weight="weight"
    )

    # Mw 는 edge 튜플들의 set
    # cardinality가 원래 최대 매칭과 같고, e가 포함되면 ordinary
    if len(Mw) == max_matching_edges and (
        (ui, vj) in Mw or (vj, ui) in Mw
    ):
        ordinary_edges.append((u_id, v_id))
    else:
        redundant_edges.append((u_id, v_id))

# ---------------- 5. 결과 저장 ---------------- #

crit_path = os.path.join(
    processed_dir, f"K27M_netctrl_critical_edges_thr{thr_label}.txt"
)
red_path = os.path.join(
    processed_dir, f"K27M_netctrl_redundant_edges_thr{thr_label}.txt"
)
ord_path = os.path.join(
    processed_dir, f"K27M_netctrl_ordinary_edges_thr{thr_label}.txt"
)

def write_edge_list(path, edge_list):
    with open(path, "w") as f:
        for u_id, v_id in edge_list:
            f.write(f"{u_id}\t{v_id}\n")

write_edge_list(crit_path, critical_edges)
write_edge_list(red_path, redundant_edges)
write_edge_list(ord_path, ordinary_edges)

print("=== 구조적 controllability 계산 완료 ===")
print("driver nodes 파일 :", driver_path)
print("critical edges   :", crit_path, f"(n={len(critical_edges)})")
print("redundant edges  :", red_path, f"(n={len(redundant_edges)})")
print("ordinary edges   :", ord_path, f"(n={len(ordinary_edges)})")
