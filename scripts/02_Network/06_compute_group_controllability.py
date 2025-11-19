#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
06_compute_group_controllability.py

05_build_group_networks.R 가 만든 그룹별 netctrl 입력 파일을 읽어서
각 그룹(OPC / AC / OC)에 대해

  - driver node 리스트
  - n_D (n_driver / n_nodes)

를 계산하고 요약 테이블로 저장한다.

사용 파일 (data/processed):
  K27M_group-<grp>_netctrl_nodes_thr0p7.tsv
  K27M_group-<grp>_netctrl_edges_thr0p7_ids.txt

산출물:
  data/processed/K27M_group-<grp>_netctrl_driver_nodes_thr0p7.txt
  results/tables/K27M_group_netctrl_summary_thr0p7.tsv
"""

import os
import pandas as pd
import networkx as nx
from networkx.algorithms import bipartite as nx_bip


def get_project_root():
    script_path = os.path.abspath(__file__)
    script_dir = os.path.dirname(script_path)            # .../scripts/02_Network
    scripts_dir = os.path.dirname(script_dir)            # .../scripts
    project_root = os.path.dirname(scripts_dir)          # 프로젝트 루트
    return project_root


def compute_driver_for_group(group_label, thr_label="0p7"):
    project_root = get_project_root()
    processed_dir = os.path.join(project_root, "data", "processed")

    prefix = f"K27M_group-{group_label}"
    nodes_path = os.path.join(
        processed_dir, f"{prefix}_netctrl_nodes_thr{thr_label}.tsv"
    )
    edges_path = os.path.join(
        processed_dir, f"{prefix}_netctrl_edges_thr{thr_label}_ids.txt"
    )

    if not (os.path.exists(nodes_path) and os.path.exists(edges_path)):
        print(f"[{group_label}] 입력 파일이 없습니다. 건너뜀.")
        return None

    print(f"=== 그룹: {group_label} ===")
    print("Nodes file :", nodes_path)
    print("Edges file :", edges_path)

    nodes = pd.read_csv(nodes_path, sep="\t")
    edges = pd.read_csv(edges_path, sep="\t", header=None, names=["u", "v"])

    n_nodes = nodes.shape[0]
    n_edges = edges.shape[0]
    print(f"노드 수: {n_nodes}")
    print(f"엣지 수: {n_edges}")

    # bipartite 그래프 구성 (Liu 모델과 동일: U -> V)
    U = [f"u{i}" for i in range(1, n_nodes + 1)]
    V = [f"v{i}" for i in range(1, n_nodes + 1)]

    B = nx.Graph()
    B.add_nodes_from(U, bipartite=0)
    B.add_nodes_from(V, bipartite=1)

    for _, row in edges.iterrows():
        ui = f"u{int(row.u)}"
        vj = f"v{int(row.v)}"
        B.add_edge(ui, vj)

    print("bipartite graph nodes:", B.number_of_nodes())
    print("bipartite graph edges:", B.number_of_edges())

    # 최대 매칭 → driver node (unmatched right nodes)
    print("최대 매칭 계산 중 ...")
    M = nx_bip.maximum_matching(B, top_nodes=U)

    matched_V = {n for n in M if n in V}
    drivers = [i for i in range(1, n_nodes + 1) if f"v{i}" not in matched_V]

    n_driver = len(drivers)
    nD = n_driver / n_nodes if n_nodes > 0 else float("nan")

    print(f"driver node 수: {n_driver}")
    print(f"n_D = {nD:.4f}")
    print()

    # driver node 파일 저장 (1-based id)
    driver_path = os.path.join(
        processed_dir, f"{prefix}_netctrl_driver_nodes_thr{thr_label}.txt"
    )
    with open(driver_path, "w") as f:
        for i in drivers:
            f.write(f"{i}\n")

    print("driver node 파일:", driver_path)

    return {
        "group": group_label,
        "threshold": float(thr_label.replace("p", ".")),
        "n_nodes": n_nodes,
        "n_edges": n_edges,
        "n_driver": n_driver,
        "nD": nD,
    }


def main():
    project_root = get_project_root()
    results_tables_dir = os.path.join(project_root, "results", "tables")
    os.makedirs(results_tables_dir, exist_ok=True)

    thr_label = "0p7"
    groups = ["OPC", "AC", "OC"]

    summary_rows = []
    for g in groups:
        res = compute_driver_for_group(g, thr_label=thr_label)
        if res is not None:
            summary_rows.append(res)

    if not summary_rows:
        print("요약할 결과가 없습니다.")
        return

    summary_df = pd.DataFrame(summary_rows)

    summary_path = os.path.join(
        results_tables_dir, f"K27M_group_netctrl_summary_thr{thr_label}.tsv"
    )
    summary_df.to_csv(summary_path, sep="\t", index=False)

    print("\n=== 그룹별 controllability 요약 완료 ===")
    print("Summary table:", summary_path)


if __name__ == "__main__":
    main()
