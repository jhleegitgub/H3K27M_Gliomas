#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
H3K27M single-cell preprocessing
- Load RSEM matrix + metadata + hierarchy (+tSNE optional)
- Robust cell-ID harmonization (auto-detect metadata column by overlap)
- Convert hierarchy score columns to numeric and label OPC/AC/OC-like states
- Export:
    * data/processed/mean_expr_by_group.tsv.gz       (CellChat/NicheNet)
    * data/processed/pseudobulk_by_sample_state.tsv.gz (WGCNA)
    * results/tables/cell_counts_by_group.tsv
    * (optional) data/processed/k27m_anndata.h5ad

Usage
  python scripts/01_CellChat_NicheNet/01_preprocess.py
  python scripts/01_CellChat_NicheNet/01_preprocess.py --config ./config.json
"""

import os, sys, json, argparse
from pathlib import Path
import pandas as pd
import numpy as np

# ──────────────────────────────────────────────────────────────────────────────
# Utils
# ──────────────────────────────────────────────────────────────────────────────
def find_repo_root(start: Path) -> Path:
    """Walk up to find a directory containing config.json or .git; fallback to cwd."""
    cur = start.resolve()
    for _ in range(8):
        if (cur/"config.json").exists() or (cur/".git").exists():
            return cur
        if cur.parent == cur:
            break
        cur = cur.parent
    return start.resolve()

def load_table(path: Path) -> pd.DataFrame | None:
    if not path or not path.exists():
        print(f"[WARN] Missing file: {path}")
        return None
    sep = "\t" if path.suffix.lower() != ".csv" else ","
    try:
        if path.suffix.lower().endswith("gz"):
            return pd.read_csv(path, sep=sep, compression="gzip")
        return pd.read_csv(path, sep=sep)
    except Exception:
        return pd.read_table(path)

def ensure_gene_rows(df: pd.DataFrame, gene_col: str, cell_axis: str) -> pd.DataFrame:
    """Return expression with genes as rows, cells as columns."""
    if df is None: return None
    if cell_axis.lower().startswith("col"):
        if gene_col not in df.columns:
            df = df.rename(columns={df.columns[0]: "gene"})
            gene_col = "gene"
        return df.set_index(gene_col)
    # cells are rows
    if gene_col in df.columns:
        return df.set_index(gene_col).T
    return df.set_index(df.columns[0]).T

def guess_cell_col_by_overlap(meta: pd.DataFrame, expr_cols) -> str:
    """Pick metadata column holding cell IDs by max overlap with expr columns."""
    candidates = ["cell","Cell","cell_id","CellID","Barcode","barcode",
                  "id","ID","name","Name","Cell.Name","CellName"]
    for c in candidates:
        if c in meta.columns: return c
    expr_set = set(map(str, expr_cols))
    best_col, best_n = None, -1
    for c in meta.columns:
        vals = set(map(str, meta[c].astype(str)))
        n = len(expr_set & vals)
        if n > best_n:
            best_col, best_n = c, n
    if best_col is not None and best_n > 0:
        print(f"[INFO] Inferred cell-id column: '{best_col}' (overlap={best_n})")
        return best_col
    print("[WARN] Could not infer a cell-id column reliably; using the first column.")
    return meta.columns[0]

def attach_hierarchy(meta: pd.DataFrame, hier: pd.DataFrame, cell_col: str, cells) -> pd.DataFrame:
    if hier is None:
        print("[INFO] No hierarchy file; state labeling may be limited.")
        return meta
    hcell = guess_cell_col_by_overlap(hier, cells)
    hier = hier.drop_duplicates(subset=[hcell]).set_index(hcell).reindex(cells).reset_index()
    return meta.merge(hier, left_on=cell_col, right_on=hcell, how="left")

def label_state_from_scores(row: pd.Series) -> str:
    """Choose max among OPC/AC/OC keywords; ignore non-numeric values."""
    scores = {}
    for k, v in row.items():
        if not isinstance(k, str): 
            continue
        lk = k.lower()
        if not any(t in lk for t in ["opc","ac","oc"]):
            continue
        try:
            val = float(v)
            if np.isnan(val): 
                continue
        except Exception:
            continue
        if "opc" in lk:
            scores["OPC-like"] = max(scores.get("OPC-like", -1e9), val)
        if "ac" in lk and "fac" not in lk:
            scores["AC-like"] = max(scores.get("AC-like", -1e9), val)
        if "oc" in lk and "doc" not in lk:
            scores["OC-like"] = max(scores.get("OC-like", -1e9), val)
    return max(scores, key=scores.get) if scores else "Unknown"

# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=str, default=None, help="path to config.json (optional)")
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    repo = find_repo_root(here)
    cfg_path = Path(args.config) if args.config else (repo / "config.json")
    if not cfg_path.exists():
        print(f"[ERROR] config.json not found at {cfg_path}. Provide with --config.")
        sys.exit(1)

    CFG = json.load(open(cfg_path, "r", encoding="utf-8"))
    def P(k: str) -> Path:
        v = CFG.get(k)
        return (repo / v) if v else None

    # output dirs
    (repo / "data/processed").mkdir(parents=True, exist_ok=True)
    (repo / "results/tables").mkdir(parents=True, exist_ok=True)

    # load
    expr = load_table(P("rsem_matrix_path"))
    meta = load_table(P("metadata_path"))
    hier = load_table(P("hierarchy_scores_path")) if CFG.get("hierarchy_scores_path") else None
    _ = load_table(P("tsne_path")) if CFG.get("tsne_path") else None

    if expr is None or meta is None:
        print("[ERROR] Expression or metadata file missing—check config.json paths.")
        sys.exit(1)

    expr = ensure_gene_rows(expr, CFG.get("gene_id_column","gene"), CFG.get("cell_id_axis","columns"))
    if expr is None:
        print("[ERROR] Could not standardize expression matrix (genes x cells).")
        sys.exit(1)

    # force expression numeric just in case
    expr = expr.apply(pd.to_numeric, errors="coerce")

    # harmonize cells
    expr.columns = expr.columns.map(str).str.strip()
    meta = meta.copy()
    cell_col = guess_cell_col_by_overlap(meta, expr.columns)

    cells = list(set(expr.columns).intersection(set(meta[cell_col].astype(str).str.strip())))
    if len(cells) == 0:
        raise RuntimeError("[ERROR] No overlapping cell IDs between expression and metadata.")
    cells_sorted = sorted(cells)

    expr = expr.loc[:, cells_sorted]
    meta[cell_col] = meta[cell_col].astype(str).str.strip()
    meta = (meta.drop_duplicates(subset=[cell_col])
                .set_index(cell_col).loc[cells_sorted].reset_index())

    # attach hierarchy and label states
    meta = attach_hierarchy(meta, hier, meta.columns[0], cells_sorted)
    if hier is not None:
        # convert potential score columns to numeric first
        score_like_cols = [
            c for c in meta.columns
            if isinstance(c, str) and any(k in c.lower()
               for k in ["opc","ac-like","oc-like","opc-variable","stem"])
        ]
        for c in score_like_cols:
            meta[c] = pd.to_numeric(meta[c], errors="coerce")
        meta["StateLabel"] = meta.apply(label_state_from_scores, axis=1)
        state_col = "StateLabel"
    else:
        state_col = None
        for c in ["State","state","cluster","CellType","celltype"]:
            if c in meta.columns: state_col = c; break
        if state_col is None:
            meta["StateLabel"] = "Unknown"; state_col = "StateLabel"

    # malignant flag
    mal_col = None
    for c in ["Malignant","malignant","Tumor","is_malignant"]:
        if c in meta.columns: mal_col = c; break
    if mal_col is None:
        meta["is_malignant"] = ~meta[state_col].str.contains(
            "microglia|oligodendro|astro", case=False, na=False
        )
        mal_col = "is_malignant"

    # group columns: sample type + state
    group_cols = []
    for c in ["SampleType","sample_type","Type","Model","ModelType","source"]:
        if c in meta.columns:
            group_cols.append(c); break
    group_cols.append(state_col)

    # mean expr per (sampleType x state) → CellChat/NicheNet
    expr_T = expr.T
    expr_T.index = cells_sorted
    expr_T[group_cols] = meta[group_cols].values
    mean_expr = expr_T.groupby(group_cols).mean(numeric_only=True).T
    mean_path = repo / "data/processed/mean_expr_by_group.tsv.gz"
    mean_expr.to_csv(mean_path, sep="\t", compression="gzip")
    print(f"[OK] wrote {mean_path}")

    # pseudobulk per (sample x state) → WGCNA
    sample_col = None
    for c in ["Sample","sample","Patient","patient","Specimen","specimen","SampleID","Case"]:
        if c in meta.columns: sample_col = c; break
    if sample_col is None: sample_col = group_cols[0]

    expr_T["__Sample"] = meta[sample_col].values
    expr_T["__State"]  = meta[state_col].values
    pseudobulk = expr_T.groupby(["__Sample","__State"]).mean(numeric_only=True).T
    pbulk_path = repo / "data/processed/pseudobulk_by_sample_state.tsv.gz"
    pseudobulk.to_csv(pbulk_path, sep="\t", compression="gzip")
    print(f"[OK] wrote {pbulk_path}")

    # counts table
    counts = meta.groupby([sample_col, state_col]).size().reset_index(name="n_cells")
    counts_path = repo / "results/tables/cell_counts_by_group.tsv"
    counts.to_csv(counts_path, sep="\t", index=False)
    print(f"[OK] wrote {counts_path}")

    # optional: AnnData
    try:
        import anndata as ad, scipy.sparse as sp
        X = sp.csr_matrix(np.asarray(expr.values, dtype=np.float32))
        adata = ad.AnnData(X=X, var=pd.DataFrame(index=expr.index), obs=meta.set_index(meta.columns[0]))
        adata_path = repo / "data/processed/k27m_anndata.h5ad"
        adata.write_h5ad(str(adata_path))
        print(f"[OK] wrote {adata_path}")
    except Exception as e:
        print("[WARN] skipped .h5ad export (install anndata/scipy to enable).", e)

if __name__ == "__main__":
    main()
