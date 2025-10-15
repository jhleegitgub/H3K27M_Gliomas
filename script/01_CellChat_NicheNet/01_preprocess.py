import os, json, gzip, sys
from pathlib import Path
import pandas as pd
import numpy as np

BASE = Path(__file__).resolve().parents[1]
CFG = json.load(open(BASE/"config.json"))
RAW = BASE/"data/raw"
OUTP = BASE/"data/processed"
RES = BASE/"results"
OUTP.mkdir(parents=True, exist_ok=True)
(RES/"tables").mkdir(parents=True, exist_ok=True)

def load_table(path):
    p = Path(path)
    if not p.exists():
        print(f"[WARN] Missing file: {p}")
        return None
    if p.suffix == ".gz":
        return pd.read_csv(p, sep="\t", compression="gzip")
    else:
        # Auto-detect delimiter
        try:
            return pd.read_csv(p, sep="\t")
        except Exception:
            return pd.read_csv(p)

def ensure_gene_rows(df, gene_col, cell_axis):
    """
    Return expression with genes as rows, cells as columns.
    """
    if df is None:
        return None
    if cell_axis == "columns":
        # Expect first column = gene, rest = cells
        if gene_col not in df.columns:
            # try the first column as gene
            df = df.rename(columns={df.columns[0]: "gene"})
            gene_col = "gene"
        df = df.set_index(gene_col)
        return df
    else:
        # cells are rows; transpose
        if gene_col not in df.columns:
            # if genes are index
            df = df.set_index(df.columns[0])
        df = df.T
        return df

# 1) Load files
expr = load_table(CFG["rsem_matrix_path"])
meta = load_table(CFG["metadata_path"])
hier = load_table(CFG["hierarchy_scores_path"])
tsne = load_table(CFG["tsne_path"])

if expr is None or meta is None:
    print("[ERROR] Expression or metadata file missing. Check data/raw and config.json.")
    sys.exit(0)

expr = ensure_gene_rows(expr, CFG["gene_id_column"], CFG["cell_id_axis"])

# Try to detect if genes are rows (index) and cells columns
if expr is None:
    print("[ERROR] Could not standardize expression table shape.")
    sys.exit(0)

# 2) Harmonize cell IDs between expr and meta
# Find intersection of columns (cells)
meta_cell_col = None
for c in ["cell", "Cell", "cell_id", "Barcode", "CellID"]:
    if c in meta.columns:
        meta_cell_col = c
        break
if meta_cell_col is None:
    # if metadata index holds cell IDs
    meta.index.name = meta.index.name or "cell"
    meta.reset_index(inplace=True)
    meta_cell_col = "cell"

cells = list(set(expr.columns).intersection(set(meta[meta_cell_col].astype(str))))
expr = expr.loc[:, cells]
meta = meta[meta[meta_cell_col].astype(str).isin(cells)].copy()
meta = meta.drop_duplicates(subset=[meta_cell_col])
meta = meta.set_index(meta_cell_col).loc[cells].reset_index()

# 3) Attach hierarchy scores if available
if hier is not None:
    # Try to find cell column
    hcell = None
    for c in ["cell","Cell","cell_id","Barcode","CellID"]:
        if c in hier.columns:
            hcell = c
            break
    if hcell is None:
        hier.index.name = hier.index.name or "cell"
        hier.reset_index(inplace=True)
        hcell = "cell"
    hier = hier.drop_duplicates(subset=[hcell])
    hier = hier.set_index(hcell).reindex(cells).reset_index()
    # Merge into meta
    meta = meta.merge(hier, left_on=meta.columns[0], right_on="cell", how="left")
else:
    print("[WARN] No hierarchy score file; state labeling will rely on metadata only.")

# 4) State labeling using thresholds if scores are present
def label_state(row):
    cols = [c for c in row.index if isinstance(c, str)]
    # try to find likely score columns
    # accept names containing 'OPC','AC','OC','stem'
    sc = {k: row[k] for k in cols if any(x in k.lower() for x in ["opc","ac","oc","stem"])}
    if not sc:
        return "Unknown"
    # example: pick the max among OPC/AC/OC (ignore stem for label)
    vals = {}
    for k,v in sc.items():
        lk = k.lower()
        if "opc" in lk:
            vals["OPC-like"] = vals.get("OPC-like", -np.inf); vals["OPC-like"] = max(vals["OPC-like"], v)
        if "ac" in lk and "fac" not in lk:
            vals["AC-like"] = vals.get("AC-like", -np.inf); vals["AC-like"] = max(vals["AC-like"], v)
        if "oc" in lk and "doc" not in lk:
            vals["OC-like"] = vals.get("OC-like", -np.inf); vals["OC-like"] = max(vals["OC-like"], v)
    if not vals:
        return "Unknown"
    return max(vals, key=vals.get)

state_col = None
if hier is not None:
    meta["StateLabel"] = meta.apply(label_state, axis=1)
    state_col = "StateLabel"
else:
    # fallback if metadata already has labels
    for c in ["State","state","cluster","CellType","celltype"]:
        if c in meta.columns:
            state_col = c
            break
    if state_col is None:
        meta["StateLabel"] = "Unknown"
        state_col = "StateLabel"

# 5) Malignant vs non-malignant if provided; otherwise infer from metadata keywords
mal_col = None
for c in ["Malignant","malignant","Tumor","is_malignant"]:
    if c in meta.columns:
        mal_col = c
        break
if mal_col is None:
    # heuristic
    meta["is_malignant"] = ~meta[state_col].str.contains(
        "microglia|oligodendro|astro", case=False, na=False
    )
    mal_col = "is_malignant"

# 6) Save per-group mean expressions (for CellChat/NicheNet)
# Define groups: sample_type (Patient/PDX/GS/DGC) x state
group_cols = []
for c in ["SampleType","sample_type","Type","Model","ModelType","source"]:
    if c in meta.columns:
        group_cols.append(c)
        break
group_cols.append(state_col)
group_df = meta[group_cols].copy()
group_df.index = cells

# compute mean by group
expr_T = expr.T  # cells x genes
expr_T.index = cells
expr_T[group_cols] = group_df
mean_expr = expr_T.groupby(group_cols).mean().T  # genes x groups
mean_expr.to_csv(OUTP/"mean_expr_by_group.tsv.gz", sep="\t", compression="gzip")
print(f"[OK] Wrote {OUTP/'mean_expr_by_group.tsv.gz'}")

# 7) Pseudobulk by sample x state (per-sample granularity if sample id exists)
sample_col = None
for c in ["Sample","sample","Patient","patient","Specimen","specimen","SampleID","Case"]:
    if c in meta.columns:
        sample_col = c
        break
if sample_col is None:
    sample_col = group_cols[0]  # fall back to sample_type

expr_T["__Sample"] = meta[sample_col].values
expr_T["__State"] = meta[state_col].values
pseudobulk = expr_T.groupby(["__Sample","__State"]).mean().T
pseudobulk.to_csv(OUTP/"pseudobulk_by_sample_state.tsv.gz", sep="\t", compression="gzip")
print(f"[OK] Wrote {OUTP/'pseudobulk_by_sample_state.tsv.gz'}")

# 8) Basic counts table
counts = meta.groupby([sample_col, state_col]).size().reset_index(name="n_cells")
counts.to_csv(RES/"tables/cell_counts_by_group.tsv", sep="\t", index=False)
print(f"[OK] Wrote {RES/'tables/cell_counts_by_group.tsv'}")

# 9) Save AnnData (optional) if anndata is available
try:
    import anndata as ad
    import scipy.sparse as sp
    X = np.asarray(expr.values, dtype=np.float32)
    X = sp.csr_matrix(X)
    adata = ad.AnnData(X=X, var=pd.DataFrame(index=expr.index), obs=meta.set_index(meta.columns[0]))
    adata.write_h5ad(str(OUTP/"k27m_anndata.h5ad"))
    print(f"[OK] Wrote {OUTP/'k27m_anndata.h5ad'}")
except Exception as e:
    print("[WARN] anndata not available or failed to write .h5ad; skipping.", e)
