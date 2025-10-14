import pandas as pd
from pathlib import Path

expr_path = "K27Mproject.RSEM.vh20170621.txt"   # 당신이 쓰는 표현행렬
meta_path = "PortalK27M_Metadata.vh20180223.txt"          # 포털 메타데이터

# 1) 발현행렬에서 cell_id 순서 가져오기
expr = pd.read_table(expr_path, sep="\t", index_col=0)
cell_ids = expr.columns.astype(str)

# 2) 메타데이터 로드
meta = pd.read_table(meta_path, sep="\t")
meta.columns = [c.strip() for c in meta.columns]

# 3) 필수 컬럼 지정
id_col   = "NAME"
type_col = "Type"

if id_col not in meta.columns or type_col not in meta.columns:
    raise ValueError(f"메타데이터에 {id_col}/{type_col} 컬럼이 없습니다. 현재 컬럼: {list(meta.columns)}")

meta[id_col] = meta[id_col].astype(str)

# 4) 발현행렬에 있는 셀만 남기고, 동일 순서로 정렬
meta_a = meta.set_index(id_col).reindex(cell_ids)

# 5) (A) 기본 cell_type.txt  (= Type)
ct = meta_a[type_col].fillna("NA")
Path("cell_type.txt").write_text(
    "\n".join(f"{cid}\t{lab}" for cid, lab in zip(cell_ids, ct)), encoding="utf-8"
)
print(f"✅ cell_type.txt saved ({len(ct)} rows)")

# 6) (B) 종양 세부상태: OC-like / AC-like / OPC-like 중 최댓값
state_cols = [c for c in ["OC-like","AC-like","OPC-like"] if c in meta.columns]
if state_cols:
    scores = meta_a[state_cols].apply(pd.to_numeric, errors="coerce")
    best_state = scores.idxmax(axis=1)  # 최댓값의 컬럼명
    # 종양이 아닌 세포는 NA로 두고 싶으면, Type이 종양인 경우에만 남기세요.
    # 예) malignant로 표기된 경우만:
    # mask_tumor = meta_a[type_col].str.lower().str.contains("malig|tumor|neoplastic", na=False)
    # best_state = best_state.where(mask_tumor, other="NA")

    Path("tumor_state.txt").write_text(
        "\n".join(f"{cid}\t{lab if isinstance(lab,str) else 'NA'}"
                  for cid, lab in zip(cell_ids, best_state)), encoding="utf-8"
    )
    print(f"✅ tumor_state.txt saved ({len(best_state)} rows)")
else:
    print("ℹ️ 종양 세부상태 컬럼(OC-like/AC-like/OPC-like)을 찾지 못해 tumor_state.txt는 생성하지 않았습니다.")
