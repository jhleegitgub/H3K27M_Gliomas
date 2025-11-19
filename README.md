# Controllability of malignant programs in H3K27M gliomas
Single-cell RNA-seq data from H3K27M diffuse midline gliomas (Filbin et al., Science 2018)을 사용해  
악성(malignant) 세포의 co-expression 네트워크를 만들고,  
OPC-like / AC-like / OC-like 프로그램 별 structural controllability를 분석한 레포지토리입니다.

---

## 1) What’s in this repo

```text
H3K27M_Gliomas/
├─ data/                     # 입력 데이터 (원본 RSEM, metadata, TF/gene set 등)
│  ├─ raw/                   # 원본: RSEM 행렬, metadata, hierarchy, tSNE, TF/gene set
│  └─ processed/             # 전처리 산출물: HVG expr, 네트워크, netctrl 입력/결과
│
├─ results/                  # 분석 결과
│  ├─ figures/               # 그림 (pdf/png)
│  └─ tables/                # 요약 테이블 (tsv/csv)
│
├─ scripts/                  # 분석 파이프라인
│  ├─ 01_Preprocess/         # malignant 필터링, 정규화, HVG2000, Group3(OPC/AC/OC) 라벨
│  ├─ 02_Network/            # co-expression 네트워크 생성 + controllability 계산
│  └─ 03_Analysis/           # netctrl 결과 요약, figure, TF/모듈 분석
│
└─ README.md
```

## 2) Data sources & provenance
Primary paper / dataset  
PubMed  
https://pubmed.ncbi.nlm.nih.gov/29674595/

Raw data  
Single-cell RNA-seq (RSEM expression, metadata 등)
https://singlecell.broadinstitute.org/single_cell/study/SCP147/single-cell-analysis-in-pediatric-midline-gliomas-with-histone-h3k27m-mutation#study-download

TF list & gene sets
Human TF list
https://humantfs.ccbr.utoronto.ca/download.php

Lineage / cell-cycle program genes, PDGFRA module genes
Filbin et al. (Science 2018) 본문 및 보조자료에서 제공된 gene set을 다운로드하여 사용.

## 3) Analysis overview
- 대상: Filbin et al. H3K27M single-cell 데이터 중 `Type == "Malignant"` 세포만 사용.
- 전처리: 정규화 후 malignant 전체에서 HVG 2,000개 선택,  
  AC-like / OPC-like / OC-like 점수 중 최댓값으로 Group3(OPC/AC/OC) 라벨 생성.
- 네트워크: HVG 2,000에서 Pearson 상관계수 계산, |r| ≥ 0.7인 gene pair를 edge로 하는  
  co-expression 네트워크 생성 (malignant 전체 + OPC/AC/OC 그룹별).

- controllability: directed 네트워크에 structural controllability(netctrl 아이디어, maximum matching)를 적용해  
  driver node와 nD( driver 비율 )를 계산.
- 추가 분석: Human TF list, cell-cycle / PDGFRA gene set을 이용해  
  TF 및 프로그램 모듈이 허브/driver에 얼마나 겹치는지 간단히 요약.


## 4) What this repo is used for
- malignant HVG 네트워크의 기본 구조(노드·엣지 수, degree 분포, 허브 존재)를 확인한다.
- OPC-like / AC-like / OC-like 세 프로그램의 네트워크 지표와 nD를 비교해  
  악성 프로그램 간 구조·controllability 차이를 본다.
- TF / cell-cycle / PDGFRA 모듈이 네트워크에서 어떤 위치(허브/driver 쪽 vs 주변)에 있는지 살펴본다.


## 5) Caveats
- co-expression 네트워크는 상관관계 기반이므로, 실제 GRN(인과 네트워크)이 아니다.
- structural controllability 결과(driver node, nD)는 그래프 topology에 기반한 이론적 지표이며,  
  실제 약물 타겟/치료 가능성과 1:1로 대응되지는 않는다.
- HVG 2,000 공통 세트를 사용해 OPC/AC/OC를 비교했기 때문에,  
  각 상태에서만 변동성이 큰 일부 유전자는 포함되지 않을 수 있다.
- 분석 대상은 환자 malignant 세포에 한정되며, PDX/GS/DGC 등 모델 비교는 다루지 않는다.


## 6) How to use with Cytoscape
### 6-1. Malignant 전체 네트워크
파일
- 네트워크(edge):  
  `data/processed/K27M_gene_network_edges_undirected_thr0p7.tsv`
- 노드 속성(node attributes):  
  `results/tables/K27M_netctrl_gene_summary_thr0p7.tsv`

사용 방법
1. Cytoscape → File → Import → Network from File  
   → `K27M_gene_network_edges_undirected_thr0p7.tsv` 선택  
   (source / target / weight 자동 인식).
2. Cytoscape → **File → Import → Table → File  
   → `K27M_netctrl_gene_summary_thr0p7.tsv` 선택  
   → “Key Column”에서 network node 이름(`name`)과 `gene`를 매칭.
3. Style 탭에서
   - 노드 크기: `deg_total` (degree 기반)
   - 노드 색: `is_driver` (driver vs non-driver), 또는 `is_TF` (TF vs non-TF)
   로 매핑하면, 허브/driver/TF 분포를 한눈에 볼 수 있다.

### 6-2. OPC / AC / OC 그룹별 네트워크
파일
- 네트워크(edge):
  - `data/processed/K27M_group-OPC_gene_network_edges_undirected_thr0p7.tsv`
  - `data/processed/K27M_group-AC_gene_network_edges_undirected_thr0p7.tsv`
  - `data/processed/K27M_group-OC_gene_network_edges_undirected_thr0p7.tsv`
- 노드 속성:
  - `data/processed/K27M_group-OPC_gene_degrees_thr0p7.tsv`
  - `data/processed/K27M_group-AC_gene_degrees_thr0p7.tsv`
  - `data/processed/K27M_group-OC_gene_degrees_thr0p7.tsv`

*용 방법
1. 위 malignant 때와 똑같이, 그룹별 edge 파일을  
   Network from File로 불러온다.
2. 해당 그룹의 `*_gene_degrees_thr0p7.tsv`를  
   Table → To Network**로 import 하고 `gene`를 key로 매칭한다.
3. 그룹별로
   - 노드 크기: `deg_total`
   - 노드 색: `group` 또는 (원하면) driver 여부
   를 지정하면, OPC/AC/OC 네트워크 구조 차이를 시각적으로 비교할 수 있다.

> netctrl용 `*_netctrl_edges_thr0p7_ids.txt`, `*_netctrl_nodes_thr0p7.tsv` 파일들은  
> Cytoscape 시각화를 위한 것이 아니라 controllability 계산용 내부 포맷이므로,  
> 네트워크 그림에는 사용하지 않는다.
