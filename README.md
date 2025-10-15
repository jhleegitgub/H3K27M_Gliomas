## 1) What’s in this repo
```text
H3K27M_Gliomas/
├── data/
│   ├── raw/                  # 원본: RSEM 행렬, metadata, hierarchy, tSNE  (git ignore)
│   └── processed/            # 전처리 산출물: AnnData, pseudobulk, mean expr
├── results/
│   ├── figures/              # 그림 (pdf/png)
│   └── tables/               # 요약 테이블 (tsv/csv)
├── scripts/
│   ├── 01_Preprocess/        # RSEM+메타데이터 통합, 라벨링
│   ├── 02_CellChat_NicheNet/ # 리간드–수용체 네트워크
│   ├── 03_SCENIC_pySCENIC/   # regulon(AUCell) 허브
│   ├── 04_CompensationNet/   # PDGFRA↔ERBB/FGFR 보상망
│   └── 05_ModulePreservation/# WGCNA 모듈 보존(환자 vs PDX/GS/DGC)
└── README.md

```
## 2) Data sources & provenance

Primary paper / dataset  
PubMed  
https://pubmed.ncbi.nlm.nih.gov/29674595/  

Raw data  
https://singlecell.broadinstitute.org/single_cell/study/SCP147/single-cell-analysis-in-pediatric-midline-gliomas-with-histone-h3k27m-mutation#study-download  
