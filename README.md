# 1) What’s in this repo
...
├─ data/
│  └─ PortalK27M_Metadata.vh20180...      # 예: 메타데이터/라벨 테이블(소용량만)
├─ script/
│  ├─ PreAnalysis.bash                    # 환경셋업·다운로드·전처리(파이프라인 진입점)
│  ├─ run_cell_type.py                    # 셀타입 라벨 처리/매핑/요약
│  └─ run_expression_network.py           # GRN 추론 + driver-node 분석(TENET/netctrl 래퍼)
└─ README.md
...
# 2) Data sources & provenance
Primary paper / dataset
PubMed: https://pubmed.ncbi.nlm.nih.gov/29674595/
GEO (GSE102130): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102130
