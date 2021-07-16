# Corestem project : CA IIT Analysis


## Aim : CA 연구자 임상 데이터를 분석을 목표로 함.

#### 1.1 load environment for analysis
1.1 script에서는 library loading, gene expression data의 loading, efficacy indices data의 loading을 포함함.

#### 2.1 cleaning process for CA data analysis
2.1 script에서는 1.1의 loading을 기반으로, gene expression 데이터를 cleaning 하는 작업을 포함함.

단순한 DEG 분석인 edgeR DEG 분석, 그리고 linear mixed model을 기반으로한 DREAM DEG 분석을 추가하였음.
이후 p-value의 분포와 volcano plot을 확인함.

마지막으로 전체 샘플에 대해서 cpm을 추출하였음.

#### 3.1 CA main analysis
3.1 script에서는 2.1의 cpm을 기반으로, efficacy index인 SARA 에 대해서 유효성 평가를 진행함.
ggpubr을 통한 각 시점별로 paired boxplot을 만들어 유효성을 비교함. 

#### 3.2 CA PCR analysis
3.2 script는 사용하지 않음.
