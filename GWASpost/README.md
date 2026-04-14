# GWASpost

`GWASpost` 是将当前仓库中的 GWAS-post 分析脚本函数化的起始 R 包。

## 已实现函数

- `run_coloc()`：对接 `Main Analyses/4.Colocalization analyses.R` 的完整共定位流程。

## run_coloc 输入格式

### GWAS 文件必需列

`SNP, CHR, POS, eaf, N, SE, BETA, P, other_allele, effect_allele`

### eQTL 文件必需列

`SNP, CHR, BP, FREQ, N, SE, BETA, P, other_allele, effect_allele`

## 示例

```r
library(GWASpost)

res <- run_coloc(
  gwas_file = "data/GWAS/GBM_GLIOMA.txt.gz",
  eqtl_dir = "data/eQTL",
  output_file = "results/coloc_summary.csv",
  type_eqtl = "quant",
  type_gwas = "cc",
  s_gwas = 12496/30686,
  min_snps = 50
)
```

