# GWASpost

`GWASpost` 是将当前仓库中的 GWAS-post 分析脚本函数化的 R 包。

## 已实现函数

- `run_coloc()`：对接 `Main Analyses/4.Colocalization analyses.R` 的共定位流程。
- `run_hdl()`：对接 `Main Analyses/17.HDL.R` 的 HDL 遗传相关分析流程。
- `run_ldsc()`：对接 `Main Analyses/9.S_LDSC.R` 的 LDSC 命令行执行封装。
- `run_Mapgen()`：对接 `Main Analyses/2.Mapgen.R` 的 mapgen+SuSiE 精细定位流程。
- `run_lcv()`：对接 `Main Analyses/18.LCV.R` 的 LCV 批量分析流程。

## run_coloc 输入格式

### GWAS 文件必需列

`SNP, CHR, POS, eaf, N, SE, BETA, P, other_allele, effect_allele`

### eQTL 文件必需列

`SNP, CHR, BP, FREQ, N, SE, BETA, P, other_allele, effect_allele`

## run_coloc 示例

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

## 其他函数最小示例

```r
# HDL
hdl_res <- run_hdl("data/GWAS/GBM_GLIOMA.txt.gz", "data/traits", "data/HDL_LDREF")

# LDSC
ldsc_job <- run_ldsc(
  sumstats_file = "data/GBM.sumstats.gz",
  ref_ld_chr = "data/baselineLD.",
  w_ld_chr = "data/weights.",
  out_prefix = "results/gbm_ldsc"
)

# mapgen
mapgen_res <- run_Mapgen(
  gwas_file = "data/GWAS/GBM_GLIOMA.txt.gz",
  ldref_dir = "data/fine_mapping"
)

# LCV
lcv_res <- run_lcv(
  gwas_file = "data/GWAS/GBM_GLIOMA.txt.gz",
  trait_dir = "data/traits",
  ldscore_file = "data/ldscores/all_chr.ldscore.tsv.gz"
)
```
