#' Run a full colocalization workflow from GWAS and local eQTL files
#'
#' This function converts the original script-based coloc workflow
#' (`Main Analyses/4.Colocalization analyses.R`) into a reusable function.
#' It reads one GWAS summary statistics file and iterates over all eQTL
#' files in a local directory, then runs `coloc::coloc.abf()` for each file.
#'
#' @param gwas_file Path to GWAS summary statistics file.
#' @param eqtl_dir Local directory containing eQTL summary statistics files.
#' @param output_file Optional output path for a CSV with coloc summaries.
#' @param eqtl_pattern Regex used to select eQTL files.
#' @param type_eqtl Trait type for eQTL data (`"quant"` or `"cc"`).
#' @param type_gwas Trait type for GWAS data (`"quant"` or `"cc"`).
#' @param s_eqtl Case proportion when `type_eqtl = "cc"`.
#' @param s_gwas Case proportion when `type_gwas = "cc"`.
#' @param min_snps Minimum number of aligned SNPs required to run coloc.
#'
#' @return A `data.frame` with one row per eQTL file and coloc posterior results.
#' @export
run_coloc <- function(
    gwas_file,
    eqtl_dir,
    output_file = NULL,
    eqtl_pattern = "\\.(txt|tsv|csv|gz)$",
    type_eqtl = "quant",
    type_gwas = "cc",
    s_eqtl = NULL,
    s_gwas = NULL,
    min_snps = 50
) {
  check_coloc_inputs(
    gwas_file = gwas_file,
    eqtl_dir = eqtl_dir,
    type_eqtl = type_eqtl,
    type_gwas = type_gwas,
    s_eqtl = s_eqtl,
    s_gwas = s_gwas,
    min_snps = min_snps
  )

  eqtl_files <- list.files(eqtl_dir, pattern = eqtl_pattern, full.names = TRUE)
  if (length(eqtl_files) == 0) {
    stop("No eQTL files found in `eqtl_dir` matching `eqtl_pattern`.")
  }

  gwas <- read_sumstats(gwas_file, source = "gwas")
  sample_size_gwas <- max(gwas$N, na.rm = TRUE)

  results <- lapply(eqtl_files, function(eqtl_file) {
    eqtl <- read_sumstats(eqtl_file, source = "eqtl")
    sample_size_eqtl <- max(eqtl$N, na.rm = TRUE)

    region_gwas <- dplyr::filter(
      gwas,
      chr == eqtl$chr[1],
      pos >= min(eqtl$pos, na.rm = TRUE),
      pos <= max(eqtl$pos, na.rm = TRUE)
    )

    merged <- align_and_merge(eqtl, region_gwas)
    if (nrow(merged) < min_snps) {
      return(NULL)
    }

    d1 <- list(
      pvalues = merged$pval_eqtl,
      N = rep(sample_size_eqtl, nrow(merged)),
      MAF = merged$maf_eqtl,
      beta = merged$beta_eqtl,
      varbeta = merged$se_eqtl^2,
      type = type_eqtl,
      snp = merged$SNP,
      chr = merged$chr,
      pos = merged$pos
    )

    d2 <- list(
      pvalues = merged$pval_gwas,
      N = rep(sample_size_gwas, nrow(merged)),
      MAF = merged$maf_gwas,
      beta = merged$beta_gwas_aligned,
      varbeta = merged$se_gwas^2,
      type = type_gwas,
      snp = merged$SNP,
      chr = merged$chr,
      pos = merged$pos
    )

    if (type_eqtl == "cc") d1$s <- s_eqtl
    if (type_gwas == "cc") d2$s <- s_gwas

    coloc_out <- coloc::coloc.abf(d1, d2)
    summary_row <- as.data.frame(as.list(coloc_out$summary))
    prior_row <- as.data.frame(as.list(coloc_out$priors))

    cbind(
      data.frame(
        eqtl_file = basename(eqtl_file),
        n_snps = nrow(merged),
        stringsAsFactors = FALSE
      ),
      prior_row,
      summary_row
    )
  })

  results <- dplyr::bind_rows(results)

  if (!is.null(output_file)) {
    data.table::fwrite(results, output_file)
  }

  results
}

check_coloc_inputs <- function(
    gwas_file,
    eqtl_dir,
    type_eqtl,
    type_gwas,
    s_eqtl,
    s_gwas,
    min_snps
) {
  if (!file.exists(gwas_file)) stop("`gwas_file` does not exist.")
  if (!dir.exists(eqtl_dir)) stop("`eqtl_dir` does not exist.")
  if (!type_eqtl %in% c("quant", "cc")) stop("`type_eqtl` must be `quant` or `cc`.")
  if (!type_gwas %in% c("quant", "cc")) stop("`type_gwas` must be `quant` or `cc`.")
  if (type_eqtl == "cc" && is.null(s_eqtl)) stop("Set `s_eqtl` when `type_eqtl = 'cc'`.")
  if (type_gwas == "cc" && is.null(s_gwas)) stop("Set `s_gwas` when `type_gwas = 'cc'`.")
  if (!is.numeric(min_snps) || min_snps < 1) stop("`min_snps` must be a positive integer.")
}

read_sumstats <- function(path, source = c("eqtl", "gwas")) {
  source <- match.arg(source)
  dt <- data.table::fread(path, data.table = FALSE)

  if (source == "eqtl") {
    needed <- c("SNP", "CHR", "BP", "FREQ", "N", "SE", "BETA", "P", "other_allele", "effect_allele")
    missing_cols <- setdiff(needed, colnames(dt))
    if (length(missing_cols) > 0) {
      stop(sprintf("eQTL file %s missing columns: %s", path, paste(missing_cols, collapse = ", ")))
    }

    dt <- dt[, needed]
    colnames(dt) <- c("SNP", "chr", "pos", "maf", "N", "se", "beta", "pval", "oa", "ea")
  } else {
    needed <- c("SNP", "CHR", "POS", "eaf", "N", "SE", "BETA", "P", "other_allele", "effect_allele")
    missing_cols <- setdiff(needed, colnames(dt))
    if (length(missing_cols) > 0) {
      stop(sprintf("GWAS file %s missing columns: %s", path, paste(missing_cols, collapse = ", ")))
    }

    dt <- dt[, needed]
    colnames(dt) <- c("SNP", "chr", "pos", "maf", "N", "se", "beta", "pval", "oa", "ea")
  }

  dt$maf <- as.numeric(dt$maf)
  dt$maf <- ifelse(dt$maf > 0.5, 1 - dt$maf, dt$maf)

  dplyr::transmute(
    dt,
    SNP = as.character(SNP),
    chr = as.character(chr),
    pos = as.numeric(pos),
    maf = as.numeric(maf),
    N = as.numeric(N),
    se = as.numeric(se),
    beta = as.numeric(beta),
    pval = as.numeric(pval),
    ea = toupper(as.character(ea)),
    oa = toupper(as.character(oa))
  ) |>
    dplyr::filter(
      is.finite(se),
      is.finite(beta),
      is.finite(pval),
      !is.na(SNP),
      !is.na(chr),
      !is.na(pos),
      !is.na(maf)
    ) |>
    dplyr::distinct(SNP, .keep_all = TRUE)
}

align_and_merge <- function(eqtl, gwas) {
  merged <- dplyr::inner_join(
    eqtl,
    gwas,
    by = "SNP",
    suffix = c("_eqtl", "_gwas")
  )

  if (nrow(merged) == 0) return(merged)

  same <- merged$ea_eqtl == merged$ea_gwas & merged$oa_eqtl == merged$oa_gwas
  flip <- merged$ea_eqtl == merged$oa_gwas & merged$oa_eqtl == merged$ea_gwas

  merged <- merged[same | flip, , drop = FALSE]
  if (nrow(merged) == 0) return(merged)

  merged$beta_gwas_aligned <- ifelse(
    merged$ea_eqtl == merged$ea_gwas,
    merged$beta_gwas,
    -merged$beta_gwas
  )

  merged$chr <- merged$chr_eqtl
  merged$pos <- merged$pos_eqtl

  merged
}
