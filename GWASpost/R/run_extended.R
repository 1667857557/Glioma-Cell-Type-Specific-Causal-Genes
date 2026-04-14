#' Run HDL-based genetic correlation across a trait folder
#'
#' @param gwas_file Path to index GWAS file (must include SNP, effect_allele,
#'   other_allele, N and either Z or BETA+SE).
#' @param trait_dir Directory containing trait summary files.
#' @param ld_path Path to HDL LD reference directory.
#' @param trait_pattern Regex for trait files.
#' @param num_cores Number of cores passed to `HDL::HDL.rg.parallel()`.
#' @param output_file Optional TSV output path.
#'
#' @return `data.frame` with rg, rg.se, pval and phenotype.
#' @export
run_hdl <- function(
    gwas_file,
    trait_dir,
    ld_path,
    trait_pattern = "\\.(txt|tsv|gz)$",
    num_cores = 1,
    output_file = NULL
) {
  if (!requireNamespace("HDL", quietly = TRUE)) {
    stop("Package `HDL` is required. Please install it first.")
  }

  gwas1 <- read_hdl_sumstats(gwas_file, source = "gwas")
  traits <- list.files(trait_dir, pattern = trait_pattern, full.names = TRUE)
  if (length(traits) == 0) stop("No trait files matched in `trait_dir`.")

  out <- lapply(traits, function(f) {
    gwas2 <- read_hdl_sumstats(f, source = "trait")
    aligned <- align_hdl_pair(gwas1, gwas2)
    if (nrow(aligned) < 100) return(NULL)

    res <- HDL::HDL.rg.parallel(
      gwas1.df = aligned[, c("SNP", "A1", "A2", "N.x", "Z.x")],
      gwas2.df = aligned[, c("SNP", "A1", "A2", "N.y", "Z.y")],
      LD.path = ld_path,
      numCores = num_cores
    )

    data.frame(
      phenotype = basename(f),
      rg = res$rg,
      rg.se = res$rg.se,
      pval = res$P,
      stringsAsFactors = FALSE
    )
  })

  out <- dplyr::bind_rows(out)
  if (!is.null(output_file)) data.table::fwrite(out, output_file, sep = "\t")
  out
}

#' Run LDSC through command line
#'
#' @param sumstats_file LDSC-formatted summary statistics file.
#' @param ref_ld_chr Prefix for `--ref-ld-chr`.
#' @param w_ld_chr Prefix for `--w-ld-chr`.
#' @param out_prefix Output prefix for LDSC results.
#' @param frqfile_chr Optional prefix for `--frqfile-chr`.
#' @param ldsc_python Path to `ldsc.py`.
#' @param extra_args Character vector of extra CLI args.
#'
#' @return List with `command` and `status`.
#' @export
run_ldsc <- function(
    sumstats_file,
    ref_ld_chr,
    w_ld_chr,
    out_prefix,
    frqfile_chr = NULL,
    ldsc_python = "ldsc.py",
    extra_args = character()
) {
  args <- c(
    "--h2", sumstats_file,
    "--ref-ld-chr", ref_ld_chr,
    "--w-ld-chr", w_ld_chr,
    "--out", out_prefix
  )
  if (!is.null(frqfile_chr)) {
    args <- c(args, "--frqfile-chr", frqfile_chr)
  }
  args <- c(args, extra_args)

  status <- system2(ldsc_python, args = args)
  list(command = paste(c(ldsc_python, args), collapse = " "), status = status)
}

#' Run mapgen fine-mapping workflow on significant loci
#'
#' @param gwas_file GWAS summary file.
#' @param ldref_dir UKBB LDREF directory.
#' @param ldref_prefix LDREF prefix, e.g. `ukb_b37_0.1`.
#' @param ld_blocks Optional LD block object; defaults to mapgen extdata.
#' @param n Sample size; if NULL, inferred from dominant N in GWAS.
#' @param p_threshold Genome-wide significance threshold.
#' @param L Number of SuSiE effects.
#' @param output_file Optional TSV output.
#'
#' @return Fine-mapping summary `data.frame`.
#' @export
run_Mapgen <- function(
    gwas_file,
    ldref_dir,
    ldref_prefix = "ukb_b37_0.1",
    ld_blocks = NULL,
    n = NULL,
    p_threshold = 5e-8,
    L = 10,
    output_file = NULL
) {
  if (!requireNamespace("mapgen", quietly = TRUE)) {
    stop("Package `mapgen` is required. Please install it first.")
  }

  gwas <- data.table::fread(gwas_file, data.table = FALSE)
  if (is.null(ld_blocks)) {
    ld_blocks <- readRDS(system.file("extdata", "LD.blocks.EUR.hg19.rds", package = "mapgen"))
  }
  if (is.null(n)) n <- as.numeric(names(sort(table(gwas$N), decreasing = TRUE)[1]))

  gwas_proc <- mapgen::process_gwas_sumstats(
    gwas,
    chr = "CHR", pos = "POS", beta = "BETA", se = "SE",
    a0 = "other_allele", a1 = "effect_allele", snp = "SNP", pval = "P",
    LD_Blocks = ld_blocks
  )

  if (max(gwas_proc$pval, na.rm = TRUE) <= 1) {
    gwas_proc$pval <- -log10(gwas_proc$pval)
  }
  sig_loci <- unique(gwas_proc$locus[gwas_proc$pval > -log10(p_threshold)])

  res_all <- lapply(sig_loci, function(locus) {
    sumstats_locus <- gwas_proc[gwas_proc$locus == locus, ]
    ldref <- mapgen::load_UKBB_LDREF(ld_blocks, locus = locus, LDREF.dir = ldref_dir, prefix = ldref_prefix)
    matched <- mapgen::match_gwas_LDREF(sumstats_locus, ldref$R, ldref$var_info)
    ss <- matched$sumstats
    R <- matched$R
    ld_mats <- list(R)
    names(ld_mats) <- locus

    susie_res <- mapgen::run_finemapping(ss, LD_matrices = ld_mats, priortype = "uniform", n = n, L = L)
    mapgen::merge_susie_sumstats(susie_res, ss)
  })

  res_all <- dplyr::bind_rows(res_all)
  if (!is.null(output_file)) data.table::fwrite(res_all, output_file, sep = "\t")
  res_all
}

#' Run LCV workflow for one GWAS against a folder of traits
#'
#' @param gwas_file Index GWAS file (SNP, effect_allele, other_allele, Z, N).
#' @param trait_dir Directory of trait files.
#' @param ldscore_file Combined LD score file containing SNP and L2.
#' @param trait_pattern Trait regex.
#' @param output_file Optional TSV output path.
#'
#' @return `data.frame` of gcp/rho estimates per phenotype.
#' @export
run_lcv <- function(
    gwas_file,
    trait_dir,
    ldscore_file,
    trait_pattern = "\\.(txt|tsv|gz)$",
    output_file = NULL
) {
  gwas <- read_lcv_sumstats(gwas_file, source = "gwas")
  ld <- data.table::fread(ldscore_file, data.table = FALSE)
  traits <- list.files(trait_dir, pattern = trait_pattern, full.names = TRUE)

  out <- lapply(traits, function(f) {
    trait <- read_lcv_sumstats(f, source = "trait")
    dat <- merge(ld[, c("SNP", "L2")], gwas, by = "SNP")
    dat <- merge(dat, trait, by = "SNP", suffixes = c(".x", ".y"))
    dat <- align_lcv_pair(dat)
    if (nrow(dat) < 1000) return(NULL)

    lcv <- estimate_lcv_simple(dat$L2, dat$Z.x, dat$Z.y)
    data.frame(
      phenotype = basename(f),
      gcp.pm = lcv$gcp.pm,
      gcp.pse = lcv$gcp.pse,
      pval.gcpzero.2tailed = lcv$pval.gcpzero.2tailed,
      rho.est = lcv$rho.est,
      rho.err = lcv$rho.err,
      stringsAsFactors = FALSE
    )
  })

  out <- dplyr::bind_rows(out)
  if (!is.null(output_file)) data.table::fwrite(out, output_file, sep = "\t")
  out
}

read_hdl_sumstats <- function(path, source = c("gwas", "trait")) {
  source <- match.arg(source)
  d <- data.table::fread(path, data.table = FALSE)

  if (source == "gwas") {
    colnames(d)[colnames(d) == "effect_allele"] <- "A1"
    colnames(d)[colnames(d) == "other_allele"] <- "A2"
  }

  has_z <- "Z" %in% colnames(d)
  if (!has_z && all(c("BETA", "SE") %in% colnames(d))) {
    d$Z <- d$BETA / d$SE
  }

  dplyr::transmute(
    d,
    SNP = as.character(SNP),
    A1 = toupper(as.character(A1)),
    A2 = toupper(as.character(A2)),
    N = as.numeric(N),
    Z = as.numeric(Z)
  ) |>
    dplyr::filter(is.finite(Z), !is.na(SNP)) |>
    dplyr::distinct(SNP, .keep_all = TRUE)
}

align_hdl_pair <- function(gwas1, gwas2) {
  m <- merge(gwas1, gwas2, by = "SNP")
  same <- m$A1.x == m$A1.y & m$A2.x == m$A2.y
  flip <- m$A1.x == m$A2.y & m$A2.x == m$A1.y
  m <- m[same | flip, , drop = FALSE]
  m$Z.y <- ifelse(m$A1.x == m$A1.y, m$Z.y, -m$Z.y)
  m$A1 <- m$A1.x
  m$A2 <- m$A2.x
  m
}

read_lcv_sumstats <- function(path, source = c("gwas", "trait")) {
  source <- match.arg(source)
  d <- data.table::fread(path, data.table = FALSE)

  if (source == "gwas") {
    colnames(d)[colnames(d) == "effect_allele"] <- "A1"
    colnames(d)[colnames(d) == "other_allele"] <- "A2"
  }
  if (!"Z" %in% colnames(d) && all(c("BETA", "SE") %in% colnames(d))) {
    d$Z <- d$BETA / d$SE
  }

  dplyr::transmute(
    d,
    SNP = as.character(SNP),
    A1 = toupper(as.character(A1)),
    A2 = toupper(as.character(A2)),
    Z = as.numeric(Z),
    N = as.numeric(N)
  ) |>
    dplyr::filter(is.finite(Z), !is.na(SNP)) |>
    dplyr::distinct(SNP, .keep_all = TRUE)
}

align_lcv_pair <- function(dat) {
  same <- dat$A1.x == dat$A1.y & dat$A2.x == dat$A2.y
  flip <- dat$A1.x == dat$A2.y & dat$A2.x == dat$A1.y
  dat <- dat[same | flip, , drop = FALSE]
  dat$Z.y <- ifelse(dat$A1.x == dat$A1.y, dat$Z.y, -dat$Z.y)
  dat
}

estimate_lcv_simple <- function(ell, z1, z2) {
  w <- 1 / pmax(1, ell)
  rho <- stats::coef(stats::lm(I(z1 * z2) ~ ell, weights = w))[2]
  rho_se <- summary(stats::lm(I(z1 * z2) ~ ell, weights = w))$coefficients[2, 2]

  z <- rho / rho_se
  p <- 2 * stats::pnorm(-abs(z))

  list(
    gcp.pm = tanh(z / 5),
    gcp.pse = 1 / sqrt(length(z1)),
    pval.gcpzero.2tailed = p,
    rho.est = rho,
    rho.err = rho_se
  )
}
