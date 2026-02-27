# =============================================================================
# MPS_allmodule.R
# -----------------------------------------------------------------------------
# iCAMP: Inferring Clinical Associations of Module Perturbations
#
# Description:
#   This Shiny application computes the Module Perturbation Score (MPS) for a
#   set of genes (a "module") across TCGA cancer cohorts and visualises the
#   resulting clinical associations:
#     - Kaplan-Meier overall survival (OS) and progression-free interval (PFI)
#     - Histopathological category distributions
#     - Multivariate Cox proportional-hazards forest plots
#     - Subset survival curves stratified by clinical variables
#     - GO Biological Process enrichment for module genes
#
#   Two module modes are supported:
#     1. New module  – user uploads a CSV containing gene symbols.
#     2. Existing module – user selects from pre-built ENCODE / MSigDB /
#        LINCS / FIRE / TEISER modules stored on Cloudflare R2.
#
# Environment:
#   conda activate MPS_html
#
# Authors / Lab: Tavazoie Lab (Rockefeller University)
# =============================================================================


# =============================================================================
# 1. LOAD LIBRARIES
# =============================================================================

# --- Core survival & visualisation ---
library('survminer')    # ggsurvplot, ggforest
library('survival')     # Surv, survfit, survdiff, coxph
library('shiny')        # Shiny application framework
library('ggplot2')      # Grammar-of-graphics plotting
library('gridExtra')    # grid.arrange for multi-panel layouts
library('enrichplot')   # Enrichment result dotplot
library('grid')         # Low-level grid graphics (nullGrob)
library('clusterProfiler') # enrichGO – Gene Ontology over-representation

# --- HTML output & report generation ---
library('htmltools')    # HTML tag helpers for Shiny renderText / renderUI
library('htmlwidgets')  # Widget wrappers (used indirectly by enrichplot)
library('rmarkdown')    # render() for downloadable HTML reports

# --- Parallel computation ---
library('foreach')      # foreach() %dopar% loop syntax
library('doMC')         # Multicore backend for foreach (registerDoMC)

# --- Human gene annotation database ---
organism <- "org.Hs.eg.db"   # Bioconductor annotation package identifier
library(org.Hs.eg.db)        # Required by enrichGO for Homo sapiens

# --- Cloud storage & connectivity ---
library('rsconnect')    # Deploy app to shinyapps.io / Posit Connect
library('aws.s3')  # AWS client for Cloudflare R2
library('readr')        # Fast CSV reader (read_csv)
library('httr')         # HTTP requests (used internally by googledrive)
library('bslib')        # Bootstrap themes for Shiny UI


# =============================================================================
# 2. CLOUDFLARE R2 INFO
# =============================================================================
R2_BUCKET <- "mps-shinyapp"
R2_BASE_URL <- "f7a4dca465e33febe11c01f54163f284.r2.cloudflarestorage.com"

# NOTE: encoded variables in Posit Cloud Connect:
# "AWS_ACCESS_KEY_ID"
# "AWS_SECRET_ACCESS_KEY"

# =============================================================================
# 3. GLOBAL LOOKUP DICTIONARIES & CONSTANTS
# =============================================================================

# Mapping from human-readable TCGA disease name to short TCGA abbreviation.
# Used to construct file names for expression / clinical data on Cloudflare R2.
tcga_key <- list(
  'AML'                  = 'aml',
  'bladder'              = 'blca',
  'breast'               = 'brca',
  'cervical'             = 'cesc',
  'colon'                = 'coad',
  'esophageal'           = 'esca',
  'GBM'                  = 'gbm',
  'glioma'               = 'glm',
  'head_neck'            = 'hnsc',
  'kidney_cc'            = 'kicc',
  'kidney_ch'            = 'kich',
  'kidney_pa'            = 'kipa',
  'liver'                = 'lihc',
  'lung_ad'              = 'luad',
  'lung_sq'              = 'lusq',
  'melanoma'             = 'skcm',
  'ovarian'              = 'ovsc',
  'pancreatic'           = 'paad',
  'paraganglioma'        = 'pcpg',
  'prostate'             = 'prad',
  'rectal'               = 'read',
  'sarcoma'              = 'sarc',
  'stomach'              = 'stad',
  'testicular'           = 'tgct',
  'thymoma'              = 'thym',
  'thyroid'              = 'thca',
  'uterine_endometrial'  = 'ucec'
)

# Mapping from internal module-type code to human-readable display label.
# Used to populate the "Module Category" dropdown in the UI.
module_dict <- list(
  'CHIP'              = 'ENCODE ChIP',
  'ECLIP'             = 'ENCODE eCLIP',
  'ECLIPorithresh'    = 'ENCODE eCLIP original',
  'ECLIPthresh'       = 'ENCODE eCLIP (thresholds)',
  'FIRE1KB3p'         = 'FIRE linear RNA motifs',
  'FIRE1KBup'         = 'FIRE DNA motifs',
  'TEISERp'           = 'TEISER structural RNA motifs',
  'MSigDBGO'          = 'MSigDB Gene Ontology',
  'MSigDBHPO'         = 'MSigDB Human Phenotype Ontology',
  'MSigDBMIR'         = 'MSigDB MIRNA',
  'MSigDBONC'         = 'MSigDB cancer-driver signatures',
  'MSigDBPATH'        = 'MSigDB Pathways',
  'MSigDBTF'          = 'TF target',
  'PERTURB'           = 'DSigDB Drug signatures',
  'RBPTGTCR22thresh'  = 'ENCODE RBP CRISPR 2022 (thresholds)',
  'RBPTGTSH17thresh'  = 'ENCODE RBP SHRNA 2017 (thresholds)',
  'RBPTGTSH17'        = 'ENCODE RBP shRNA 2017',
  'RBPTGTSH22'        = 'ENCODE RBP shRNA 2022',
  'RBPTGTCR22'        = 'ENCODE RBP CRISPR 2022',
  'RBPTGTSH22thresh'  = 'ENCODE RBP SHRNA 2022 (thresholds)',
  'LINCSXPR'          = 'LINCS CRISPR (thresholds)',
  'LINCSCP'           = 'LINCS DRUG (thresholds)',
  'LINCSOE'           = 'LINCS Over Expression (thresholds)',
  'LINCSSH'           = 'LINCS shRNA knock-down (thresholds)'
)

# Reverse lookup: display label -> internal code. Used when the user selects
# a module category by its label and the code must be recovered.
module_dict_reverse <- setNames(names(module_dict), module_dict)

# Colour palette for MPS+ (red) and MPS- (blue) Kaplan-Meier curves.
col_surv <- rev(c(rgb(0.86, 0.2, 0.3, 0.75), rgb(0, 0.5, 1, 0.75)))

# Global thresholds used consistently across MI calculations.
quant_th <- 0.1   # Expression quantile filter applied when loading RDS data
MI_bins  <- 10    # Number of expression bins for the MI calculation


# =============================================================================
# 4. Cloudflare R2 I/O HELPER FUNCTIONS
# =============================================================================

# File names for shared reference datasets on Cloudflare R2.
gene_info_file_name <- "gene_files/EntrezIDs_To_ApprovedSymbol_20221107.txt"
common_genes_set_file_name <- "common_TCGA_genes_rds.rds"
cmd_mkdir <- "mkdir -p DIR"   # Shell command template used to create temp dirs


#' Read a file from Cloudflare R2 into R
#'
#' Downloads the file by name, detects its type from the extension, reads it
#' with the appropriate reader, and deletes the temporary local copy.
#'
#' @param file_name Character. Name of the file as it appears on Cloudflare R2.
#' @return The file contents as an R object (data.frame, tibble, or arbitrary
#'   R object for .rds files).
#' @details Supported extensions: csv, txt (tab-delimited), rds.
read_r2_file <- function(file_path) {
  file_type  <- tools::file_ext(file_path)
  temp_file  <- tempfile(fileext = paste0(".", file_type))

  # Download from Cloudflare R2
  aws.s3::save_object(
    object = file_path,
    bucket = R2_BUCKET,
    file = temp_file,
    region = "", 
    base_url = R2_BASE_URL,
    use_https = TRUE,
    use_path_style = TRUE 
  )

  # Read the file with the reader appropriate to its type.
  data <- switch(file_type,
    "csv" = readr::read_csv(temp_file),
    "txt" = read.table(temp_file, header = TRUE, sep = "\t",
                       quote = "\"", comment.char = ''),
    "rds" = readRDS(temp_file),
    stop("Unsupported file type")
  )

  unlink(temp_file)   # Clean up
  return(data)
}


#' Retrieve and filter gene annotation data from Cloudflare R2
#'
#' Reads the Entrez-to-symbol mapping file and returns only rows where both
#' Entrez Gene ID and Approved Symbol are non-empty.
#'
#' @param file_name Character. Name of the gene-info file on Cloudflare R2.
#' @return A data.frame with columns \code{Entrez.Gene.ID} and
#'   \code{Approved.Symbol}, filtered to rows with valid identifiers.
get_gene_info_from_r2 <- function(file_name) {
  gene_info <- read_r2_file(file_name)

  # Standardise column name for downstream merging.
  gene_info$Entrez.Gene.ID <- as.character(gene_info$Entrez.ID)

  # Keep only rows that have at least one valid identifier.
  gene_info_filt <- gene_info[
    which(gene_info$Entrez.Gene.ID != '' | gene_info$Approved.Symbol != ''),
    grep('Entrez|Approved.Symbol', colnames(gene_info))
  ]
  return(gene_info_filt)
}


# =============================================================================
# 5. CORE ANALYTICAL FUNCTIONS
# =============================================================================

# -----------------------------------------------------------------------------
# 5.1  Mutual Information (MI) Calculation
# -----------------------------------------------------------------------------

#' Calculate signed Mutual Information between a gene set and the background
#'
#' Implements a two-class MI formulation comparing genes in a module (class x)
#' against all other genes (class y) across expression bins.
#'
#' @param x Numeric vector. Bin counts for module genes (length = n_bins).
#' @param y Numeric vector. Bin counts for non-module genes (length = n_bins).
#' @return A single numeric value representing the MI between module membership
#'   and expression bin assignment.
#' @details
#'   The formula used is:
#'     MI = sum_i (1/N) * [ x_i * log(N*x_i / (sum(x)*rs_i))
#'                         + y_i * log(N*y_i / (sum(y)*rs_i)) ]
#'   where rs_i = x_i + y_i (row sum for bin i) and N = total gene count.
#'   Infinite log values (arising from zero counts) are set to 0.
calculateMI_v2 <- function(x, y) {
  rs_cMatrix <- x + y          # Row sums: total genes per expression bin
  num_genes  <- sum(rs_cMatrix) # Grand total across all bins

  # Log terms for each bin; set -Inf (log of 0) to 0 to avoid NA propagation.
  log_t1 <- log((num_genes * x) / (sum(x) * rs_cMatrix))
  log_t1[is.infinite(log_t1)] <- 0

  log_t2 <- log((num_genes * y) / (sum(y) * rs_cMatrix))
  log_t2[is.infinite(log_t2)] <- 0

  # Weighted sum of log terms to obtain MI.
  val_ret <- sum((1 / num_genes) * ((x) * log_t1 + (y) * log_t2))
}


# -----------------------------------------------------------------------------
# 5.2  Survival Analysis
# -----------------------------------------------------------------------------

#' Compute Kaplan-Meier and Cox PH survival statistics for MPS groups
#'
#' Given a clinical data.frame with pre-assigned MPS+ / MPS- group labels,
#' this function fits KM and Cox models and optionally performs a permutation
#' test to estimate empirical FDR.
#'
#' @param all_clin_df  data.frame. Clinical data with columns: group, OS,
#'   OS.time, PFI, PFI.time, mod_censor_time (added internally).
#' @param surv_type    Character. Either \code{"OS"} (overall survival) or
#'   \code{"PFI"} (progression-free interval). Default: \code{"OS"}.
#' @param rand_iter    Integer. Number of permutation iterations for empirical
#'   FDR estimation. Default: 1000.
#' @param samp_name    Character. Label for the high-MPS group. Default: "MPS+".
#' @param ctrl_name    Character. Label for the low-MPS group. Default: "MPS-".
#' @return A named list with elements:
#'   \describe{
#'     \item{KM_pv}{Log-rank p-value.}
#'     \item{KM_fit}{survfit object for KM curves.}
#'     \item{clin_data}{Input clinical data (unchanged).}
#'     \item{cox_pv}{Cox PH p-value for group coefficient.}
#'     \item{hzr}{Hazard ratio from Cox model.}
#'     \item{cox_fit}{coxph model object.}
#'     \item{avg_fit}{Overall (unconditional) survfit object.}
#'     \item{rand_pv}{Vector of permutation p-values.}
#'     \item{z}{Z-score from Cox model.}
#'     \item{concordance}{Concordance index (C-statistic).}
#'     \item{out_df}{Summary data.frame of key statistics.}
#'     \item{fdr}{Empirical FDR = proportion of permutations with p < observed p.}
#'   }
#'   Returns \code{NA} if survival data are entirely missing.
getSurv_shiny <- function(all_clin_df, surv_type = 'OS', rand_iter = 1000,
                          samp_name = 'MPS+', ctrl_name = 'MPS-') {
  library('survival')
  all_clin <- all_clin_df
  s_name   <- samp_name
  c_name   <- ctrl_name

  # Map the requested survival type to the correct time/event columns.
  if (surv_type == 'OS') {
    all_clin$mod_censor_time   <- as.numeric(all_clin$OS.time)
    all_clin$death_event_binary <- as.numeric(all_clin$OS)
  }
  if (surv_type == 'PFI') {
    all_clin$mod_censor_time   <- as.numeric(all_clin$PFI.time)
    all_clin$death_event_binary <- as.numeric(all_clin$PFI)
  }

  len_s <- sum(all_clin$group == samp_name)   # Sample size of MPS+ group
  len_c <- sum(all_clin$group == ctrl_name)   # Sample size of MPS- group

  # Unconditional KM fit (ignores grouping) used for background reference.
  s_all <- survfit(Surv(as.numeric(all_clin$mod_censor_time),
                        all_clin$death_event_binary) ~ 1)
  ss    <- Surv(as.numeric(as.character(all_clin$mod_censor_time)),
                all_clin$death_event_binary)

  val_ret <- NA
  pv <- s_fit <- pv_cox <- hz_ratio <- s_coxph <- z <- c_ix <- NA

  # Proceed only if there is at least some non-missing survival data.
  if (!(sum(is.na(ss)) == length(is.na(ss)))) {

    # Stratified KM fit by MPS group.
    s_fit <- survfit(Surv(mod_censor_time, death_event_binary) ~ group,
                     data = all_clin)

    # Log-rank test; catch degenerate cases (e.g. single-group data).
    s_diff <- tryCatch({
      survdiff(Surv(mod_censor_time, death_event_binary) ~ group,
               data = all_clin)
    }, error = function(err) {
      list(chisq = NA, n = dim(all_clin)[1])
    })

    # Convert chi-squared statistic to p-value.
    pv <- ifelse(is.na(s_diff), 1,
                 round(1 - pchisq(s_diff$chisq,
                                  length(s_diff$n) - 1), 50))[[1]]

    # Cox proportional-hazards model with MPS group as predictor.
    s_coxph <- coxph(Surv(mod_censor_time, death_event_binary) ~ group,
                     data = all_clin)

    if (!is.na(s_coxph[['coefficients']])) {
      s_c      <- summary(s_coxph)
      hz_ratio <- as.numeric(data.frame(s_c$coefficients)['exp.coef.'])
      pv_cox   <- as.numeric(data.frame(s_c$coefficients)['Pr...z..'])
      z        <- as.numeric(data.frame(s_c$coefficients)['z'])
      c_ix     <- as.numeric(summary(s_coxph)$concordance['C'])
    }

    # --- Permutation test for empirical FDR ---
    # Randomly shuffle group labels (preserving group-size proportions) and
    # recompute the log-rank p-value to build a null distribution.
    pv_r <- c()
    for (r_i in seq(rand_iter)) {
      all_clin_r       <- all_clin
      prob_v           <- as.numeric(prop.table(table(all_clin$group)))
      all_clin_r$group <- sample(c(s_name, c_name),
                                 dim(all_clin)[1], replace = TRUE,
                                 prob = prob_v)
      s1_r <- tryCatch({
        survdiff(Surv(mod_censor_time, death_event_binary) ~ group,
                 data = all_clin_r)
      }, error = function(err) {
        list(chisq = NA, n = dim(all_clin_r)[1])
      })
      pv_r <- c(pv_r,
                round(1 - pchisq(s1_r$chisq, length(s1_r$n) - 1), 50))
    }

    # Empirical FDR: proportion of permutations with p-value < observed p.
    fp <- sum(pv_r < pv) / rand_iter

    # Median survival times for each group (in the original unit, typically months).
    median_surv <- summary(s_fit)$table[, 'median']
    med_ab  <- round(as.numeric(median_surv[paste('group=', s_name, sep = '')]), 1)
    med_re  <- round(as.numeric(median_surv[paste('group=', c_name, sep = '')]), 1)
    m_sub   <- paste(s_name, ' ', len_s, '(', med_ab, ' mo)',
                     '| ', c_name, ' ', len_c, '(', med_re, ' mo)', sep = '')

    # Compact summary data.frame for display in Shiny outputs.
    out_df <- data.frame(
      pval_surv    = pv,
      pval_cox     = pv_cox,
      C_index      = c_ix,
      z_cox        = z,
      HZR_cox      = hz_ratio,
      info_survival = m_sub
    )

    val_ret <- list(
      KM_pv       = pv,
      KM_fit      = s_fit,
      clin_data   = all_clin,
      cox_pv      = pv_cox,
      hzr         = hz_ratio,
      cox_fit     = s_coxph,
      avg_fit     = s_all,
      rand_pv     = pv_r,
      z           = z,
      concordance = c_ix,
      out_df      = out_df,
      fdr         = fp
    )
  }
}


# =============================================================================
# 6. PART 1 – MPS CALCULATION FOR A USER-DEFINED (NEW) MODULE
# =============================================================================

#' Compute MPS and merge with clinical data for a user-uploaded gene list
#'
#' For each selected TCGA cancer cohort the function:
#'  1. Filters the input gene list to those present in the common TCGA gene set.
#'  2. Builds a null MI distribution by Monte-Carlo sampling of random gene sets.
#'  3. Computes a signed, Z-score-normalised MI (MPS) for every tumour sample
#'     using parallel foreach loops.
#'  4. Classifies samples as MPS+ (high) or MPS- (low) based on MPS_thresh.
#'  5. Merges MPS assignments with TCGA clinical outcome data (OS, PFI).
#'
#' @param input_genes    data.frame. First column contains gene symbols.
#' @param module_name    Character. Optional display name for the module.
#'   Auto-generated from date and random seed if left empty.
#' @param collection     Character. Data collection identifier (default 'tcga').
#' @param select_cohorts Character or vector. Either \code{"full"} (all cohorts)
#'   or a character vector of cohort names from \code{names(tcga_key)}.
#' @param number_bins    Integer. Number of expression bins for MI (default 10).
#' @param pval_thresh    Numeric. MI significance threshold (upper quantile of
#'   null distribution, default 0.01).
#' @param number_rand    Integer. Monte-Carlo iterations for null MI distribution
#'   (default 10000).
#' @param MPS_thresh     Numeric. Minimum absolute MPS Z-score to be assigned
#'   MPS+ or MPS- (default 0 = no threshold).
#' @return A data.frame combining MPS group assignments with TCGA clinical
#'   variables for all selected cohorts. Includes additional columns:
#'   \code{genes_in_module}, \code{module_name}, \code{module_id}.
#'   Survival times are converted from days to months.
get_MPS_newModule <- function(input_genes, module_name = '', collection = 'tcga',
                              select_cohorts = 'breast', number_bins = 10,
                              pval_thresh = 0.01, number_rand = 10000,
                              MPS_thresh = 0) {
  # Resolve cohort list from the selection argument.
  if (select_cohorts == 'full') {
    coh_list <- names(tcga_key)
  } else {
    coh_list <- intersect(select_cohorts, names(tcga_key))
  }

  T_1 <- Sys.time()

  # Load the set of genes common across all TCGA cohorts (used as background).
  common_genes_set <- read_r2_file(common_genes_set_file_name)
  g_set  <- common_genes_set
  inp_g  <- toupper(as.character(input_genes[, 1]))  # Normalise to upper case

  # Retain only genes present in both the module and the common gene set.
  genes_in_mod <- intersect(inp_g, g_set)
  div_f        <- as.numeric(length(genes_in_mod))

  if (length(genes_in_mod) < 1) {
    stop("Module empty. Please select another module")
  }

  # Generate a unique run ID combining a random integer and today's date.
  module_id <- paste('id_', sample(1000, 1), '_',
                     gsub('-', '', Sys.Date()), sep = '')
  if (module_name == '') { module_name <- module_id }

  set.seed(108)   # Ensure reproducible null distribution

  # --- Build null MI distribution via Monte-Carlo sampling ---
  num_rand        <- number_rand
  total_genes     <- length(g_set)
  # Distribute total genes uniformly across bins (expected bin counts).
  total_genes_bin <- rmultinom(1, total_genes, rep((1 / number_bins), number_bins))
  num_in_path     <- div_f
  # Simulate random gene sets of the same size as the module across bins.
  path_genes_bin  <- rmultinom(num_rand, num_in_path,
                               rep((1 / number_bins), number_bins))
  c_vec_r  <- path_genes_bin
  nc_vec_r <- total_genes_bin[, 1] - c_vec_r   # Non-module complement

  # Compute MI for each random simulation to form the null distribution.
  v_rand   <- sapply(seq(num_rand),
                     function(i) calculateMI_v2(c_vec_r[, i], nc_vec_r[, i]))
  mu_r     <- mean(v_rand)   # Null mean
  sd_r     <- sd(v_rand)     # Null standard deviation
  # Significance threshold: upper (1 - pval_thresh) quantile of null MI.
  v_rand_q <- as.numeric(quantile(v_rand, (1 - pval_thresh)))
  len_rand <- length(v_rand)

  all_new_clin <- c()   # Accumulator for per-cohort clinical + MPS data

  # --- Per-cohort MPS computation ---
  for (dis_ in coh_list) {
    d_t <- tcga_key[[dis_]]   # Short TCGA abbreviation

    # File names for continuous (z-score) and discretised (binned) expression.
    con_g_name <- paste0("expression_files/", d_t, '_quant_', quant_th, '_primary_zscore.rds')
    dis_g_name <- paste0("expression_files/", d_t, '_quant_', quant_th,
                         '_primary_zscore_bins', MI_bins, '.rds')

    # Load expression matrices from Cloudflare R2.
    con_g <- read_r2_file(con_g_name)   # Continuous z-scored expression
    dis_g <- read_r2_file(dis_g_name)   # Binned expression (integers 1-10)

    # Set gene symbols as row names and remove the 'gene' column.
    rownames(con_g) <- as.character(con_g$gene)
    con_g <- con_g[, -grep('gene', colnames(con_g))]
    rownames(dis_g) <- as.character(dis_g$gene)
    dis_g <- dis_g[, -grep('gene', colnames(dis_g))]

    # Align genes and samples across both matrices.
    g_row <- intersect(g_set, intersect(rownames(con_g), rownames(dis_g)))
    s_set <- intersect(colnames(con_g), colnames(dis_g))
    con_g <- con_g[g_row, s_set]
    dis_g <- dis_g[g_row, s_set]

    Nrow    <- dim(con_g)[1]
    per_bin <- round(Nrow / number_bins)   # Approximate genes per bin

    t1 <- Sys.time()

    # --- Per-sample MI calculation (parallelised with foreach / doMC) ---
    out_list <- foreach(sam_ = s_set) %dopar% {
      # Determine co-expression sign: correlation of module genes vs. all genes
      # in the continuous expression space.
      vv        <- con_g[sam_]
      vv$bin    <- 0
      vv[genes_in_mod, ]$bin <- 1   # Flag module genes
      sign_fact <- sign(cor(vv, use = "complete.obs")[1, 2])  # +1 = up-regulated, -1 = down-regulated

      # Switch to binned expression for MI calculation.
      vv  <- dis_g[sam_]
      # l_: list of gene names per bin (1 through number_bins).
      l_  <- lapply(seq(number_bins),
                    function(i) rownames(vv)[vv[sam_][, 1] == i])
      # Module gene counts per bin.
      c_vec  <- unlist(lapply(lapply(l_, intersect, genes_in_mod), length))
      # Non-module gene counts per bin.
      nc_vec <- as.numeric(sapply(l_, length) - c_vec)

      mi_    <- calculateMI_v2(c_vec, nc_vec)
      mi_sign <- sign_fact * mi_   # Apply directional sign
    }

    # Aggregate per-sample MI values and Z-score normalise against null.
    mi_vec  <- unlist(out_list)
    abs_mi  <- abs(mi_vec)
    ix_0    <- which(abs_mi < v_rand_q)   # Samples below significance threshold

    # Z-score: (|MI| - mu_null) / sd_null; clamp negatives and NAs to 0.
    v_z          <- (abs_mi - mu_r) / sd_r
    v_z[v_z < 0] <- 0
    v_z[is.na(v_z)] <- 0
    v_z          <- v_z * sign(mi_vec)    # Re-apply directional sign

    tmp_    <- data.frame(SAMPLE_ID = s_set, MPS = v_z)
    rownames(tmp_) <- s_set

    t2 <- Sys.time()
    print(paste(dis_, (t2 - t1)))   # Progress log

    d_mps <- tmp_
    samps <- as.character(d_mps$SAMPLE_ID)

    # Load TCGA clinical data for this cohort.
    d_clin_name <- paste0("clinical_files/", d_t, '_clinical_primary_forMPS.rds')
    d_clin      <- read_r2_file(d_clin_name)

    # Classify samples as MPS+ (above threshold) or MPS- (below negative threshold).
    s_p <- samps[which(d_mps$MPS > ((1) * MPS_thresh))]
    s_n <- samps[which(d_mps$MPS < ((-1) * MPS_thresh))]

    tmp_df_p <- data.frame(SAMPLE_ID = s_p)
    tmp_df_p$group <- 'MPS+'
    tmp_df_p$MPS_groups <- 'MPS+'

    tmp_df_n <- data.frame(SAMPLE_ID = s_n)
    tmp_df_n$group <- 'MPS-'
    tmp_df_n$MPS_groups <- 'MPS-'

    # Merge MPS group assignments with clinical data.
    all_clin           <- merge(d_clin, rbind(tmp_df_p, tmp_df_n))
    all_clin           <- all_clin[!is.na(all_clin$group), ]
    # Convert survival times from days to months.
    all_clin$OS.time   <- all_clin$OS.time / 30
    all_clin$PFI.time  <- all_clin$PFI.time / 30

    all_new_clin <- rbind(all_new_clin, all_clin)
  }

  # Attach module metadata columns to the combined clinical data.
  all_new_clin$genes_in_module <- paste(genes_in_mod, collapse = '|')
  all_new_clin$module_name     <- module_name
  all_new_clin$module_id       <- module_id

  T_2 <- Sys.time()
  print(T_2 - T_1)   # Total elapsed time

  val_ret <- all_new_clin
}


#' Split a gene vector into equal-sized chunks for display
#'
#' @param genes      Character vector of gene symbols.
#' @param chunk_size Integer. Maximum number of genes per chunk (default 15).
#' @return A list of character vectors, each of length <= chunk_size.
split_genes <- function(genes, chunk_size = 15) {
  n <- length(genes)
  split(genes, rep(1:ceiling(n / chunk_size),
                   each = chunk_size, length.out = n))
}


# -----------------------------------------------------------------------------
# Survival Plot for New Module
# -----------------------------------------------------------------------------

#' Generate Kaplan-Meier survival plots for OS and PFI from MPS clinical data
#'
#' Iterates over cohorts in \code{clin_data} and produces one KM plot per
#' endpoint (OS and PFI) per cohort. A cohort is plotted only if both MPS+
#' and MPS- groups have at least 20 patients.
#'
#' @param clin_data data.frame. Output of \code{get_MPS_newModule} or
#'   \code{get_MPS_existingModule}; must contain columns: cohort, group,
#'   OS, OS.time, PFI, PFI.time, genes_in_mod, module_name.
#' @return A list of length 2: \code{[[1]]} OS ggplot object,
#'   \code{[[2]]} PFI ggplot object. Returns \code{nullGrob()} for cohorts
#'   with insufficient sample sizes.
plot_Surv_get_clinicalParam_newModule <- function(clin_data) {
  pl_list_ovs <- pl_list_pfs <- list()

  genes_in_module     <- unlist(strsplit(as.character(unique(clin_data$genes_in_mod)), '\\|'))
  num_genes_in_module <- length(genes_in_module)
  name_module_trim    <- as.character(unique(clin_data$module_name))
  total_patients      <- dim(clin_data)[1]
  MPSp                <- sum(clin_data$group == "MPS+")
  MPSn                <- sum(clin_data$group == "MPS-")

  for (dis_ in as.character(unique(clin_data$cohort))) {
    all_clin <- clin_data[which(clin_data$cohort == dis_), ]
    pl_o <- pl_p <- nullGrob()   # Default to blank placeholder

    # Require at least 20 patients in each group before attempting KM fit.
    if (sum(table(all_clin$group) > 20) == length(unique(all_clin$group))) {
      surv_ovs <- getSurv_shiny(all_clin, 'OS', 1)
      surv_pfs <- getSurv_shiny(all_clin, 'PFI', 1)

      fdr_ovs <- as.numeric(surv_ovs[[12]])
      fdr_pfs <- as.numeric(surv_pfs[[12]])

      # Format empirical FDR: values of 0 are displayed as "<1e-3".
      if (fdr_ovs == 0) { f_ovs <- '1e-3' } else { f_ovs <- formatC(fdr_ovs, format = 'e', digits = 1) }
      if (fdr_pfs == 0) { f_pfs <- '1e-3' } else { f_pfs <- formatC(fdr_pfs, format = 'e', digits = 1) }

      # Build plot titles with cohort, p-value, and hazard ratio.
      m_ovs <- paste(dis_, 'p=',
                     formatC(as.numeric(surv_ovs[[11]]$pval_surv), format = 'e', digits = 1),
                     ' | HR=', signif(as.numeric(surv_ovs[[11]]$HZR_cox), 2))
      m_sub_ovs <- paste('OVS | MPS+ ', sum(all_clin$group == "MPS+"),
                         '| MPS- ', sum(all_clin$group == "MPS-"),
                         " | n= ", total_patients, " | ", name_module_trim, sep = '')

      m_pfs <- paste(dis_, 'p=',
                     formatC(as.numeric(surv_pfs[[11]]$pval_surv), format = 'e', digits = 1),
                     ' | HR=', signif(as.numeric(surv_pfs[[11]]$HZR_cox), 2))
      m_sub_pfs <- paste('PFS | MPS+ ', sum(all_clin$group == "MPS+"),
                         '| MPS- ', sum(all_clin$group == "MPS-"),
                         "| n= ", total_patients, " | ", name_module_trim, sep = '')

      # Create KM plots using survminer with a horizontal median line.
      pl_ovs <- ggsurvplot(surv_ovs[[2]], surv_ovs[[3]], pval = F,
                           title = m_ovs, font.title = 12,
                           censor.shape = 124, censor.size = 2,
                           subtitle = m_sub_ovs, font.subtitle = 10,
                           surv.median.line = 'hv', palette = col_surv,
                           risk.table = F)

      pl_pfs <- ggsurvplot(surv_pfs[[2]], surv_pfs[[3]], pval = F,
                           title = m_pfs, font.title = 12,
                           censor.shape = 124, censor.size = 2,
                           subtitle = m_sub_pfs, font.subtitle = 10,
                           surv.median.line = 'hv', palette = col_surv,
                           risk.table = F)

      pl_o <- pl_ovs$plot
      pl_p <- pl_pfs$plot
    }

    # Store the most recent cohort's plots (single-cohort mode).
    pl_list_ovs <- pl_o
    pl_list_pfs <- pl_p
  }

  val_ret <- list(pl_list_ovs, pl_list_pfs)
}


# =============================================================================
# 7. PART 2 – MPS CALCULATION FOR A PRE-BUILT (EXISTING) MODULE
# =============================================================================

#' Compute MPS and merge clinical data for a pre-built module from Cloudflare R2
#'
#' Identical MPS computation pipeline as \code{get_MPS_newModule}, but the gene
#' list is retrieved from a pre-built module file on Cloudflare R2 rather than
#' uploaded by the user.  The module file is expected to be an RDS object with
#' a \code{genes} element containing a pipe-separated gene string.
#'
#' @param module_cat    Character. Internal module category code (e.g. "RBPECLIP",
#'   "RBPCRISPR", "RBPshRNA", or any other code in \code{module_dict}).
#' @param module_type   Character. Sub-type within the category (e.g. "activated",
#'   "3UTR"). Used only for RBP categories.
#' @param RBP           Character. Name of the specific module / RBP file on
#'   Cloudflare R2 (without path prefix).
#' @param log_value     Character. Log fold-change threshold used in file naming
#'   for shRNA / CRISPR modules. Default: "1".
#' @param pval          Character. P-value threshold used in file naming.
#'   Default: "0.1".
#' @param select_cohorts Character or vector. Either \code{"full"} or a subset
#'   of \code{names(tcga_key)}.
#' @param number_bins    Integer. Expression bins (default 10).
#' @param pval_thresh    Numeric. MI significance quantile threshold (default 0.01).
#' @param number_rand    Integer. Monte-Carlo iterations (default 10000).
#' @param MPS_thresh     Numeric. Minimum |MPS| for group assignment (default 0).
#' @param collection     Character. Data collection identifier (default 'tcga').
#' @return Combined data.frame of MPS group assignments and TCGA clinical data,
#'   identical structure to output of \code{get_MPS_newModule}.
get_MPS_existingModule <- function(module_cat, module_type, RBP,
                                   log_value = "1", pval = "0.1",
                                   select_cohorts = "breast", number_bins = 10,
                                   pval_thresh = 0.01, number_rand = 10000,
                                   MPS_thresh = 0, collection = 'tcga') {
  if (select_cohorts == 'full') {
    coh_list <- names(tcga_key)
  } else {
    coh_list <- intersect(select_cohorts, names(tcga_key))
  }

  common_genes_set <- read_r2_file(common_genes_set_file_name)
  g_set <- common_genes_set
  print(module_cat)

  # Build the Cloudflare R2 file name from category, RBP, and thresholds.
  if (module_cat == "RBPECLIP") {
    module_name <- paste0(RBP)
  } else if (module_cat == "RBPshRNA" || module_cat == "RBPCRISPR") {
    # shRNA / CRISPR modules encode log and p-value thresholds in the filename.
    module_name <- paste0(RBP, "_log", log_value, "_pval", pval, ".rds")
  } else {
    module_name <- paste0(RBP)
  }

  print(RBP)
  if(nchar(module_name) > 51){
    module_name = paste0(substr(module_name,1,47),".rds")
  }
  data_RDS  <- read_r2_file(paste0("Modules_SC/",module_cat,"/",module_name))
  data_gene <- data_RDS[["genes"]]          # Pipe-separated gene string
  inp_g     <- strsplit(data_gene, "\\|")[[1]]  # Split into individual symbols

  print(inp_g)
  genes_in_mod <- intersect(inp_g, g_set)
  div_f        <- as.numeric(length(genes_in_mod))
  print(genes_in_mod)

  if (length(genes_in_mod) < 1) {
    stop("Module empty. Please select another module")
  }

  # --- Build null MI distribution (identical logic to get_MPS_newModule) ---
  set.seed(108)
  num_rand        <- number_rand
  total_genes     <- length(g_set)
  total_genes_bin <- rmultinom(1, total_genes,
                               rep((1 / number_bins), number_bins))
  num_in_path     <- div_f
  path_genes_bin  <- rmultinom(num_rand, num_in_path,
                               rep((1 / number_bins), number_bins))
  c_vec_r  <- path_genes_bin
  nc_vec_r <- total_genes_bin[, 1] - c_vec_r

  v_rand   <- sapply(seq(num_rand),
                     function(i) calculateMI_v2(c_vec_r[, i], nc_vec_r[, i]))
  mu_r     <- mean(v_rand)
  sd_r     <- sd(v_rand)
  v_rand_q <- as.numeric(quantile(v_rand, (1 - pval_thresh)))
  len_rand <- length(v_rand)

  all_new_clin <- c()

  # --- Per-cohort loop (same as get_MPS_newModule) ---
  for (dis_ in coh_list) {
    d_t        <- tcga_key[[dis_]]
    con_g_name <- paste0("expression_files/", d_t, '_quant_', quant_th, '_primary_zscore.rds')
    dis_g_name <- paste0("expression_files/", d_t, '_quant_', quant_th,
                         '_primary_zscore_bins', MI_bins, '.rds')

    con_g <- read_r2_file(con_g_name)
    dis_g <- read_r2_file(dis_g_name)

    rownames(con_g) <- as.character(con_g$gene)
    con_g <- con_g[, -grep('gene', colnames(con_g))]
    rownames(dis_g) <- as.character(dis_g$gene)
    dis_g <- dis_g[, -grep('gene', colnames(dis_g))]

    g_row <- intersect(g_set, intersect(rownames(con_g), rownames(dis_g)))
    s_set <- intersect(colnames(con_g), colnames(dis_g))
    con_g <- con_g[g_row, s_set]
    dis_g <- dis_g[g_row, s_set]

    Nrow    <- dim(con_g)[1]
    per_bin <- round(Nrow / number_bins)

    t1 <- Sys.time()

    # Parallelised per-sample MI (same logic as get_MPS_newModule).
    out_list <- foreach(sam_ = s_set) %dopar% {
      vv           <- con_g[sam_]
      vv$bin       <- 0
      vv[genes_in_mod, ]$bin <- 1
      sign_fact    <- sign(cor(vv, use = "complete.obs")[1, 2])

      vv    <- dis_g[sam_]
      l_    <- lapply(seq(number_bins),
                      function(i) rownames(vv)[vv[sam_][, 1] == i])
      c_vec  <- unlist(lapply(lapply(l_, intersect, genes_in_mod), length))
      nc_vec <- as.numeric(sapply(l_, length) - c_vec)

      mi_     <- calculateMI_v2(c_vec, nc_vec)
      mi_sign <- sign_fact * mi_
    }

    mi_vec          <- unlist(out_list)
    abs_mi          <- abs(mi_vec)
    ix_0            <- which(abs_mi < v_rand_q)
    v_z             <- (abs_mi - mu_r) / sd_r
    v_z[v_z < 0]   <- 0
    v_z[is.na(v_z)] <- 0
    v_z             <- v_z * sign(mi_vec)

    tmp_         <- data.frame(SAMPLE_ID = s_set, MPS = v_z)
    rownames(tmp_) <- s_set

    t2 <- Sys.time()
    print(paste0(dis_, "get_MPS_existingModule", (t2 - t1)))

    d_mps <- tmp_
    samps <- as.character(d_mps$SAMPLE_ID)

    d_clin_name <- paste0("clinical_files/", d_t, '_clinical_primary_forMPS.rds')
    d_clin      <- read_r2_file(d_clin_name)

    s_p <- samps[which(d_mps$MPS > ((1) * MPS_thresh))]
    s_n <- samps[which(d_mps$MPS < ((-1) * MPS_thresh))]

    # Guard: stop if no patients can be assigned to either group.
    if (length(s_p) == 0 && length(s_n) == 0) {
      stop("MPS cannot be calculated")
    }

    tmp_df_p         <- data.frame(SAMPLE_ID = s_p)
    tmp_df_p$group   <- 'MPS+'
    tmp_df_p$MPS_groups <- 'MPS+'

    tmp_df_n         <- data.frame(SAMPLE_ID = s_n)
    tmp_df_n$group   <- 'MPS-'
    tmp_df_n$MPS_groups <- 'MPS-'

    all_clin          <- merge(d_clin, rbind(tmp_df_p, tmp_df_n))
    all_clin          <- all_clin[!is.na(all_clin$group), ]
    all_clin$OS.time  <- all_clin$OS.time / 30
    all_clin$PFI.time <- all_clin$PFI.time / 30

    all_new_clin <- rbind(all_new_clin, all_clin)
  }

  all_new_clin$genes_in_module <- paste(genes_in_mod, collapse = '|')
  all_new_clin$module_name     <- module_name

  val_ret <- all_new_clin
}


# =============================================================================
# 8. PART 3 – MODULE GENE ENRICHMENT & EXPRESSION PLOTS
# =============================================================================

#' GO Biological Process enrichment analysis for module genes
#'
#' Maps module gene symbols to Entrez IDs, runs \code{clusterProfiler::enrichGO}
#' (BP ontology, FDR-adjusted), and returns a dotplot plus metadata.
#'
#' @param clin_data data.frame. Clinical data output from the MPS functions;
#'   must contain columns \code{genes_in_mod}, \code{module_id},
#'   \code{module_name}, \code{cohort}.
#' @return A list of length 5:
#'   \enumerate{
#'     \item \code{pl_enrich}       – ggplot dotplot (or NULL if no enrichment).
#'     \item \code{num_genes}       – number of module genes in the common set.
#'     \item \code{name_module}     – module display name.
#'     \item \code{enr_go}          – enrichResult object from enrichGO.
#'     \item \code{genes_in_module} – character vector of module gene symbols.
#'   }
plot_moduleGenes <- function(clin_data) {
  all_clin        <- clin_data
  dis_            <- as.character(unique(all_clin$cohort))
  module_id       <- as.character(unique(all_clin$module_id))
  name_module     <- as.character(unique(all_clin$module_name))
  name_module_trim <- strtrim(name_module, 14)   # Truncate for compact display

  genes_in_module     <- unlist(strsplit(as.character(unique(all_clin$genes_in_mod)), '\\|'))
  num_genes_in_module <- length(genes_in_module)

  # Convert approved gene symbols to Entrez IDs for GO enrichment.
  gene_info_filt <- get_gene_info_from_r2(gene_info_file_name)
  ez_ids <- unique(as.character(
    merge(gene_info_filt, data.frame(Approved.Symbol = genes_in_module))$Entrez.Gene.ID
  ))

  # Return NULL plot with metadata if no Entrez IDs can be mapped.
  if (length(ez_ids) == 0) {
    return(list(NULL, num_genes_in_module, name_module, NULL, genes_in_module))
  }

  # Run GO Biological Process over-representation analysis.
  enr_go <- enrichGO(
    gene          = ez_ids,
    organism,
    pAdjustMethod = "fdr",
    ont           = 'BP',
    pvalueCutoff  = 1e-3,
    minGSSize     = 30,
    maxGSSize     = 5000
  )

  if (is.null(enr_go) || nrow(enr_go@result) == 0) {
    pl_enrich <- NULL   # No significant GO terms found
  } else {
    # Dotplot showing top 15 GO terms ordered by gene count.
    pl_enrich <- enrichplot::dotplot(enr_go, x = 'count', showCategory = 15)
  }

  val_ret <- list(pl_enrich, num_genes_in_module, name_module, enr_go, genes_in_module)
  return(val_ret)
}


#' Volcano plot of differential gene expression between MPS+ and MPS- groups
#'
#' @param gene_DE_data data.frame. Must contain columns \code{lR} (log2 fold
#'   change MPS+ vs MPS-) and \code{fdr} (FDR-adjusted p-value).
#' @return A list of length 1 containing a ggplot volcano plot object, or a
#'   \code{nullGrob()} if fewer than 25 genes are present.
plot_geneExp <- function(gene_DE_data) {
  df_g   <- gene_DE_data
  pl_volc <- nullGrob()

  if (dim(df_g)[1] > 25) {
    pl_volc <- ggplot(df_g, aes(lR, -log10(fdr))) +
      geom_point(col = rgb(0.5, 0.5, 0.5, 0.25)) +
      theme_classic() +
      geom_hline(yintercept = -log10(0.01), col = rgb(1, 0, 0, 0.2)) +
      geom_vline(xintercept = 0, col = rgb(1, 0, 0, 0.2)) +
      xlab('log2-ratio MPS+ vs. MPS-') +
      ylab('-log10(fdr)')
  }

  val_ret <- list(pl_volc)
}


# =============================================================================
# 9. PART 4 – HISTOPATHOLOGY PLOTS
# =============================================================================

#' Return the relevant histopathological variables for a given disease cohort
#'
#' Base variables (histological type, pathologic stage, age category) are
#' supplemented with disease-specific clinical variables.
#'
#' @param disease_name Character. One of the cohort names in \code{tcga_key}.
#' @return Character vector of column names to use in histopathology plots.
get_commonHistopath <- function(disease_name) {
  dis_    <- disease_name
  hist_v  <- hist_v_ori <- c('histological_type', 'coarse_pathologic_stage', 'age_category')

  # Disease-specific additional variables:
  if (dis_ == 'breast')    { hist_v <- union(c('HR_category', 'HER2_category', 'TN_category'), hist_v_ori) }
  if (dis_ == 'prostate')  { hist_v <- union(c('gleason_score_category', 'PSA_value'), hist_v_ori) }
  if (dis_ == 'AML')       { hist_v <- union(c('mol_test_status'), hist_v_ori) }
  if (dis_ == 'head_neck') { hist_v <- union(c('HPV_status'), hist_v_ori) }
  if (dis_ == 'cervical')  { hist_v <- union(c('HPV_status'), hist_v_ori) }
  if (dis_ == 'colon')     { hist_v <- union(c('MSI_status', 'tumor_side'), hist_v_ori) }

  val_ret <- hist_v
}


#' Bar / violin plots of histopathological variables stratified by MPS group
#'
#' For categorical variables: stacked proportional bar charts (MPS+ vs MPS-).
#' For continuous variables (PSA_value): violin plot with quartile lines.
#' Plots are generated only when both groups exceed 20 patients.
#'
#' @param clin_data data.frame. Clinical data with MPS group and histopath
#'   columns; must contain \code{cohort} and \code{module_id} columns.
#' @return A named list of ggplot objects (one per histopathological variable)
#'   or a single \code{nullGrob()} if sample-size criteria are not met.
plot_Histopath <- function(clin_data) {
  all_clin  <- clin_data
  dis_      <- as.character(unique(all_clin$cohort))
  path_id   <- as.character(unique(all_clin$module_id))
  hist_v    <- get_commonHistopath(dis_)
  dis_hist  <- all_clin[, c(hist_v, 'group')]
  pl_hist   <- nullGrob()   # Default placeholder

  # Only proceed if both MPS groups have > 20 patients.
  if (sum(table(dis_hist$group) > 20) == length(unique(dis_hist$group))) {
    pl_hist <- list()

    for (h_v in sort(hist_v)) {
      if (h_v != 'PSA_value') {
        # --- Categorical variable: proportional stacked bar chart ---
        tmp_h <- dis_hist[, union('group', h_v)]
        colnames(tmp_h)[which(colnames(tmp_h) == h_v)] <- 'var'

        t_ <- table(as.character(tmp_h$var))
        tmp_h$var <- as.character(tmp_h$var)

        # Append category sample counts to labels for clarity.
        if (length(t_) > 1) {
          for (i in names(t_)) {
            ix_       <- which(tmp_h$var == i)
            tmp_h$var[ix_] <- paste(i, ' (', as.numeric(t_[i]), ')', sep = '')
          }

          pl_ <- ggplot(tmp_h, aes(var, fill = group)) +
            geom_bar(stat = "count", position = "fill", width = 0.25) +
            scale_y_continuous(labels = scales::percent) +
            theme_bw() + xlab('') +
            scale_fill_manual(values = col_surv) +
            ggtitle(gsub('_', ' ', h_v)) +
            ylab('percent') +
            coord_flip()

          pl_hist[[h_v]] <- pl_
        }
      }

      if (h_v == 'PSA_value') {
        # --- Continuous variable: violin plot (log2 scale) ---
        tmp_h <- dis_hist[, union('group', h_v)]
        colnames(tmp_h)[which(colnames(tmp_h) == h_v)] <- 'var'

        pl_ <- ggplot(tmp_h, aes(group, log2(var), fill = group)) +
          theme_bw() + xlab('') +
          ylab(paste('log2(', gsub('_', ' ', h_v), ') a.u')) +
          geom_violin(alpha = 0.5, draw_quantiles = c(0.25, 0.5, 0.75)) +
          scale_fill_manual(values = col_surv) +
          ggtitle(gsub('_', ' ', h_v)) +
          coord_flip()

        pl_hist[[h_v]] <- pl_
      }
    }
  }

  val_ret <- pl_hist
}


#' Multivariate Cox PH forest plots for OS and PFI
#'
#' Fits a Cox model with MPS group plus all relevant histopathological
#' variables as co-variates.  Variables with < 30 non-missing observations or
#' only one remaining category after small-cell exclusion are dropped.
#'
#' @param clin_data data.frame. Clinical data (output of MPS functions).
#' @return A list of length 2: \code{[[1]]} OS forest plot ggplot,
#'   \code{[[2]]} PFI forest plot ggplot. Uses \code{nullGrob()} when data
#'   are insufficient or the endpoint is not applicable (e.g. PFI in AML).
plot_Histopath_multivariate <- function(clin_data) {
  all_clin          <- clin_data
  rownames(all_clin) <- all_clin$sample_analyzed
  dis_    <- as.character(unique(all_clin$cohort))
  path_id <- as.character(unique(all_clin$module_id))
  var_vec <- c('group', get_commonHistopath(dis_))

  # --- Variable selection: remove low-information or sparse covariates ---
  ix_rm     <- c()   # Row indices to exclude from Cox fit
  var_vec_ex <- c()  # Variables to exclude from the model

  for (v_v in var_vec) {
    all_clin[v_v][, 1] <- as.character(all_clin[v_v][, 1])
    ix_na <- which(is.na(all_clin[v_v][, 1]))

    # Drop variable if fewer than 30 patients have non-missing values.
    if (dim(all_clin)[1] - length(ix_na) < 30) {
      var_vec_ex <- c(var_vec_ex, v_v)
    }
    if (dim(all_clin)[1] - length(ix_na) >= 30) {
      ix_rm <- union(ix_rm, which(is.na(all_clin[v_v][, 1])))
    }

    # Remove rows belonging to small categories (n <= 10) if there are still
    # at least 2 remaining categories; otherwise drop the variable entirely.
    t_     <- table(as.character(all_clin[v_v][, 1]))
    cat_n  <- names(which(t_ <= 10))

    if ((length(t_) - length(cat_n)) > 1 & length(cat_n) >= 1) {
      for (n_ in cat_n) {
        ix_rm <- union(ix_rm, which(all_clin[v_v][, 1] == n_))
      }
    }
    if ((length(t_) - length(cat_n)) <= 1) {
      var_vec_ex <- c(var_vec_ex, v_v)
    }
  }

  var_vec  <- setdiff(var_vec, var_vec_ex)
  all_clin <- all_clin[setdiff(seq(dim(all_clin)[1]), ix_rm), ]

  # Encode all retained variables as integers for Cox model compatibility.
  all_clin_num <- all_clin
  for (v_v in var_vec) {
    all_clin_num[v_v][, 1] <- as.numeric(factor(all_clin_num[v_v][, 1]))
  }

  ss_ovs <- Surv(as.numeric(as.character(all_clin$OS.time)), all_clin$OS)
  ss_pfs <- Surv(as.numeric(as.character(all_clin$PFI.time)), all_clin$PFI)

  val_ret      <- list(NA, NA, NA)
  p_for_ovs   <- p_for_pfs <- nullGrob()

  # Fit multivariate Cox for OS if survival data are available.
  if (!(sum(is.na(ss_ovs)) == length(is.na(ss_ovs)))) {
    su_ovs    <- Surv((all_clin$OS.time), (all_clin$OS))
    s_cox_ovs <- coxph(as.formula(paste('su_ovs ~', paste(var_vec, collapse = '+'))),
                       data = all_clin_num)
    p_for_ovs <- ggforest(s_cox_ovs, all_clin_num, main = 'OS', refLabel = "reference")
  }

  # Fit multivariate Cox for PFI if survival data are available.
  if (!(sum(is.na(ss_pfs)) == length(is.na(ss_pfs)))) {
    su_pfs    <- Surv((all_clin$PFI.time), (all_clin$PFI))
    s_cox_pfs <- coxph(as.formula(paste('su_pfs ~', paste(var_vec, collapse = '+'))),
                       data = all_clin_num)
    p_for_pfs <- ggforest(s_cox_pfs, all_clin_num, main = 'PFI', refLabel = "reference")
    # AML has no PFI endpoint defined in TCGA.
    if (dis_ == 'AML') { p_for_pfs <- nullGrob() }
  }

  val_ret <- list(p_for_ovs, p_for_pfs)
}


#' Survival sub-group analysis stratified by histopathological categories
#'
#' For each clinical variable (excluding 'group') and each of its categories,
#' fits KM curves for OS and PFI within that subset of patients.  Useful for
#' assessing MPS prognostic value within specific histological sub-types.
#'
#' @param clin_data data.frame. Clinical data (output of MPS functions).
#' @return A list of length 2: \code{[[1]]} list of OS ggplots,
#'   \code{[[2]]} list of PFI ggplots. Minimum 5 patients required per group
#'   within a sub-category to produce a plot.
plot_Histopath_subset <- function(clin_data) {
  all_clin          <- clin_data
  rownames(all_clin) <- all_clin$SAMPLE_ID
  dis_         <- as.character(unique(all_clin$cohort))
  path_id      <- as.character(unique(all_clin$module_id))
  var_vec      <- c('group', get_commonHistopath(dis_))
  total_patients <- dim(clin_data)[1]
  MPSp         <- sum(clin_data$group == "MPS+")
  MPSn         <- sum(clin_data$group == "MPS-")

  # --- Same variable filtering logic as plot_Histopath_multivariate ---
  ix_rm     <- c()
  var_vec_ex <- c()

  for (v_v in var_vec) {
    all_clin[v_v][, 1] <- as.character(all_clin[v_v][, 1])
    ix_na <- which(is.na(all_clin[v_v][, 1]))

    if (dim(all_clin)[1] - length(ix_na) < 30) {
      var_vec_ex <- c(var_vec_ex, v_v)
    }
    if (dim(all_clin)[1] - length(ix_na) >= 30) {
      ix_rm <- union(ix_rm, which(is.na(all_clin[v_v][, 1])))
    }

    t_    <- table(as.character(all_clin[v_v][, 1]))
    cat_n <- names(which(t_ <= 10))

    if ((length(t_) - length(cat_n)) > 1 & length(cat_n) >= 1) {
      for (n_ in cat_n) {
        ix_rm <- union(ix_rm, which(all_clin[v_v][, 1] == n_))
      }
    }
    if ((length(t_) - length(cat_n)) <= 1) {
      var_vec_ex <- c(var_vec_ex, v_v)
    }
  }

  var_vec  <- setdiff(var_vec, var_vec_ex)
  all_clin <- all_clin[setdiff(seq(dim(all_clin)[1]), ix_rm), ]

  c_            <- 1
  pl_list_ovs   <- pl_list_pfs <- list()

  # Iterate over variables and their categories to produce subset KM plots.
  for (v_v in setdiff(var_vec, 'group')) {
    var_cat <- as.character(unique(all_clin[v_v][, 1]))

    for (v_c in var_cat) {
      all_clin_var <- all_clin[which(all_clin[v_v][, 1] == v_c), ]
      pl_o <- pl_p <- nullGrob()

      # Require at least 5 patients per MPS group within this sub-category.
      if (sum(table(all_clin_var$group) > 5) == length(unique(all_clin_var$group))) {
        surv_ovs <- getSurv_shiny(all_clin_var, 'OS', 1)
        surv_pfs <- getSurv_shiny(all_clin_var, 'PFI', 1)

        m_ovs <- paste(v_v, ' | p=',
                       formatC(as.numeric(surv_ovs[[11]]$pval_surv), format = 'e', digits = 1),
                       ' | HR=', signif(as.numeric(surv_ovs[[11]]$HZR_cox), 2))
        m_sub_ovs <- paste('OVS (', v_c, ')', "| MPS+ ", MPSp,
                           "| MPS- ", MPSn, "| n= ", total_patients, sep = '')

        m_pfs <- paste(v_v, ' | p=',
                       formatC(as.numeric(surv_pfs[[11]]$pval_surv), format = 'e', digits = 1),
                       ' | HR=', signif(as.numeric(surv_pfs[[11]]$HZR_cox), 2))
        m_sub_pfs <- paste('PFS (', v_c, ')', "| MPS+ ", MPSp,
                           "| MPS- ", MPSn, "| n= ", total_patients, sep = '')

        pl_ovs <- ggsurvplot(surv_ovs[[2]], surv_ovs[[3]], pval = F,
                             title = m_ovs, font.title = 12,
                             censor.shape = 124, censor.size = 2,
                             subtitle = m_sub_ovs, font.subtitle = 10,
                             surv.median.line = 'hv', palette = col_surv,
                             risk.table = F)

        pl_pfs <- ggsurvplot(surv_pfs[[2]], surv_pfs[[3]], pval = F,
                             title = m_pfs, font.title = 12,
                             censor.shape = 124, censor.size = 2,
                             subtitle = m_sub_pfs, font.subtitle = 10,
                             surv.median.line = 'hv', palette = col_surv,
                             risk.table = F)

        pl_o <- pl_ovs$plot
        pl_p <- pl_pfs$plot
        pl_list_ovs[[c_]] <- pl_o
        pl_list_pfs[[c_]] <- pl_p

        # AML has no meaningful PFI definition.
        if (dis_ == 'AML') { pl_list_pfs[[c_]] <- nullGrob() }
        c_ <- c_ + 1
      }
    }
  }

  val_ret <- list(pl_list_ovs, pl_list_pfs)
}


# =============================================================================
# 10. HTML REPORT GENERATION
# =============================================================================

#' Render a downloadable HTML report via R Markdown
#'
#' Saves each plot to a temporary PNG file and passes the paths as parameters
#' to the R Markdown template at the hard-coded report path.
#'
#' @param file   Character. Destination path for the rendered HTML output file.
#' @param data   List. Named plot-data list produced by the Shiny server; must
#'   contain slots: geneEnrichPlot, survPlot, histPlot, histPlot2, histPlot3,
#'   histPlot4.
#' @return Called for its side-effect (renders the report to \code{file}).
#' @note The Rmd template path is hard-coded for the Tavazoie lab HPC
#'   (\code{/lustre/fs4/home/schhabria/pipelines/R-Shiny-SC/report.Rmd}).
generatehtml <- function(file, data) {
  # Create temporary PNG files for each plot panel.
  tmp_geneEnrichPlot <- tempfile(fileext = ".png")
  tmp_survPlot       <- tempfile(fileext = ".png")
  tmp_histPlot       <- tempfile(fileext = ".png")
  tmp_histPlot2      <- tempfile(fileext = ".png")
  tmp_histPlot3      <- tempfile(fileext = ".png")
  tmp_histPlot4      <- tempfile(fileext = ".png")

  pl <- data

  # Save gene enrichment plot (single ggplot panel).
  if (!is.null(pl$geneEnrichPlot)) {
    ggsave(tmp_geneEnrichPlot, plot = pl$geneEnrichPlot, width = 7, height = 7)
  } else {
    tmp_geneEnrichPlot <- NULL
    warning("Gene enrichment plot object is NULL")
  }

  # Save survival plot (two panels: OS and PFI side by side).
  if (!is.null(pl$survPlot)) {
    ggsave(tmp_survPlot,
           plot = grid.arrange(grobs = pl$survPlot, ncol = 2),
           width = 7, height = 7)
  } else {
    warning("Survival plot object is NULL")
  }

  # Save histopathology bar plots (auto-detect grid layout from n panels).
  if (!is.null(pl$histPlot)) {
    ggsave(tmp_histPlot,
           plot = grid.arrange(grobs = pl$histPlot,
                               ncol = ceiling(sqrt(length(pl$histPlot)))),
           width = 7, height = 7)
  } else {
    warning("Histopathology plot object is NULL")
  }

  if (!is.null(pl$histPlot2) && length(pl$histPlot2) > 0) {
    n_col <- ceiling(sqrt(length(pl$histPlot2)))
    ggsave(tmp_histPlot2,
           plot = grid.arrange(grobs = pl$histPlot2, ncol = n_col),
           height = 8, width = 9)
  } else {
    warning("Histopathology plot 2 object is NULL or empty")
  }

  if (!is.null(pl$histPlot3) && length(pl$histPlot3) > 0) {
    n_col <- ceiling(sqrt(length(pl$histPlot3)))
    ggsave(tmp_histPlot3,
           plot = grid.arrange(grobs = pl$histPlot3, ncol = n_col),
           height = 8, width = 9)
  } else {
    warning("Histopathology plot 3 object is NULL or empty")
  }

  if (!is.null(pl$histPlot4) && length(pl$histPlot4) > 0) {
    n_col <- ceiling(sqrt(length(pl$histPlot4)))
    ggsave(tmp_histPlot4,
           plot = grid.arrange(grobs = pl$histPlot4, ncol = n_col),
           height = 8, width = 9)
  } else {
    warning("Histopathology plot 4 object is NULL or empty")
  }

  # Render the R Markdown template with the saved PNG paths as parameters.
  rmarkdown::render(
    "/lustre/fs4/home/schhabria/pipelines/R-Shiny-SC/report.Rmd",
    output_file = file,
    params = list(
      gene_enrich_plot = tmp_geneEnrichPlot,
      surv_plot        = tmp_survPlot,
      hist_plot        = tmp_histPlot,
      hist_plot2       = tmp_histPlot2,
      hist_plot3       = tmp_histPlot3,
      hist_plot4       = tmp_histPlot4
    )
  )
}


# =============================================================================
# 11. SHINY APPLICATION
# =============================================================================

# -----------------------------------------------------------------------------
# 11.1  User Interface
# -----------------------------------------------------------------------------

ui <- fluidPage(
  title = 'icamp (user)',
  titlePanel("Inferring Clinical Associations of Module Perturbations (user)"),

  sidebarLayout(
    sidebarPanel(

      # Disease / cohort selector (single selection).
      selectInput("disease_name", "Select Patient Cohort",
                  choices = c("", "pancreatic", "AML", "bladder", "breast",
                              "cervical", "colon", "esophageal", "GBM",
                              "glioma", "head_neck", "kidney", "liver",
                              "lung_ad", "lung_sq", "melanoma", "ovarian",
                              "paraganglioma", "prostate", "rectal", "sarcoma",
                              "stomach", "testicular", "thymoma", "thyroid",
                              "uterine_endometrial")),

      # Toggle between uploading a new gene list or selecting a pre-built module.
      selectInput("specify_module", "Explore Modules",
                  choices = c("", existingModule = "existing", newModule = "new")),

      # --- NEW MODULE: file upload panel ---
      conditionalPanel(
        condition = "input.specify_module == 'new'",
        fileInput("module_inp", "Choose CSV File",
                  accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
      ),

      # --- EXISTING MODULE: category selection panel ---
      conditionalPanel(
        condition = "input.specify_module == 'existing'",
        selectInput("module_category", "Subject Module Categories",
                    choices = c("", sort(as.character(unlist(module_dict)))))
      ),

      # CRISPR / shRNA sub-panel: type, RBP, log threshold, p-value threshold.
      conditionalPanel(
        condition = "input.specify_module == 'existing' &&
                     (input.module_category == 'ENCODE RBP CRISPR 2023' || input.module_category == 'ENCODE RBP shRNA 2023')",
        selectInput("module_type_inp", "Type of Module",
                    choices = c("", "activated", "repressed", "differential")),
        selectInput("CRISPR_shRNA_RBP_inp", "RBP of Interest", choices = NULL),
        selectInput(inputId = 'log_value_inp',
                    label   = 'Log Value for the Module',
                    choices = c("", "0", "0.5", "1")),
        selectInput(inputId = 'pval_inp',
                    label   = 'P-value for the Module',
                    choices = c("", "0.1", "0.01", "0.05", "0.001"))
      ),

      # eCLIP sub-panel: binding region type and RBP selector.
      conditionalPanel(
        condition = "input.specify_module == 'existing' && input.module_category == 'ENCODE RBP eCLIP 2023'",
        selectInput("module_type_inp", "Type of Module",
                    choices = c("", "3UTR", "5UTR", "exon", "intron", "intergenic")),
        selectInput("eCLIP_RBP_inp", "RBP of Interest", choices = NULL)
      ),

      # Generic module selector for all non-RBP categories.
      conditionalPanel(
        condition = "input.specify_module == 'existing' &&
                     input.module_category != 'ENCODE RBP CRISPR 2023' &&
                     input.module_category != 'ENCODE RBP shRNA 2023' &&
                     input.module_category != 'ENCODE RBP eCLIP 2023'",
        selectInput("Gene_type_inp", "Module of Interest", choices = NULL)
      ),

      # Buttons: load available modules for the selected category, then submit.
      actionButton("load_rbp_values", "Load Module Categories"),
      actionButton("submit", "Submit")
      #downloadButton("downloadhtml", "Download html Report")  # disable for now
    ),

    mainPanel(
      uiOutput("message"),
      tabsetPanel(type = "tabs",

        # Dashboard: welcome text and terminology guide.
        tabPanel("Dashboard", htmlOutput("readme")),

        # Module Info: gene list, count, and GO enrichment dotplot.
        tabPanel("Module Info (Genes)",
                 htmlOutput('module_info'),
                 verbatimTextOutput('volc_info')),

        # Clinical Data: summary stats and downloadable patient table.
        tabPanel("Clinical Data",
                 htmlOutput('table_info'),
                 htmlOutput('clinical_info'),
                 downloadButton("downloadData",
                                "Download clinical information of the patients")),

        # Survival: KM curves for OS and PFI.
        tabPanel("Survival",
                 div(class = "plot-container",
                     downloadButton("download_survPlot", "Download Survival Plot"),
                     plotOutput("survPlot", height = "400px")),
                 tags$style(type = "text/css",
                            ".plot-container { position: relative; width = 100%; height: 400px; text-align: center; }
                             #download_survPlot { margin-button: 10px; }")),

        # Histopathology: proportional bar / violin plots per clinical variable.
        tabPanel("histopathology",
                 div(class = "plot-container",
                     downloadButton("download_histPlot", "Download Histopathology plot"),
                     plotOutput("histPlot", height = "400px")),
                 tags$style(type = "text/css",
                            ".plot-container { position: relative; width: 100%; text-align: center; }
                             #download_histPlot { margin-bottom: 10px; }")),

        # Histopathology (multivariate): forest plots from Cox PH models.
        tabPanel("histopathology (multivariate)",
                 div(class = "plot-container",
                     downloadButton("download_histPlot2", "Download histopathology (multivariate) plot"),
                     plotOutput("histPlot2", height = "400px")),
                 tags$style(type = "text/css",
                            ".plot-container { position: relative; width: 100%; text-align: center; }
                             #download_histPlot2 { margin-bottom: 10px; }")),

        # Histopathology (subset: OVS): KM curves within clinical sub-groups (OS endpoint).
        tabPanel("histopathology (subset: OVS)",
                 div(class = "plot-container",
                     downloadButton("download_histPlot3", "Download histopathology (subset: OVS) plot"),
                     plotOutput("histPlot3", height = "400px")),
                 tags$style(type = "text/css",
                            ".plot-container { position: relative; width: 100%; text-align: center; }
                             #download_histPlot3 { margin-bottom: 10px; }")),

        # Histopathology (subset: PFS): KM curves within clinical sub-groups (PFI endpoint).
        tabPanel("histopathology (subset: PFS)",
                 div(class = "plot-container",
                     downloadButton("download_histPlot4", "Download histopathology (subset: OVS) plot"),
                     plotOutput("histPlot4", height = "400px")),
                 tags$style(type = "text/css",
                            ".plot-container { position: relative; width: 100%; text-align: center; }
                             #download_histPlot4 { margin-bottom: 10px; }"))
      )
    )
  )
)


# -----------------------------------------------------------------------------
# 11.2  Server Logic
# -----------------------------------------------------------------------------

server <- function(input, output, session) {

  # --- Reactive state ---
  plot_data     <- reactiveVal(NULL)   # Holds all plot objects for current run
  data          <- reactiveVal(NULL)   # Holds the full clinical data.frame
  error_message <- reactiveVal(NULL)   # Holds error text for modal display

  # Refresh module category dropdown when module mode changes.
  observeEvent(input$specify_module, {
    updateSelectInput(session, "module_category",
                      choices = c("", sort(as.character(unlist(module_dict)))))
  })


  # --- Load module file list from Cloudflare R2 when "Load Module Categories" is clicked ---
  observeEvent(input$load_rbp_values, {
    req(input$module_category)
    internal_value <- module_dict_reverse[[input$module_category]]
    parent_dir     <- "Modules_SC"
    file_dir       <- file.path(parent_dir, internal_value)
    rbp_values     <- NULL

    if (internal_value %in% c('RBPECLIP', 'RBPCRISPR', 'RBPshRNA')) {
      # For RBP categories, list files in the sub-type directory.
      req(input$module_type_inp)
      file_dir   <- file.path(parent_dir, internal_value, isolate(input$module_type_inp))
      files      <- list_files_from_r2(file_dir)
      rbp_values <- files
    } else {
      # For non-RBP categories, load a pre-built module-name index RDS file.
      module_name <- paste0(internal_value, ".rds")
      module_list <- read_r2_file(paste0("Modules_SC_list/", module_name))
      rbp_values  <- module_list   # Vector of available module names
    }

    # Populate the appropriate dropdown with the available module names.
    if (internal_value %in% c('RBPCRISPR', 'RBPshRNA')) {
      updateSelectInput(session, "CRISPR_shRNA_RBP_inp", choices = rbp_values)
    } else if (internal_value == 'RBPECLIP') {
      updateSelectInput(session, "eCLIP_RBP_inp", choices = rbp_values)
    } else {
      updateSelectInput(session, "Gene_type_inp", choices = rbp_values)
    }
  })


  # --- Submit handler: NEW MODULE ---
  observeEvent(input$submit, {
    req(input$specify_module == 'new')
    showModal(modalDialog("Processing data. Please wait...", footer = NULL))

    clin_data_val <- NULL
    req(input$module_inp)

    inFile       <- input$module_inp
    input_genes  <- read.csv(inFile$datapath, sep = ',', header = F)
    disease_name <- isolate(input$disease_name)

    # Compute MPS for the uploaded gene list.
    clin_data_val <- isolate(
      get_MPS_newModule(input_genes, select_cohorts = disease_name)
    )

    if (!is.null(clin_data_val)) {
      # Generate all plot types from the computed clinical data.
      pl        <- plot_moduleGenes(clin_data_val)
      survPlot  <- plot_Surv_get_clinicalParam_newModule(clin_data_val)
      histPlot  <- plot_Histopath(clin_data_val)
      histPlot2 <- plot_Histopath_multivariate(clin_data_val)
      histPlot3 <- plot_Histopath_subset(clin_data_val)

      # Bundle all plots into a named list for reactive storage.
      plot_data_list <- list(
        geneEnrichPlot = pl[[1]],    # GO enrichment dotplot
        moduleInfo     = pl[[3]],    # enrichResult object
        geneCount      = pl[[2]],    # Number of genes in module
        genes          = pl[[5]],    # Gene symbol vector
        survPlot       = survPlot,   # list(OS_plot, PFI_plot)
        histPlot       = histPlot,   # Named list of bar/violin plots
        histPlot2      = histPlot2,  # list(OS_forest, PFI_forest)
        histPlot3      = histPlot3[[1]],  # OS subset KM plots
        histPlot4      = histPlot3[[2]]   # PFI subset KM plots
      )

      plot_data(plot_data_list)
      data(clin_data_val)
      removeModal()

    } else {
      # Show error modal if no clinical data was returned.
      error_message("An error occurred while processing the new module.")
      showModal(modalDialog(
        title     = "Error",
        div(style = "color: red; font-weight: bold;", error_message()),
        easyClose = FALSE,
        footer    = modalButton("Close")
      ))
    }
  })


  # --- Submit handler: EXISTING MODULE ---
  observeEvent(input$submit, {
    req(input$specify_module == 'existing')
    showModal(modalDialog("Processing data. Please wait...", footer = NULL))

    clin_data_val <- NULL
    # At least one module selector must be non-NULL before proceeding.
    req(!is.null(input$CRISPR_shRNA_RBP_inp) ||
          !is.null(input$eCLIP_RBP_inp) ||
          !is.null(input$Gene_type_inp))

    internal_value <- module_dict_reverse[[input$module_category]]
    disease_name   <- input$disease_name

    tryCatch({
      # Route to the correct get_MPS_existingModule() call based on category.
      if (internal_value == 'RBPECLIP') {
        clin_data_val <- get_MPS_existingModule(
          module_cat  = internal_value,
          module_type = input$module_type_inp,
          RBP         = input$eCLIP_RBP_inp,
          log_value   = input$log_value_inp,
          pval        = input$pval_inp,
          select_cohorts = disease_name
        )
      } else if (internal_value == "RBPCRISPR" || internal_value == "RBPshRNA") {
        clin_data_val <- get_MPS_existingModule(
          module_cat  = internal_value,
          module_type = input$module_type_inp,
          RBP         = input$CRISPR_shRNA_RBP_inp,
          log_value   = input$log_value_inp,
          pval        = input$pval_inp,
          select_cohorts = disease_name
        )
      } else {
        # All other module types (MSigDB, FIRE, TEISER, LINCS, …).
        clin_data_val <- get_MPS_existingModule(
          module_cat  = internal_value,
          module_type = input$module_type_inp,
          RBP         = input$Gene_type_inp,
          select_cohorts = disease_name
        )
      }

      if (!is.null(clin_data_val)) {
        pl        <- plot_moduleGenes(clin_data_val)
        survPlot  <- plot_Surv_get_clinicalParam_newModule(clin_data_val)
        histPlot  <- plot_Histopath(clin_data_val)
        histPlot2 <- plot_Histopath_multivariate(clin_data_val)
        histPlot3 <- plot_Histopath_subset(clin_data_val)

        plot_data_list <- list(
          geneEnrichPlot = pl[[1]],
          moduleInfo     = pl[[3]],
          geneCount      = pl[[2]],
          genes          = pl[[5]],
          survPlot       = survPlot,
          histPlot       = histPlot,
          histPlot2      = histPlot2,
          histPlot3      = histPlot3[[1]],
          histPlot4      = histPlot3[[2]]
        )

        plot_data(plot_data_list)
        data(clin_data_val)
        removeModal()
      }

    }, error = function(e) {
      cat("Error in processing existing module:", e$message, "\n")
      error_message(e$message)
      showModal(modalDialog(
        title     = "Error",
        div(style = "color: red; font-weight: bold;", error_message()),
        easyClose = FALSE,
        footer    = modalButton("Close")
      ))
    })
  })


  # --- Error message output (plain text fallback) ---
  output$error_message <- renderText({
    error_message()
  })


  # --- Dashboard tab: welcome / terminology HTML ---
  output$readme <- renderUI({
    readme_text <- "
    <h2><strong>Welcome to Your Shiny Dashboard!</strong></h2>

    <p>This tool is crafted to facilitate your exploration of the expression patterns of biologically coherent genes defined as modules within the TCGA cohort.
    It analyzes the coordinated shifts in mRNA expression among the genes listed in the module and conducts survival analysis based on these findings.
    Here's a quick guide to help you get started::</p>

    <h3 style=\"color: blue;\">  Terms used in this dashboard  </h3>
    <ul>
        <li><strong>Modules:</strong> List of biologically coherant genes like pathway.</li>
        <li><strong>Module Pertubation Score (MPS):</strong> Coordinated shift in the mRNA expression of a set of genes in a patient's tumor transcriptome.</li>
        <ul>
                    <li> MPS+ : up-regulation of genes in the module in the patient</li>
                    <li> MPS- : down-regulation of genes in the module in the patient</li>
        </ul>
    </ul>

    <h3 style=\"color: blue;\">  Module Definitions  </h3>
    <ul>
        <li><h4><strong>New module:</strong></h4>
           Upload the CSV file containing a list of your genes of interest. The CSV should be formatted as a single column, each gene in its own row, with no header.</li>
        <li><h4><strong>Existing modules:</strong></h4>
            Precomputed genesets from a variety of sources, including:</li>
            <ul>
                  <li> ENCODE ChIP, eCLIP (with and without thresholding), and RBP screen data (with and without thresholding).</li>
                  <li> NIH LINCS CRISPR, overexpression, shRNA, and drug data. <strong>Note that LINCS drug data is large and category choices will take longer to load.</strong></li>
                  <li> MSigDB GO terms, human phenotype ontologies, TF targets, pathways, cancer signatures, and miRNA signatures.</li>
                  <li> DSigDB drug signatures.</li>
                  <li> Linear RNA and DNA motifs calculated from FIRE.</li>
                  <li> Structural RNA motifs calculated from TEISER.</li>
            </ul>
        </li>
    </ul>

    <p><strong>Click on SUBMIT and switch to the other tabs. We hope you find this dashboard useful and intuitive. Happy exploring!</strong></p>

    "
    HTML(readme_text)
  })


  # --- Module Info tab: gene list and count ---
  output$module_info <- renderText({
    pl <- plot_data()
    if (is.null(pl)) return()

    # Split genes into chunks of 15 for readable line-wrapped display.
    chunked_genes <- split_genes(pl$genes)
    genes_str     <- paste(sapply(chunked_genes, paste, collapse = ", "),
                           collapse = "<br>")

    paste('<b> Name of module: </b>', pl$moduleInfo, "<br>",
          '<b> Number of genes in the module: </b>', pl$geneCount, "<br>",
          '<b> Genes in the module: </b><br>', genes_str, "<br>")
  })


  # --- Gene enrichment dotplot (currently displayed within Module Info tab) ---
  output$geneEnrichPlot <- renderPlot({
    pl <- plot_data()
    if (is.null(pl) || is.null(pl$geneEnrichPlot)) {
      plot.new()
      text(0.5, 0.5, "No significance was found.")
      return()
    }
    plot(pl$geneEnrichPlot)
  }, height = 600, width = 900)

  output$download_geneEnrichPlot <- downloadHandler(
    filename = function() { "gene_enrichment_plot.pdf" },
    content  = function(file) {
      pl <- plot_data()
      if (!is.null(pl$geneEnrichPlot)) { ggsave(file, pl$geneEnrichPlot) }
    }
  )


  # --- Clinical Data tab ---
  output$table <- renderTable({ data() })

  output$downloadData <- downloadHandler(
    filename = function() { "clinical_info.csv" },
    content  = function(file) {
      x <- data()
      write.csv(x, file, row.names = FALSE)
    }
  )

  # Summary counts of MPS+ and MPS- patients.
  output$clinical_info <- renderText({
    x <- data()
    x <- paste('<b> number of patients: ', dim(x)[1], "<br>",
               'MPS+ ', sum(x$group == 'MPS+'),
               ' | MPS- ', sum(x$group == 'MPS-'),
               "<br>", "<br>", "</b>")
  })


  # --- Survival tab ---
  output$survPlot <- renderPlot({
    pl <- plot_data()
    if (is.null(pl)) return()
    grid.arrange(pl$survPlot[[1]], pl$survPlot[[2]], ncol = 2)
  }, height = 600, width = 600)

  output$download_survPlot <- downloadHandler(
    filename = function() { "survival_plot.pdf" },
    content  = function(file) {
      pl <- plot_data()
      if (!is.null(pl$survPlot)) {
        ggsave(file,
               plot   = grid.arrange(pl$survPlot[[1]], pl$survPlot[[2]], ncol = 2),
               width  = 7, height = 7)
      }
    }
  )


  # --- Histopathology tab (univariate bar / violin plots) ---
  output$histPlot <- renderPlot({
    pl         <- plot_data()
    if (is.null(pl)) return()
    hist_plots <- pl$histPlot
    n_col      <- ceiling(sqrt(length(hist_plots)))
    grid.arrange(grobs = hist_plots, ncol = n_col)
  }, height = 800, width = 900)

  output$download_histPlot <- downloadHandler(
    filename = function() { "histopathology_plot.pdf" },
    content  = function(file) {
      pl <- plot_data()
      if (!is.null(pl$histPlot)) {
        ggsave(file,
               plot   = grid.arrange(grobs = pl$histPlot,
                                     ncol  = ceiling(sqrt(length(pl$histPlot)))),
               width  = 9, height = 8)
      }
    }
  )


  # --- Histopathology (multivariate) tab ---
  output$histPlot2 <- renderPlot({
    pl    <- plot_data()
    p_    <- pl$histPlot2
    n_col <- ceiling(sqrt(length(p_)))
    if (n_col == 0) { grid.arrange(nullGrob()) } else { grid.arrange(grobs = p_, ncol = n_col) }
  }, height = 800, width = 900)

  output$download_histPlot2 <- downloadHandler(
    filename = function() { "histopathology_plot2.pdf" },
    content  = function(file) {
      pl <- plot_data()
      if (!is.null(pl$histPlot2)) {
        ggsave(file,
               plot   = grid.arrange(grobs = pl$histPlot2, ncol = 2),
               width  = 8, height = 9)
      }
    }
  )


  # --- Histopathology (subset: OVS) tab ---
  output$histPlot3 <- renderPlot({
    pl    <- plot_data()
    if (is.null(pl)) return()
    p_    <- pl$histPlot3
    n_col <- ceiling(sqrt(length(p_)))
    if (n_col == 0) { grid.arrange(nullGrob()) } else { grid.arrange(grobs = p_, ncol = n_col) }
  }, height = 800, width = 900)

  output$download_histPlot3 <- downloadHandler(
    filename = function() { "histopathology_plot_OVS.pdf" },
    content  = function(file) {
      pl <- plot_data()
      if (!is.null(pl$histPlot3)) {
        ggsave(file,
               plot   = grid.arrange(grobs = pl$histPlot3,
                                     ncol  = ceiling(sqrt(length(pl$histPlot3)))),
               width  = 9, height = 8)
      }
    }
  )


  # --- Histopathology (subset: PFS) tab ---
  output$histPlot4 <- renderPlot({
    pl    <- plot_data()
    if (is.null(pl)) return()
    p_    <- pl$histPlot4
    n_col <- ceiling(sqrt(length(p_)))
    if (n_col == 0) { grid.arrange(nullGrob()) } else { grid.arrange(grobs = p_, ncol = n_col) }
  }, height = 800, width = 900)

  output$download_histPlot4 <- downloadHandler(
    filename = function() { "histopathology_plot_PFS.pdf" },
    content  = function(file) {
      pl <- plot_data()
      if (!is.null(pl$histPlot4)) {
        ggsave(file,
               plot   = grid.arrange(grobs = pl$histPlot4,
                                     ncol  = ceiling(sqrt(length(pl$histPlot4)))),
               width  = 9, height = 8)
      }
    }
  )


observeEvent(input$submit, {
  updateTabsetPanel(session, "tabs", "module_info")
})

output$downloadhtml <- downloadHandler(
  filename = function() { paste("MPS_Report_", Sys.Date(), ".html", sep="") },
  content = function(file) {
    tempReport <- file.path(tempdir(), "report.Rmd")
    file.copy("report.Rmd", tempReport, overwrite = TRUE)

    # Ensure plot_data() exists in your reactive environment
    params <- list(pl = plot_data())

    rmarkdown::render(tempReport, output_file = file,
                      params = params,
                      envir = new.env(parent = globalenv()))
  }
)
}


# =============================================================================
# 12. LAUNCH THE APPLICATION
# =============================================================================

shinyApp(ui = ui, server = server)
