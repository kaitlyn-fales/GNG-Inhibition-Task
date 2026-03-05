########### Preprocess data for inhibition task DCM ##################

# Packages
library(RNifti)
library(nlme)
library(dplyr)

# Directory for results
data_dir <- "Data" # to store .RData versions

################ List of subject IDs to loop through #######################
subjects <- c(110,111,120,121,124,128,134,143,152,160,171,172,173,181,
              184,196,199,214,215,223,227,230,247,252,256,258,265,266,268,
              271,275,276,277) # sub 176 weird
############################################################################

# Load mask
mask <- readNifti("mask.nii")

# Base directory for raw data
base_dir <- "/storage/group/alh98/default/VLN_BIDS/fmri.analysis.bai.alh/derivatives/derfmrinf.23_trimmed"
timing_dir <- "/storage/group/alh98/default/VLN_BIDS/participant_data_complete"

# Specify TR
TR <- 2

########## Functions ##################################################
# Creation of stimulus indicator function
stim.indicator <- function(os,dur,times){
  stim <- rep(0,length(times)-1)
  for (i in 1:length(os)){
    stim[os[i]:(os[i] + dur[i] - 1)] <- 1
  }
  return(stim)
}

# Regress out motion in each voxel
MotionRegress <- function(y, X){
  if (var(y) != 0){
    df <- data.frame(y, X)
    residuals <- as.numeric(gls(y ~ ., data = df, correlation = corAR1())$residuals)
  } else {
    residuals <- NA
  }
  return(residuals)
}

# Function to extract first PC as representative BOLD time course
get_BOLD_eigenvariate <- function(BOLD) {
  
  # Do PCA and pull out relevant quantities
  pca <- prcomp(BOLD, center = TRUE, scale. = FALSE)
  
  V <- pca$rotation[,1]   
  U <- pca$x[,1]  
  
  # Enforce sign convention: mean positive voxel loadings
  d <- sign(sum(V)); if (d == 0) d <- 1
  Y <- (U * d) / sqrt(ncol(BOLD)) 
  
  return(Y)
}

# Regress only 6 motion parameters + spike regressors + DCT for drift correction
regress_param_motion <- function(confounds) {
  # --- Sanitize columns to numeric ---
  confounds[] <- lapply(confounds, function(col) {
    if (is.factor(col)) col <- as.character(col)
    if (is.character(col)) {
      col[col %in% c("", "n/a", "NA", "None")] <- NA
      suppressWarnings(as.numeric(col))
    } else {
      as.numeric(col)
    }
  })
  
  # --- 1. Six motion parameters ---
  mot6_names <- c("trans_x", "trans_y", "trans_z",
                  "rot_x", "rot_y", "rot_z")
  mot6 <- if (all(mot6_names %in% colnames(confounds))) {
    as.matrix(confounds[, mot6_names, drop = FALSE])
  } else {
    stop("Missing one or more motion parameter columns.")
  }
  
  # --- 2. Motion outlier (spike) regressors ---
  spike_cols <- grep("^motion_outlier", colnames(confounds), value = TRUE)
  spikes <- if (length(spike_cols) > 0) {
    tmp <- as.matrix(confounds[, spike_cols, drop = FALSE])
    tmp <- tmp[, colSums(tmp, na.rm = TRUE) > 0, drop = FALSE]
    if (ncol(tmp) > 0) tmp else NULL
  } else NULL
  
  # --- 3. First three cosine (drift) regressors from fMRIPrep ---
  cos_cols_all <- grep("^cosine\\d{2}$", colnames(confounds), value = TRUE)
  cos_use <- head(cos_cols_all, 3)
  cosines <- if (length(cos_use) > 0) {
    as.matrix(confounds[, cos_use, drop = FALSE])
  } else NULL
  
  # --- 4. Combine all parts and clean ---
  parts <- list(mot6, spikes, cosines)
  parts <- parts[!vapply(parts, function(x)
    is.null(x) || nrow(x) == 0 || ncol(x) == 0, logical(1))]
  
  if (length(parts) == 0) stop("No valid confound columns found.")
  
  X <- do.call(cbind, parts)
  X[is.na(X)] <- 0
  
  # Drop constant columns (just in case)
  keep <- apply(X, 2, function(col) sd(col, na.rm = TRUE) > 1e-8)
  X <- X[, keep, drop = FALSE]
  
  return(as.matrix(X))
}
########################################################################

########## Loop to process subject data ################################
# Process subjects and runs with 2-column u (general task + inhibition)
for (k in 1:length(subjects)) {
  sub <- subjects[k]
  cat("\nProcessing subject:", sub, "\n")
  
  # Detect available runs
  func_dir <- file.path(base_dir, paste0("sub-", sub), "ses-1", "func")
  bold_files <- list.files(func_dir, pattern = "_masked\\.nii\\.gz$", full.names = TRUE)
  bold_runs <- as.numeric(gsub(".*run-([0-9]+)_masked.*", "\\1", bold_files))
  
  # Restrict to first 3 runs only
  if (length(bold_runs) == 0) {
    cat("No BOLD runs found for subject", sub, "- skipping.\n")
    next
  }
  bold_runs <- sort(unique(bold_runs))
  use_runs <- head(bold_runs, 3)
  
  # Check timing files
  timing_path <- file.path(timing_dir, sub, "Timing_Files_Pre")
  timing_files <- if (dir.exists(timing_path))
    list.files(timing_path, pattern = "_TIMING\\.txt$", full.names = TRUE)
  else character(0)
  timing_runs <- as.numeric(stringr::str_extract(timing_files, "(?<=run[-_]?)(\\d+)"))
  valid_runs <- intersect(use_runs, timing_runs)
  if (length(valid_runs) == 0) {
    cat("No valid runs for subject", sub, "- skipping.\n")
    next
  }
  
  cat("Using runs:", paste(valid_runs, collapse = ", "), "\n")
  
  for (m in valid_runs) {
    cat("Processing run", m, "for subject", sub, "\n")
    
    bold_path <- file.path(func_dir, paste0("sub-", sub, "_ses-1_task-GNG_run-", m, "_masked.nii.gz"))
    conf_path <- file.path(func_dir, paste0("sub-", sub, "_ses-1_task-GNG_acq-epi_run-", m, "_desc-confounds_timeseries.tsv"))
    
    if (!file.exists(bold_path) | !file.exists(conf_path)) {
      cat("Missing BOLD or confound for run", m, "- skipping run.\n")
      next
    }
    
    # Load BOLD and motion-only confounds
    dat <- readNifti(bold_path)
    confounds <- read.table(conf_path, header = TRUE)
    X <- regress_param_motion(confounds)
    
    # Extract VOIs
    VOI_list <- list()
    for (i in 1:3) {
      Vxyz <- t(which(mask == i, arr.ind = TRUE))
      BOLD <- matrix(NA, nrow = dim(dat)[4], ncol = ncol(Vxyz))
      for (j in 1:ncol(BOLD)) BOLD[, j] <- dat[Vxyz[1, j], Vxyz[2, j], Vxyz[3, j], ]
      BOLD <- apply(BOLD, 2, function(y) MotionRegress(y, X))
      sds <- apply(BOLD, 2, sd)
      BOLD <- BOLD[, sds > 0, drop = FALSE]
      
      if (ncol(BOLD) == 0) {
        cat("No valid voxels for VOI", i, "in run", m, "- skipping run.\n")
        next
      }
      VOI_list[[i]] <- get_BOLD_eigenvariate(BOLD)
    }
    
    y_obs <- dplyr::bind_cols(VOI_list)
    colnames(y_obs) <- c("MFG", "Insula", "Precuneus/PCC")
    y_obs[, 3] <- -1 * y_obs[, 3]  # negative for task− region
    scale <- max(y_obs) - min(y_obs)
    scale <- 4 / max(scale, 4)
    y_obs <- y_obs * scale
    
    # --- Stimulus files ---
    timing_prefix <- file.path(timing_dir, sub, "Timing_Files_Pre", paste0("sub-", sub, "_GoNoGo_"))
    NI  <- paste0(timing_prefix, "NI_run",  m, "-1_TIMING.txt")
    SNI <- paste0(timing_prefix, "SNI_run", m, "-1_TIMING.txt")
    SI  <- paste0(timing_prefix, "SI_run",  m, "-1_TIMING.txt")
    NNI <- paste0(timing_prefix, "NNI_run", m, "-1_TIMING.txt")
    
    task_files <- c(NI, SNI, SI, NNI)
    task_list <- lapply(task_files[file.exists(task_files)], read.table)
    task <- dplyr::bind_rows(task_list) %>% dplyr::arrange(V1)
    
    nscan <- nrow(y_obs)
    times <- seq(0, TR * nscan, TR)
    
    # Adjust for dropped volumes at start
    drop_vols <- 3  # first 3 TRs dropped
    time_shift <- drop_vols * TR
    
    # Adjust event onset times
    inhibition_all <- round((task$V1 - time_shift) / TR)
    dur_all <- round(task$V2 / TR)
    inhibition_all[inhibition_all < 0] <- 0  # clip negative onsets
    
    u_all <- stim.indicator(inhibition_all, dur_all, times)
    
    # Similarly for inhibition-only blocks
    inhibit_task_list <- lapply(list(NI, SI), function(f) {
      if (file.exists(f)) read.table(f) else NULL
    })
    inhibit_task <- dplyr::bind_rows(inhibit_task_list)
    
    inhibition_inhibit <- round((inhibit_task$V1 - time_shift) / TR)
    dur_inhibit <- round(inhibit_task$V2 / TR)
    inhibition_inhibit[inhibition_inhibit < 0] <- 0  # clip negative onsets
    
    u_inhibit <- stim.indicator(inhibition_inhibit, dur_inhibit, times)
    
    # Combine into scans x 2 matrix
    u <- cbind(u_all, u_inhibit)
    
    # Check that DMN region is negatively correlated with inhibition in case sign flip from PCA didn't line up right
    if (cor(y_obs[, "Precuneus/PCC"], u[,2]) > 0) {
      y_obs[, "Precuneus/PCC"] <- -y_obs[, "Precuneus/PCC"]
    }
    
    # Save .RData
    dat <- list(times = times[-1], u = u, y_obs = y_obs, scale = scale)
    save(dat, file = file.path(data_dir, paste0("sub-", sub, "_ses1_run", m, "_full_roi.RData")))
    
    # Plotting
    fig_dir <- file.path(data_dir, "Figures")
    if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
    png(file.path(fig_dir, paste0("sub-", sub, "_ses1_run", m, "_full_roi.png")),
        width = 900, height = 900, type = "cairo")
    
    par(mfrow = c(4, 1), mar = c(4, 4, 2, 1), oma = c(0, 0, 3, 0))
    # VOIs
    plot(y_obs[, 1], type = "l", lwd = 1.5, col = "darkred",
         xlab = "Scans", ylab = "MFG (Task+)", ylim = c(min(y_obs), max(y_obs)))
    abline(h = 0, col = "gray60", lty = 2)
    
    plot(y_obs[, 2], type = "l", lwd = 1.5, col = "darkred",
         xlab = "Scans", ylab = "Insula (Task+)", ylim = c(min(y_obs), max(y_obs)))
    abline(h = 0, col = "gray60", lty = 2)
    
    plot(y_obs[, 3], type = "l", lwd = 1.5, col = "darkblue",
         xlab = "Scans", ylab = "Precuneus/PCC (Task−)", ylim = c(min(y_obs), max(y_obs)))
    abline(h = 0, col = "gray60", lty = 2)
    
    # Stimulus panel
    plot(u[,1], type = "s", col = "black", lwd = 1.5, xlab = "Scans",
         ylab = "Stimulus", ylim = c(0,1), yaxt = "n", lty = 2) # general task dotted
    lines(u[,2]*0.5, type = "s", col = "red", lwd = 1.5, lty = 1)  # inhibition solid, half height
    axis(2, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"))
    abline(h = 0, col = "gray70")
    
    mtext(paste0("Subject ", sub, " — Run ", m, " (Full ROI)"), line = 1, outer = TRUE, cex = 1.3)
    dev.off()
    
    cat("Saved run", m, "for subject", sub, "\n")
  }
  
  cat("Finished subject:", sub, "\n")
}

########################################################################
