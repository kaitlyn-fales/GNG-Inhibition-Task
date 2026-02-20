########### Preprocess data for inhibition task DCM ##################

# Packages
library(RNifti)
library(nlme)
library(dplyr)

# Directory for results
data_dir <- "Data" # to store .RData versions

################ List of subject IDs to loop through #######################
subjects <- c(120,121,124,128,134,143,160,172,173,176,181,
              184,196,199,214,215,223,227,230,247,252,256,258,265,266,268,
              271,275,276,277) # sub 110, 111, 152, 171 no sess 2 data
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

# Extract motion param for post-processing regression
regress_param <- function(x){
  names <- c("trans_x","trans_y","trans_z","rot_x","rot_y","rot_z")
  regress <- cbind(x$trans_x,x$trans_y,x$trans_z,x$rot_x,x$rot_y,x$rot_z)
  colnames(regress) <- names
  return(regress)
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

# Generate FSL-style motion regressors from fMRIPrep confounds
regress_param <- function(confounds) {
  # Extract 6 motion params
  motion6 <- confounds[, c("trans_x","trans_y","trans_z",
                           "rot_x","rot_y","rot_z"), drop = FALSE]
  
  # Fill first row of derivatives with 0
  deriv_cols <- c("trans_x_derivative1","trans_y_derivative1","trans_z_derivative1",
                  "rot_x_derivative1","rot_y_derivative1","rot_z_derivative1")
  confounds[1, deriv_cols] <- 0
  
  # Extract derivatives
  deriv6 <- confounds[, deriv_cols, drop = FALSE]
  
  # Combine motion + derivatives
  X <- cbind(motion6, deriv6)
  
  # Extract spike regressors if present
  spike_cols <- grep("^motion_outlier", colnames(confounds), value = TRUE)
  if(length(spike_cols) > 0){
    spikes <- confounds[, spike_cols, drop = FALSE]
    X <- cbind(X, spikes)
  }
  
  # Ensure numeric
  X[] <- lapply(X, function(x) as.numeric(as.character(x)))
  
  return(as.matrix(X))
}
########################################################################

########## Loop to process subject data ################################
for (k in 1:length(subjects)) {
  
  sub <- subjects[k]
  cat("\nProcessing subject:", sub, "\n")
  
  ########### Detect available runs dynamically ###########
  bold_files <- list.files(
    path = paste0(base_dir, "/sub-", sub, "/ses-2/func/"),
    pattern = "_masked\\.nii\\.gz$",
    full.names = TRUE
  )
  bold_runs <- as.numeric(gsub(".*run-([0-9]+)_masked.*", "\\1", bold_files))
  
  # ---- Timing detection ----
  timing_path <- paste0(timing_dir, "/", sub, "/Timing_Files_Post/")
  
  if (dir.exists(timing_path)) {
    timing_files <- list.files(
      path = timing_path,
      pattern = "_TIMING\\.txt$",
      full.names = TRUE
    )
    
    if (length(timing_files) == 0) {
      cat("No timing files found for subject", sub, "\n")
      timing_runs <- numeric(0)
    } else {
      # Extract run numbers robustly
      timing_runs <- stringr::str_extract(timing_files, "(?<=run[-_]?)(\\d+)")
      timing_runs <- as.numeric(timing_runs)
      
      if (all(is.na(timing_runs))) {
        cat("Timing files found, but none matched 'run' pattern for subject", sub, "\n")
        cat("Example filenames:\n")
        print(head(basename(timing_files)))
        timing_runs <- numeric(0)
      }
    }
  } else {
    cat("Timing directory missing for subject", sub, "\n")
    timing_runs <- numeric(0)
  }
  
  valid_runs <- intersect(bold_runs, timing_runs)
  if (length(valid_runs) == 0) {
    cat("No valid runs found for subject", sub, "- skipping.\n")
    next
  }
  
  cat("Valid runs:", paste(valid_runs, collapse = ", "), "\n")
  
  ########### Extract first PC from each VOI ###########
  VOI_list <- list()
  
  for (i in 1:2) {
    BOLD_list <- list()
    
    for (m in valid_runs) {
      bold_path <- paste0(base_dir, "/sub-", sub, "/ses-2/func/sub-", sub,
                          "_ses-2_task-GNG_run-", m, "_masked.nii.gz")
      conf_path <- paste0(base_dir, "/sub-", sub, "/ses-2/func/sub-", sub,
                          "_ses-2_task-GNG_acq-epi_run-", m,
                          "_desc-confounds_timeseries.tsv")
      
      if (!file.exists(bold_path) | !file.exists(conf_path)) {
        cat("Skipping missing run", m, "for subject", sub, "\n")
        next
      }
      
      dat <- readNifti(bold_path)
      confounds <- read.table(conf_path, header = TRUE)
      X <- regress_param(confounds)
      
      # Extract coordinates and build BOLD matrix
      Vxyz <- t(which(mask == i, arr.ind = TRUE))
      BOLD <- matrix(NA, nrow = dim(dat)[4], ncol = ncol(Vxyz))
      for (j in 1:ncol(BOLD)) {
        BOLD[, j] <- dat[Vxyz[1, j], Vxyz[2, j], Vxyz[3, j], ]
      }
      
      # Motion regression and filter out zero-variance voxels
      BOLD <- apply(BOLD, 2, function(y) MotionRegress(y, X = X))
      sds <- apply(BOLD, 2, sd)
      BOLD <- BOLD[, sds > 0]
      
      BOLD_list[[as.character(m)]] <- as.data.frame(BOLD)
    }
    
    if (length(BOLD_list) == 0) {
      cat("No valid BOLD runs for VOI", i, "in subject", sub, "\n")
      next
    }
    
    BOLD <- dplyr::bind_rows(BOLD_list)
    VOI_list[[i]] <- get_BOLD_eigenvariate(BOLD)
  }
  
  if (length(VOI_list) < 2) {
    cat("Skipping subject", sub, "due to insufficient VOIs.\n")
    next
  }
  
  ########### Prepare y_obs ###########
  y_obs <- dplyr::bind_cols(VOI_list)
  colnames(y_obs) <- c("MFG", "PCC")
  y_obs[, 2] <- -1 * y_obs[, 2]
  
  scale <- max(y_obs) - min(y_obs)
  scale <- 4 / max(scale, 4)
  y_obs <- y_obs * scale
  
  ########### Generate stimulus indicator ###########
  u_list <- list()
  
  for (q in valid_runs) {
    
    get_timing_file <- function(prefix, run, type) {
      file1 <- paste0(prefix, type, "_run", run, "-1_TIMING.txt")
      file2 <- paste0(prefix, type, "_run", run, "-2_TIMING.txt")
      file3 <- paste0(prefix, type, "_run", run, "-4_TIMING.txt")
      if (file.exists(file1)) return(file1)
      if (file.exists(file2)) return(file2)
      if (file.exists(file3)) return(file3)
      stop(paste0("Timing file not found for run ", run, " type ", type))
    }
    
    timing_prefix <- paste0(timing_dir, "/", sub, "/Timing_Files_Post/sub-", sub, "_GoNoGo_")
    
    NI  <- get_timing_file(timing_prefix, q, "NI")
    SNI <- get_timing_file(timing_prefix, q, "SNI")
    SI  <- get_timing_file(timing_prefix, q, "SI")
    NNI <- get_timing_file(timing_prefix, q, "NNI")
    
    if (!any(file.exists(c(NI, SNI, SI, NNI)))) {
      cat("No timing files for run", q, "— skipping.\n")
      next
    }
    
    timing_list <- list()
    for (file in c(NI, SNI, SI, NNI)) {
      if (file.exists(file)) timing_list[[file]] <- read.table(file)
    }
    
    task <- dplyr::bind_rows(timing_list) %>% dplyr::arrange(V1)
    inhibition <- ceiling(task$V1 / TR)
    dur <- ceiling(task$V2 / TR)
    
    # Define times per run
    nscan <- nrow(BOLD_list[[as.character(q)]])
    max_time <- TR * nscan
    times <- seq(0, max_time, TR)
    
    u_list[[as.character(q)]] <- as.data.frame(stim.indicator(inhibition, dur, times))
  }
  
  if (length(u_list) == 0) {
    cat("No valid stimulus files for subject", sub, "- skipping.\n")
    next
  }
  
  u <- dplyr::bind_rows(u_list)
  u <- as.numeric(u[, 1])
  
  ########### Save data ###########
  nscan <- nrow(y_obs)
  times <- seq(0, TR * nscan, TR)
  dat <- list(times = times[-1], u = u, y_obs = y_obs, scale = scale)
  
  ########### Plot data for sanity check ###########
  fig_dir <- paste0(data_dir, "/Figures")
  
  png(paste0(fig_dir, "/sub-", sub, "_ses2.png"), 
      width = 900, height = 800, type = "cairo")
  
  par(mfrow = c(3, 1), mar = c(4, 4, 2, 1), oma = c(0, 0, 3, 0))
  
  # Panel 1: MFG
  plot(y_obs[, 1], type = "l", lwd = 1.5, col = "darkred",
       xlab = "Scans", ylab = "MFG (Task+)",
       ylim = c(min(y_obs), max(y_obs)))
  abline(h = 0, col = "gray60", lty = 2)
  
  # Panel 2: PCC
  plot(y_obs[, 2], type = "l", lwd = 1.5, col = "darkblue",
       xlab = "Scans", ylab = "PCC (Task−)",
       ylim = c(min(y_obs), max(y_obs)))
  abline(h = 0, col = "gray60", lty = 2)
  
  # Panel 3: Stimulus timing
  plot(u, type = "s", col = "black", lwd = 1.5,
       xlab = "Scans", ylab = "Stimulus (u)",
       ylim = c(-0.1, 1.1), yaxt = "n")
  axis(2, at = c(0, 1), labels = c("Off", "On"))
  abline(h = 0, col = "gray70")
  
  mtext(paste0("Subject ", sub, " — Session 2"), line = 1, outer = TRUE, cex = 1.3)
  
  dev.off()
  cat("Saved figure for subject:", sub, "\n")
  
  save(dat, file = paste0(data_dir, "/sub-", sub, "_ses2.RData"))
  cat("Finished subject:", sub, "\n")
}

########################################################################
