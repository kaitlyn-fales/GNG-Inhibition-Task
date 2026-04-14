########### Preprocess data for inhibition task DCM with GLM voxel selection ##################

library(RNifti)
library(nlme)
library(dplyr)

# ---------------- User settings ----------------
subjects <- c(120,121,124,128,134,143,152,160,171,172,173,181,
              184,196,199,214,215,223,227,230,247,252,256,258,265,266,268,
              271,275,276,277) # need to use p = 0.1 for sub 199 and 215, run 2

base_dir <- "/storage/group/alh98/default/VLN_BIDS/fmri.analysis.bai.alh/derivatives/derfmrinf.23_trimmed"
timing_dir <- "/storage/group/alh98/default/VLN_BIDS/participant_data_complete"
data_dir <- "Data"   # where to save .RData
mask_file <- "mask.nii"
TR <- 2
pval_thresh <- 0.05   # uncorrected threshold for voxel selectionx

# Exploratory voxel-count settings (does NOT affect current VOI selection)
count_p_grid <- c(0.05, 0.10, 0.15)
min_count_vox <- 1
roi_names <- c("MFG", "Insula", "Precuneus/PCC")

# Collect one combined df across all processed subject/run/ROI combinations
voxel_count_summary <- list()
# -------------------------------------------------

# ---------------- Functions --------------------
stim.indicator <- function(os, dur, times){
  stim <- rep(0,length(times)-1)
  for(i in seq_along(os)){
    stim[os[i]:(os[i]+dur[i]-1)] <- 1
  }
  stim
}

MotionRegress <- function(y, X){
  if(var(y)!=0){
    df <- data.frame(y, X)
    residuals <- as.numeric(gls(y~., data=df, correlation=corAR1())$residuals)
  } else residuals <- NA
  residuals
}

get_BOLD_eigenvariate <- function(BOLD){
  pca <- prcomp(BOLD, center=TRUE, scale.=FALSE)
  V <- pca$rotation[,1]
  U <- pca$x[,1]
  d <- sign(sum(V)); if(d==0) d <- 1
  (U*d)/sqrt(ncol(BOLD))
}

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

spm_hrf <- function(TR, dt=TR/16, p=c(6,16,1,1,6,0,32)){
  t <- seq(0, p[7], by=dt)
  h1 <- (t^(p[1]/p[3]) * exp(-(t-p[1])/p[3]))
  h2 <- (t^(p[2]/p[4]) * exp(-(t-p[2])/p[4]))
  h <- h1 - h2/p[5]; h[t<0]<-0
  h/sum(h)
}

convolve_hrf <- function(stim_onoff, TR, hrf){
  rep_factor <- 16
  stim_hi <- rep(stim_onoff, each=rep_factor)
  conv_hi <- convolve(stim_hi, rev(hrf), type="open")[1:length(stim_hi)]
  conv_down <- tapply(conv_hi, rep(1:length(stim_onoff), each=rep_factor), mean)
  as.numeric(conv_down)
}

fit_voxel <- function(y, design, coef_name = "Inhibit") {
  # safe voxelwise fit: returns beta, t-value, and p-value for coef_name
  if (any(is.na(y)) || var(y) == 0) return(c(beta = NA, tval = NA, pval = NA))
  
  df <- data.frame(y = y, design)
  fmla <- as.formula(paste("y ~ -1 +", paste(colnames(design), collapse = " + ")))
  
  # Try GLS with AR1; fallback to lm on error
  res_try <- try(gls(fmla, data = df, correlation = corAR1(form = ~1), method = "ML"),
                 silent = TRUE)
  if (!inherits(res_try, "try-error")) {
    ttab <- summary(res_try)$tTable
    if (coef_name %in% rownames(ttab)) {
      b <- ttab[coef_name, "Value"]
      t <- ttab[coef_name, "t-value"]
      p <- ttab[coef_name, "p-value"]
    } else {
      b <- NA; t <- NA; p <- NA
    }
  } else {
    # fallback to ordinary lm
    res_lm <- try(lm(fmla, data = df), silent = TRUE)
    if (inherits(res_lm, "try-error")) {
      b <- NA; t <- NA; p <- NA
    } else {
      sm <- summary(res_lm)
      if (coef_name %in% rownames(sm$coefficients)) {
        b <- sm$coefficients[coef_name, "Estimate"]
        # lm uses "t value" and "Pr(>|t|)"
        t <- sm$coefficients[coef_name, "t value"]
        p <- sm$coefficients[coef_name, "Pr(>|t|)"]
      } else {
        b <- NA; t <- NA; p <- NA
      }
    }
  }
  
  c(beta = b, tval = t, pval = p)
}

get_timing_file <- function(prefix, run, type){
  file1 <- paste0(prefix,type,"_run",run,"-1_TIMING.txt")
  file2 <- paste0(prefix,type,"_run",run,"-2_TIMING.txt")
  if(file.exists(file1)) return(file1)
  if(file.exists(file2)) return(file2)
  stop(paste0("Timing file not found for run ",run," type ",type))
}

summarize_voxel_counts <- function(res_df,
                                   p_grid = c(0.05, 0.10, 0.15),
                                   min_vox = 1) {
  # Use t-value sign for consistent activation/deactivation labeling
  pos_idx_all <- which(!is.na(res_df$p_uncorrected) &
                         !is.na(res_df$tval) &
                         res_df$tval > 0)
  
  neg_idx_all <- which(!is.na(res_df$p_uncorrected) &
                         !is.na(res_df$tval) &
                         res_df$tval < 0)
  
  used_p <- tail(p_grid, 1)
  n_selected <- 0
  
  for (pth in p_grid) {
    pos_now <- which(!is.na(res_df$p_uncorrected) &
                       !is.na(res_df$tval) &
                       res_df$p_uncorrected <= pth &
                       res_df$tval > 0)
    
    neg_now <- which(!is.na(res_df$p_uncorrected) &
                       !is.na(res_df$tval) &
                       res_df$p_uncorrected <= pth &
                       res_df$tval < 0)
    
    n_now <- length(pos_now) + length(neg_now)
    used_p <- pth
    n_selected <- n_now
    
    if (n_now >= min_vox) break
  }
  
  activated_idx <- which(!is.na(res_df$p_uncorrected) &
                           !is.na(res_df$tval) &
                           res_df$p_uncorrected <= used_p &
                           res_df$tval > 0)
  
  deactivated_idx <- which(!is.na(res_df$p_uncorrected) &
                             !is.na(res_df$tval) &
                             res_df$p_uncorrected <= used_p &
                             res_df$tval < 0)
  
  noisy_idx <- setdiff(seq_len(nrow(res_df)),
                       union(activated_idx, deactivated_idx))
  
  n_total <- nrow(res_df)
  
  data.frame(
    n_total = n_total,
    n_activated = length(activated_idx),
    n_deactivated = length(deactivated_idx),
    n_noisy = length(noisy_idx),
    prop_activated = length(activated_idx) / n_total,
    prop_deactivated = length(deactivated_idx) / n_total,
    prop_noisy = length(noisy_idx) / n_total,
    count_p_used = used_p
  )
}
# -------------------------------------------------

mask <- readNifti(mask_file)

for(sub in subjects){
  cat("\nProcessing subject:", sub,"\n")
  
  func_dir <- file.path(base_dir, paste0("sub-",sub),"ses-1","func")
  bold_files <- list.files(func_dir, pattern="_masked\\.nii\\.gz$", full.names=TRUE)
  bold_runs <- sort(unique(as.numeric(gsub(".*run-([0-9]+)_masked.*","\\1",bold_files))))
  use_runs <- head(bold_runs,3)
  
  timing_path <- file.path(timing_dir, sub, "Timing_Files_Pre")
  timing_files <- if(dir.exists(timing_path)) list.files(timing_path, pattern="_TIMING\\.txt$", full.names=TRUE) else character(0)
  timing_runs <- as.numeric(stringr::str_extract(timing_files, "(?<=run[-_]?)(\\d+)"))
  valid_runs <- intersect(use_runs, timing_runs)
  if(length(valid_runs)==0){cat("No valid runs - skip.\n"); next}
  
  for(m in valid_runs){
    cat("Processing run", m,"\n")
    
    bold_path <- file.path(func_dir, paste0("sub-",sub,"_ses-1_task-GNG_run-",m,"_masked.nii.gz"))
    conf_path <- file.path(func_dir, paste0("sub-",sub,"_ses-1_task-GNG_acq-epi_run-",m,"_desc-confounds_timeseries.tsv"))
    if(!file.exists(bold_path) || !file.exists(conf_path)){cat("Missing files - skip\n"); next}
    
    dat4d <- readNifti(bold_path)
    confounds <- read.table(conf_path, header=TRUE)
    Xconf <- regress_param_motion(confounds)
    
    # --- Timing files ---
    timing_prefix <- file.path(timing_dir, sub, "Timing_Files_Pre", paste0("sub-",sub,"_GoNoGo_"))
    NI  <- get_timing_file(timing_prefix,m,"NI")
    SNI <- get_timing_file(timing_prefix,m,"SNI")
    SI  <- get_timing_file(timing_prefix,m,"SI")
    NNI <- get_timing_file(timing_prefix,m,"NNI")
    
    # --- u matrix (0/1) ---
    nscan <- dim(dat4d)[4]
    times <- seq(0, TR*nscan, TR)
    drop_vols <- 3
    time_shift <- drop_vols*TR
    
    # General task
    task_list <- lapply(list(NI,SNI,SI,NNI), function(f){read.table(f)})
    task_all <- dplyr::bind_rows(task_list)
    inhibition_all <- round((task_all$V1 - time_shift)/TR); dur_all <- round(task_all$V2/TR)
    inhibition_all[inhibition_all<0] <- 0
    u_all <- stim.indicator(inhibition_all,dur_all,times)
    
    # Inhibition blocks only (NI + SI)
    inhibit_list <- lapply(list(NI,SI), function(f){read.table(f)})
    task_inhibit <- dplyr::bind_rows(inhibit_list)
    inhibition_inhibit <- round((task_inhibit$V1 - time_shift)/TR); dur_inhibit <- round(task_inhibit$V2/TR)
    inhibition_inhibit[inhibition_inhibit<0] <- 0
    u_inhibit <- stim.indicator(inhibition_inhibit,dur_inhibit,times)
    
    u <- cbind(u_all, u_inhibit)
    
    # --- GLM: only inhibition contrast ---
    stim_hrf <- convolve_hrf(u_inhibit, TR, spm_hrf(TR))
    stim_hrf <- scale(stim_hrf, center=TRUE, scale=FALSE)
    design <- cbind(Intercept=1, Inhibit=stim_hrf, Xconf)
    colnames(design)[2] <- "Inhibit"
    
    # --- ROI PCA based on significant voxels ---
    # ---------- ROI loop: one-sided p with sign check, fallback to two-sided p only ----------
    y_obs_list <- list()
    for (i in 1:3) {
      Vxyz <- t(which(mask == i, arr.ind = TRUE))
      nvox <- if (ncol(Vxyz) == 0) 0 else ncol(Vxyz)
      results <- matrix(NA, nrow = nvox, ncol = 6)
      colnames(results) <- c("x", "y", "z", "beta", "tval", "p_uncorrected")
      
      if (nvox == 0) {
        cat("Warning: ROI", i, "has zero voxels in mask\n")
        y_obs_list[[i]] <- rep(NA, nscan)
        next
      }
      
      for (j in seq_len(nvox)) {
        v <- Vxyz[, j]
        y <- dat4d[v[1], v[2], v[3], ]
        fit <- fit_voxel(y, design, coef_name = "Inhibit")
        results[j, ] <- c(v[1], v[2], v[3], fit["beta"], fit["tval"], fit["pval"])
      }
      
      res_df <- as.data.frame(results, stringsAsFactors = FALSE)
      res_df$beta <- as.numeric(res_df$beta)
      res_df$tval <- as.numeric(res_df$tval)
      res_df$p_uncorrected <- as.numeric(res_df$p_uncorrected)
      
      # Exploratory voxel-count summary (separate from main VOI selection)
      count_row <- summarize_voxel_counts(
        res_df = res_df,
        p_grid = count_p_grid,
        min_vox = min_count_vox
      )
      
      voxel_count_summary[[length(voxel_count_summary) + 1]] <- data.frame(
        subject = sub,
        run = m,
        roi = roi_names[i],
        count_row
      )
      
      # Compute one-sided p: p_one = p_two / 2 for voxels with expected sign, otherwise set large
      p_one <- rep(Inf, nrow(res_df))
      if (i %in% c(1, 2)) {
        # ROIs 1 and 2 expected positive
        ok_pos <- which(!is.na(res_df$p_uncorrected) & !is.na(res_df$beta) & res_df$beta > 0)
        p_one[ok_pos] <- res_df$p_uncorrected[ok_pos] / 2
        roi_label <- "positive"
      } else {
        # ROI 3 expected negative
        ok_neg <- which(!is.na(res_df$p_uncorrected) & !is.na(res_df$beta) & res_df$beta < 0)
        p_one[ok_neg] <- res_df$p_uncorrected[ok_neg] / 2
        roi_label <- "negative"
      }
      
      # Primary selection: one-sided p <= threshold AND correct sign (already encoded in p_one)
      pass_one <- which(p_one <= pval_thresh)
      
      # Fallback 1: if none pass one-sided rule, use two-sided p-only voxels
      if (length(pass_one) == 0) {
        pass_two <- which(!is.na(res_df$p_uncorrected) & res_df$p_uncorrected <= pval_thresh)
        if (length(pass_two) == 0) {
          cat(sprintf("ROI %d: no voxels pass one-sided or two-sided p<=%.4f; skipping ROI\n", i, pval_thresh))
          y_obs_list[[i]] <- rep(NA, nscan)
          next
        } else {
          sig_idx <- pass_two
          cat(sprintf("ROI %d: no voxels passed one-sided %s test; using two-sided p-only voxels (n=%d)\n",
                      i, roi_label, length(sig_idx)))
        }
      } else {
        sig_idx <- pass_one
        cat(sprintf("ROI %d: %d voxels passed one-sided p<=%.4f (direction=%s)\n",
                    i, length(sig_idx), pval_thresh, roi_label))
      }
      
      # Extract selected voxel coordinates
      Vsig <- Vxyz[, sig_idx, drop = FALSE]
      BOLD <- matrix(NA, nrow = nscan, ncol = ncol(Vsig))
      for (j in 1:ncol(Vsig)) BOLD[, j] <- dat4d[Vsig[1, j], Vsig[2, j], Vsig[3, j], ]
      
      # Regress motion (and cosines if you included them earlier in Xconf)
      BOLD <- apply(BOLD, 2, function(y) MotionRegress(y, Xconf))
      sds <- apply(BOLD, 2, sd)
      BOLD <- BOLD[, sds > 0, drop = FALSE]
      
      if (ncol(BOLD) == 0) {
        cat(sprintf("ROI %d: all selected voxels dropped after motion regression; skipping ROI\n", i))
        y_obs_list[[i]] <- rep(NA, nscan)
        next
      }
      
      VOI <- get_BOLD_eigenvariate(BOLD)
      y_obs_list[[i]] <- VOI
    }
    # -------------------------------------------------------------------------------
    
    y_obs <- dplyr::bind_cols(y_obs_list)
    colnames(y_obs) <- c("MFG","Insula","Precuneus/PCC")
    
    # PCC sanity check
    if(cor(y_obs[,"Precuneus/PCC"], u[,2]) > 0){
      y_obs[,"Precuneus/PCC"] <- -y_obs[,"Precuneus/PCC"]
    }
    
    # Scale VOIs like before
    scale <- max(y_obs) - min(y_obs)
    scale <- 4 / max(scale, 4)
    y_obs <- y_obs * scale
    
    # Save
    dat <- list(times = times[-1], u = u, y_obs = y_obs, scale = scale)
    save(dat, file = file.path(data_dir, paste0("sub-", sub, "_ses1_run", m, "_thresh_roi.RData")))
    
    # Plotting
    fig_dir <- file.path(data_dir, "Figures")
    if(!dir.exists(fig_dir)) dir.create(fig_dir, recursive=TRUE)
    png(file.path(fig_dir, paste0("sub-",sub,"_ses1_run",m,"_thresh_roi.png")),
        width=900, height=900, type="cairo")
    par(mfrow=c(4,1), mar=c(4,4,2,1), oma=c(0,0,3,0))
    plot(y_obs[,1], type="l", lwd=1.5, col="darkred", xlab="Scans", ylab="MFG (Task+)", ylim=c(min(y_obs), max(y_obs)))
    abline(h=0,col="gray60",lty=2)
    plot(y_obs[,2], type="l", lwd=1.5, col="darkred", xlab="Scans", ylab="Insula (Task+)", ylim=c(min(y_obs), max(y_obs)))
    abline(h=0,col="gray60",lty=2)
    plot(y_obs[,3], type="l", lwd=1.5, col="darkblue", xlab="Scans", ylab="Precuneus/PCC (Task−)", ylim=c(min(y_obs), max(y_obs)))
    abline(h=0,col="gray60",lty=2)
    plot(u[,1], type="s", col="black", lwd=1.5, xlab="Scans", ylab="Stimulus", ylim=c(0,1), yaxt="n", lty=2)
    lines(u[,2]*0.5, type="s", col="red", lwd=1.5, lty=1)
    axis(2, at=c(0,0.5,1), labels=c("0","0.5","1"))
    abline(h=0,col="gray70")
    mtext(paste0("Subject ",sub," — Run ",m," (GLM-selected Voxels)"), line=1, outer=TRUE, cex=1.3)
    dev.off()
    
    cat("Saved run", m,"for subject",sub,"\n")
  }
  
  cat("Finished subject", sub,"\n")
}

voxel_count_summary_df <- dplyr::bind_rows(voxel_count_summary)

save(voxel_count_summary_df,
     file = file.path(data_dir, "voxel_count_summary_df.RData"))
