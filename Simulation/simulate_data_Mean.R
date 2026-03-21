# Packages
library(deSolve)
library(nlme)
library(tidyverse)

############ Setting parameters #######################################
# 2 nodes, 2 experimental inputs
m = 2; n_u = 2

# Number of voxels per ROI
nvoxel = 100

# Number of simulation replicates
nreps <- 50

# Set parameters
nu = list(nu_A = c(0.1,-0.3,0.05), 
          nu_B = c(0.15,-0.1), 
          nu_C = c(0.9))

# Indices of parameters
A_idxs <- matrix(c(1,1,
                   2,1,
                   2,2), byrow = T, ncol = 2)
B_idxs <- matrix(c(2,1,1,
                   2,2,1), byrow = T, ncol = 3)
C_idxs <- matrix(c(1,1), byrow = T, ncol = 2)

idxs <- list(A_idxs = A_idxs,
             B_idxs = B_idxs,
             C_idxs = C_idxs)

# Initial condition
z0 <- rep(0.1,m)

# Set vector of SNRs
SNR_vals <- seq(0.1,2,length.out = 8)

# Set pval threshold for the cherry picking task approach (one-sided)
pval_thresh = 0.01
########################################################################

############ Functions #################################################
# Function for structuring parameters
struct_paramMats = function(m, n_u, idxs, nu){
  
  A = matrix(data = 0,nrow = m, ncol = m)
  for(i in 1:nrow(idxs$A_idxs)){
    A[idxs$A_idxs[i,1], idxs$A_idxs[i,2]] = with(nu,nu_A[i])
  }
  diag(A) = -0.5*exp(diag(A))
  
  B = lapply(1:n_u, function(x){matrix(data = 0,nrow = m, ncol = m)})
  for(i in 1:nrow(idxs$B_idxs)){
    B[[idxs$B_idxs[i,1]]][idxs$B_idxs[i,2], idxs$B_idxs[i,3]] = with(nu,nu_B[i])
  }
  
  C = matrix(data=0, nrow = m, ncol = n_u)
  for(i in 1:nrow(idxs$C_idxs)){
    C[idxs$C_idxs[i,1], idxs$C_idxs[i,2]] = with(nu,nu_C[i])
  }
  
  return(list(A=A,B=B,C=C))
  
}

# ODE for neural activation that needs to be solved
linear <- function(t, z, params, input) {
  # t is the current time point in the integration, 
  # z is the current estimate of the variables in the ODE system
  A = params[["A"]]
  B = params[["B"]]
  C = params[["C"]]
  u = matrix(input(t),ncol=1)
  B_all = Reduce("+",lapply(1:nrow(u), function(i) u[i,1]*B[[i]]))
  dz <- (A+B_all)%*%z + C%*%u
  return(list(dz))
}

# Putting together to form U matrix of experimental inputs
input_u = function(t){
  c(approxfun(x = times,y = u[,1],rule = 2)(t),
    approxfun(x = times,y = u[,2],rule = 2)(t))
}

# hrf function (using the difference of two gammas - canonical)
HRF = function(t){dgamma(x = t, shape = 6,rate = 1) - (1/6)*dgamma(x = t,shape = 16,rate = 1)}
HRF_mu = function(mu,tp){
  convolve(mu, rev(HRF(tp)),type="open")[1:length(tp)]
}

# Function for simple voxelwise task glm and extract coefficient and pval
task_glm <- function(y, X){
  df <- data.frame(y, X)
  mdl <- gls(y~., data=df, correlation = corAR1())
  
  tval <- summary(mdl)$tTable[,"t-value"][2]
  pval <- summary(mdl)$tTable[,"p-value"][2]
  
  return(c(tval,pval))
}
########################################################################

# Randomly select a design u matrix from the real data (unique to each replicate)
subjects <- c(110,111,120,121,124,128,134,143,152,160,171,172,173,181,
              184,196,199,214,215,223,227,230,247,252,256,258,265,266,268,
              271,275,276,277)

# Global seed
set.seed(12345)
rep_seeds <- sample.int(1e7, nreps)

for (i in 1:nreps){
  
  set.seed(rep_seeds[i])
  sub <- sample(subjects, 1)
  run <- sample(1:3, 1)
  
  # Load data to extract u and times
  load(paste0("Data/sub-",sub,"_ses1_run",run,"_full_roi.RData"))
  u <- rbind(0,dat$u) # unique design
  times <- c(0,dat$times) # same for all subjects
  
  rm(dat) # get rid of data object
  
  # Model params
  paramMats = struct_paramMats(m = m,n_u = n_u,idxs = idxs,nu = nu)
  
  # Solve the ODE for z
  out_z <- ode(y = z0,
               times = times,
               func = linear,
               parms= with(paramMats, list(A=A, B=B, C=C)),
               input = input_u,
               atol = 1e-6, 
               rtol = 1e-6
  )
  
  # Get rid of first column of z (same as times)
  out_z <- out_z[,-1]
  
  # Hemodynamic model 
  y_signal = sapply(1:m, function(i) HRF_mu(out_z[-1,i],times[-1]))
  
  # Target signals
  target_MFG <- y_signal[,1]
  target_PCC <- y_signal[,2]
  
  # Generate voxel-specific scaling for replicate
  scale_factors <- 0.3 + (1.3-0.3)*rbeta(nvoxel*2,2,3)
  scale_factors_MFG <- scale_factors[1:100]
  scale_factors_PCC <- scale_factors[101:200]
  
  # Seeds for reproducible noise
  snr_seeds <- sample.int(1e7,length(SNR_vals))
  
  for (j in 1:length(SNR_vals)){
    
    # Seed only affects noise
    set.seed(rep_seeds[i] + snr_seeds[j])
    
    # Set SNR
    SNR = SNR_vals[j]
    
    MFG <- list()
    PCC <- list()
    
    # Simulate 200 voxels (100 per ROI)
    for (k in 1:nvoxel){
      
      # Add dampening/amplification - mimic a wide variety of voxels
      scale_factor <- c(scale_factors_MFG[k],scale_factors_PCC[k])
      
      # Multiply the voxel signal by scale_factor (different for each voxel)
      y_scaled <- sweep(y_signal, MARGIN = 2, STATS = scale_factor, FUN = "*")
      
      # Add noise to the voxel in each ROI
      variance <-  numeric()
      y_obs <- matrix(NA, nrow = length(times)-1, ncol = m)
      for (l in 1:m){
        variance[l] <- (var(y_scaled[,l])+(mean(y_scaled[,l]))^2)/SNR
        y_obs[,l] <- y_scaled[,l] + rnorm(length(times)-1, mean = 0, sd = sqrt(variance[l]))
      }
      
      # Take the first column and add to region 1 list (adding "voxels" to ROI)
      MFG[[k]] <- y_obs[,1]
      
      # Take the second column and add to region 2 list (adding "voxels" to ROI)
      PCC[[k]] <- y_obs[,2]
      
    }
    
    MFG <- suppressMessages(bind_cols(MFG))
    PCC <- suppressMessages(bind_cols(PCC))
    
    ### Using full ROI approach ###
    # Get mean time course for each ROI - full ROI approach using all voxels
    MFG_VOI_full <- rowMeans(MFG)
    PCC_VOI_full <- rowMeans(PCC) 
    
    # PCC sanity check, make sure negatively correlated with Inhibition stimulus
    if(cor(PCC_VOI_full, u[-1,2]) > 0){
      PCC_VOI_full <- -PCC_VOI_full
    }
    
    # Extract pre-model quantities
    # Pre-Stan metrics: PCA full ROI
    cor_MFG_mean_full <- cor(MFG_VOI_full, target_MFG)
    cor_PCC_mean_full <- cor(PCC_VOI_full, target_PCC)
    
    snr_eff_MFG_mean_full <- var(target_MFG) / var(MFG_VOI_full - target_MFG)
    snr_eff_PCC_mean_full <- var(target_PCC) / var(PCC_VOI_full - target_PCC)
    ###
    
    ### Using the thresholded approach ###
    # Convolve inhibition stimulus with HRF
    inhibit_stim <- HRF_mu(u[-1,2],times[-1])
    
    # Apply task glm for inhibition stimulus for each ROI
    task_result_MFG <- apply(MFG, 2, function(y) task_glm(y, inhibit_stim))
    task_result_PCC <- apply(PCC, 2, function(y) task_glm(y, inhibit_stim))
    
    # Get task significant thresholded voxels (one-sided test)
    sig_idx_MFG <- which(task_result_MFG[1,] > 0 & task_result_MFG[2,]/2 < pval_thresh, arr.ind = T)
    sig_idx_PCC <- which(task_result_PCC[1,] < 0 & task_result_PCC[2,]/2 < pval_thresh, arr.ind = T) # task deactivate
    
    # Proportion of voxels selected
    prop_selected_MFG <- length(sig_idx_MFG) / nvoxel
    prop_selected_PCC <- length(sig_idx_PCC) / nvoxel
    
    # Get first PC for each ROI - full ROI approach using all voxels
    MFG_VOI_thresh <- rowMeans(MFG[,sig_idx_MFG])
    PCC_VOI_thresh <- rowMeans(PCC[,sig_idx_PCC]) 
    
    # PCC sanity check, make sure negatively correlated with Inhibition stimulus
    if(cor(PCC_VOI_thresh, u[-1,2]) > 0){
      PCC_VOI_thresh <- -PCC_VOI_thresh
    }
    
    # Pre-Stan metrics: PCA thresholded ROI
    cor_MFG_mean_thresh <- cor(MFG_VOI_thresh, target_MFG)
    cor_PCC_mean_thresh <- cor(PCC_VOI_thresh, target_PCC)
    
    snr_eff_MFG_mean_thresh <- var(target_MFG) / var(MFG_VOI_thresh - target_MFG)
    snr_eff_PCC_mean_thresh <- var(target_PCC) / var(PCC_VOI_thresh - target_PCC)
    ###
    
    # Make as proper data objects and export
    y_obs_full <- cbind(MFG_VOI_full,PCC_VOI_full)
    dat <- list(times = times[-1], u = u[-1,], y_obs = y_obs_full, SNR = SNR)
    save(dat, file = paste0("Simulation/Data_Mean/full_roi_snr",j,"_",i,".RData"))
    
    y_obs_thresh <- cbind(MFG_VOI_thresh,PCC_VOI_thresh)
    dat <- list(times = times[-1], u = u[-1,], y_obs = y_obs_thresh, SNR = SNR)
    save(dat, file = paste0("Simulation/Data_Mean/thresh_roi_snr",j,"_",i,".RData"))
    
    # Save simulation metrics
    metrics_dat <- list(
      rep = i,
      snr_index = j,
      SNR = SNR,
      sub = sub,
      run = run,
      
      cor_MFG_mean_full = cor_MFG_mean_full,
      cor_PCC_mean_full = cor_PCC_mean_full,
      cor_MFG_mean_thresh = cor_MFG_mean_thresh,
      cor_PCC_mean_thresh = cor_PCC_mean_thresh,
      
      snr_eff_MFG_mean_full = snr_eff_MFG_mean_full,
      snr_eff_PCC_mean_full = snr_eff_PCC_mean_full,
      snr_eff_MFG_mean_thresh = snr_eff_MFG_mean_thresh,
      snr_eff_PCC_mean_thresh = snr_eff_PCC_mean_thresh,
      
      prop_selected_MFG = prop_selected_MFG,
      prop_selected_PCC = prop_selected_PCC
    )
    
    save(metrics_dat,
         file = paste0("Simulation/Metrics/mean_metrics_snr", j, "_", i, ".RData"))
    
  }
  
}



