# Packages
library(cdcm)

dir.create("Simulation/Data_NoAgg",showWarnings = F,recursive = T)

############ Setting parameters #######################################
# 2 nodes, 2 experimental inputs
m = 2; n_u = 2

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
  
  # Seeds for reproducible noise
  snr_seeds <- sample.int(1e7,length(SNR_vals))
  
  for (j in 1:length(SNR_vals)){
    
    # Seed only affects noise
    sim_seed <- set.seed(rep_seeds[i] + snr_seeds[j])
    
    # Set SNR
    SNR = SNR_vals[j]
    
    # Simulate
    sim <- simulate_cdcm(
      nu = nu,
      hypothesis_idxs = idxs,
      nscan = length(times)-1,
      m = m,
      n_u = n_u,
      TR = (times[2] - times[1]),
      U = u,
      SNR = SNR,
      seed = sim_seed,
      z0 = z0
    )
    
    dat <- sim$simulated_data
    
    # Make as proper data objects and export
    dat$SNR <- SNR
    save(dat, file = paste0("Simulation/Data_NoAgg/snr",j,"_",i,".RData"))
    
  }
  
}



