# Script to take the .RData files and make them Matlab compatible for SPM

library(R.matlab)

# List all RData files
files <- list.files(path = paste(getwd(),"Data",sep = "/"), pattern = "\\.RData$")

for (i in 1:length(files)){
  
  # Subject specs
  file_name <- basename(files[i])
  
  # Match sub ID and run number, ignoring anything between run number and .RData
  matches <- str_match(file_name, "sub-(\\d+)_ses1_run(\\d+)_(full|thresh)_roi\\.RData")
  
  sub <- matches[2]        # subject number
  run <- as.numeric(matches[3])  # run number
  type <- matches[4]       # "full" or "thresh"
  
  ses <- 1
  
  # Load file
  load(paste0("Data/",files[i]))
  
  # Version for MFG
  region <- "MFG"
  y_obs_MFG <- dat$y_obs[,c(region,"Precuneus/PCC")]
  
  writeMat(paste0("Data/sub_",sub,"_ses1_run",run,"_",region,"_",type,"_roi.mat"),
           U = dat$u,
           Y = y_obs_MFG)
  
  # Version for Insula
  region <- "Insula"
  y_obs_Insula <- dat$y_obs[,c(region,"Precuneus/PCC")]
  
  writeMat(paste0("Data/sub_",sub,"_ses1_run",run,"_",region,"_",type,"_roi.mat"),
           U = dat$u,
           Y = y_obs_Insula)
  
}

