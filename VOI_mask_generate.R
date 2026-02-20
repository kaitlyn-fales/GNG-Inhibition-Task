# Creation of VOI mask for GNG Inhibition Task with DCM

# Load package
library(RNifti)

# Initialize array of appropriate size
mask <- array(0, dim = c(67,80,46))

# Store coordinates for centers of VOIs (converted from MNI coordinates)
MFG <- c(47,46,31)
PCC <- c(34,29,27)

# Add MFG into mask with a radius of ~8mm 
center <- MFG
size <- c(6, 6, 4)  # number of voxels along X, Y, Z

x_range <- (center[1] - floor((size[1]-1)/2)):(center[1] + ceiling((size[1]-1)/2))
y_range <- (center[2] - floor((size[2]-1)/2)):(center[2] + ceiling((size[2]-1)/2))
z_range <- (center[3] - floor((size[3]-1)/2)):(center[3] + ceiling((size[3]-1)/2))

cube_voxels <- expand.grid(x=x_range, y=y_range, z=z_range)

for (i in 1:nrow(cube_voxels)) {
  mask[cube_voxels$x[i], cube_voxels$y[i], cube_voxels$z[i]] <- 1
}

# Add PCC into mask with a radius of ~8mm 
center <- PCC
size <- c(6, 6, 4)  # number of voxels along X, Y, Z

x_range <- (center[1] - floor((size[1]-1)/2)):(center[1] + ceiling((size[1]-1)/2))
y_range <- (center[2] - floor((size[2]-1)/2)):(center[2] + ceiling((size[2]-1)/2))
z_range <- (center[3] - floor((size[3]-1)/2)):(center[3] + ceiling((size[3]-1)/2))

cube_voxels <- expand.grid(x=x_range, y=y_range, z=z_range)

for (i in 1:nrow(cube_voxels)) {
  mask[cube_voxels$x[i], cube_voxels$y[i], cube_voxels$z[i]] <- 2
}

# Load in sample data file (from GNG task data, might need to update paths)
file_path <- "/storage/group/alh98/default/VLN_BIDS/fmri.analysis.bai.alh/derivatives/derfmrinf.23_trimmed/sub-120/ses-1/func"
dat <- readNifti(paste(file_path,"sub-120_ses-1_task-GNG_run-1_masked.nii.gz",sep = "/"))

# Save mask and use data as template for header file
writeNifti(mask, "mask.nii", template = dat)






