# 08/01/2021 
# Lucy Vanes
# NNMF code for: The effects of neurobiology and environment on childhood outcomes following very preterm birth
# This code runs NNMF on voxelwise MRI data for a specific rank (identified using 2a_NNMF_test_rank.R)
# And outputs components back to .nii files for visualisation

library(NMF)
library(pracma)
library(dplyr)
library(neurobase)

# mask should be in mask_dir and final images to be used should be in data_dir
# mask and images must be in the same space with the same image dimensions!

mask_dir <- "C:/Users/vanes/OneDrive/Documents/GitHub/preterm-outcomes/data/NNMF_data/masks"
data_dir <- "C:/Users/vanes/OneDrive/Documents/GitHub/preterm-outcomes/data/NNMF_data/preprocessed"
networks_output_dir1 <- "C:/Users/vanes/OneDrive/Documents/GitHub/preterm-outcomes/data/NNMF_data/Network_vis_thresh25" 
networks_output_dir2 <- "C:/Users/vanes/OneDrive/Documents/GitHub/preterm-outcomes/data/NNMF_data/Network_vis_times100000" 

setwd(mask_dir)
mask <- readnii("AAL_atlas_mask_2mm.nii.gz") # read mask into R

mask_vec <- c(mask)                          # vectorise mask, this is now a vector of 0s and 1s
outside_mask <- which(mask_vec==0)           # define what is outside and inside mask
inside_mask <- which(mask_vec==1)
length(mask[mask==1])

setwd(data_dir)
scans <- dir(pattern="smoothedJacs_2mm_AAL_masked")   # some pattern to find your images
example_scan <- readnii(scans[1])
example_scan_vec <- c(example_scan)
scandata <- as.data.frame(matrix(nrow=length(inside_mask), ncol=length(scan)))

for (s in 1:length(scans)){
  print(scans[s])
  t <- readnii(scans[s])
  tvec <- c(t)
  tvec_masked <- tvec[mask==1]
  scandata[,s] <- tvec_masked
}

dat <- scandata
dat <- exp(dat) # exponentiate jacobians to get positive values
dat <- as.matrix(dat)
dim(dat)     # this should be voxels X subjects

#================================================================

best_rank <- 15   # choose rank (ascertained with previous script)

#=================================================================

# rerun NNMF with chosen number of ranks, using nndsdv
res_nmf <- nmf(dat, best_rank, method="lee",seed="nndsvd", .options="v")

# results NMF
V.hat <- fitted(res_nmf)
w_nmf <- basis(res_nmf)         # w is the "basis matrix" containing voxel weights for each component
h_nmf <- coef(res_nmf)          # h is the coefficient matrix containing subject weights for each component
h <- as.data.frame(t(h_nmf))
h$id <- substr(scans, 3, 6)
h$id <- factor(h$id)
dim(w_nmf)
dim(h_nmf)

#============================================
# write components back to nii 
#============================================
featuresMax <- extractFeatures(res_nmf, method="max") # see ?extractFeatures for details
features03 <- extractFeatures(res_nmf, 0.3)            # see ?extractFeatures for details

mask_vec <- c(mask)
outside_mask <- which(mask_vec==0)
inside_mask <- which(mask_vec==1)

for (r in 1:best_rank){
  #thresholding at 25%
  #===============================
  W <- as.data.frame(w_nmf)[,r]
  setwd(networks_output_dir1)
  W_zero <- which(almost.zero(W))
  W[W_zero] <- 0
  W_thresh <- range(W[W>0])[1] + (0.25*((range(W[W>0])[2]-range(W[W>0])[1])))
  W[W < W_thresh] <- 0
  new_nifti <- mask_vec
  new_nifti[inside_mask] <- W
  new_nifti <- niftiarr(mask, new_nifti)
  nii_name <- paste("NNMF_network_",r, ".nii.gz", sep="")
  writenii(new_nifti, nii_name)
  
  # unthresholded, times 100000
  #===============================
  W <- as.data.frame(w_nmf)[,r]
  setwd(networks_output_dir2)
  W <- W*100000
  new_nifti <- mask_vec
  new_nifti[inside_mask] <- W
  new_nifti <- niftiarr(mask, new_nifti)
  nii_name <- paste("NNMF_network_",r, ".nii.gz", sep="")
  writenii(new_nifti, nii_name)
}

#======================================================================================
# calculate weighted mean jacobian of each network
# using apply_weights.sh 
#======================================================================================
