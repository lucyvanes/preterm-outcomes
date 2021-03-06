# 08/01/2021 
# Lucy Vanes
# NNMF code for: The effects of neurobiology and environment on childhood outcomes following very preterm birth
# This code runs NNMF on voxelwise MRI data for a specified range of ranks
# and outputs data to assess the optimal rank

library(NMF)
library(NNLM)
library(pracma)
library(dplyr)
library(neurobase)
library(parallel)

#=====================================
#           prepare data
#=====================================
mask_dir <- "C:/Users/vanes/OneDrive/Documents/GitHub/preterm-outcomes/data/NNMF_data/masks"
data_dir <- "C:/Users/vanes/OneDrive/Documents/GitHub/preterm-outcomes/data/NNMF_data/preprocessed"
output_dir <- "C:/Users/vanes/OneDrive/Documents/GitHub/preterm-outcomes/data/NNMF_data/output"
permuted_output_dir <- "C:/Users/vanes/OneDrive/Documents/GitHub/preterm-outcomes/data/NNMF_data/permuted_output"
  
setwd(mask_dir)
mask <- readnii("AAL_atlas_mask_2mm.nii.gz")       # name of mask (can be whole-brain mask)
mask_vec <- c(mask)
outside_mask <- which(mask_vec==0)
inside_mask <- which(mask_vec==1)
length(mask[mask==1])

setwd(data_dir)
scans <- dir(pattern="smoothedJacs")
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
dat <- exp(dat) # exponentiate jacobians
dat <- as.matrix(dat)

# random permutation of matrix
dat_permut <- dat
set.seed(12345)                     
for (c in 1:ncol(dat)){
  dat_permut[,c] <- dat[sample(nrow(dat)),c]
}

#===================================================
# run NNMF over number of ranks - ORIGINAL DATA
#===================================================

ranks <- seq(2,20,1) # choose ranks to run over
numCores <- 5 # choose number of cores for parallelising; or use detectCores()

setwd(output_dir)
source("C:/Users/vanes/OneDrive/Documents/GitHub/preterm-outcomes/2_NNMF/nnmf_test_ranks.R")
rank_res <- mclapply(ranks,nnmf_test_ranks, data=dat,  WH_true=1, mc.cores=numCores)

# save results
#=========================================
rank_res_data <- data.frame("rank"=ranks)
rank_res_data$recon_error_frob <- NA
rank_res_data$mean_WH_recon_error_frob <- NA
rank_res_data$RMSE_train <- NA
rank_res_data$RMSE_test <- NA

for (n in 1:length(rank_res)){
  print(n)
  rank_res_data$recon_error_frob[n] <- rank_res[[n]]$recon_error_frob
  rank_res_data$mean_WH_recon_error_frob[n] <- rank_res[[n]]$mean_WH_recon_error_frob
  rank_res_data$RMSE_train[n] <- rank_res[[n]]$RMSE_train
  rank_res_data$RMSE_test[n] <- rank_res[[n]]$RMSE_test
}

for (i in 1:length(ranks)){
  r <- ranks[i]
  print(r)
  wold_holdouts <- rank_res[[i]]$wold_holdouts
  write.csv(wold_holdouts, paste("wold_holdouts_rank_",r,".csv",sep=""), row.names=F, quote=F)
}

write.csv(rank_res_data,"rank_res_data_nnmf_WH_2to20.csv", row.names=F, quote=F)

#===================================================
# run NNMF over number of ranks - PERMUTED DATA
#===================================================
ranks <- seq(2,20,1)
numCores <- 5

setwd(permuted_output_dir)
source("C:/Users/vanes/OneDrive/Documents/GitHub/preterm-outcomes/2_NNMF/nnmf_test_ranks.R")
rank_res <- mclapply(ranks,nnmf_test_ranks, data=dat_permut,  WH_true=1, mc.cores=numCores)

# save results
#=========================================
rank_res_data <- data.frame("rank"=ranks)
rank_res_data$recon_error_frob <- NA
rank_res_data$mean_WH_recon_error_frob <- NA
rank_res_data$RMSE_train <- NA
rank_res_data$RMSE_test <- NA

for (n in 1:length(rank_res)){
  print(n)
  rank_res_data$recon_error_frob[n] <- rank_res[[n]]$recon_error_frob
  rank_res_data$mean_WH_recon_error_frob[n] <- rank_res[[n]]$mean_WH_recon_error_frob
  rank_res_data$RMSE_train[n] <- rank_res[[n]]$RMSE_train
  rank_res_data$RMSE_test[n] <- rank_res[[n]]$RMSE_test
}

for (i in 1:length(ranks)){
  r <- ranks[i]
  print(r)
  wold_holdouts <- rank_res[[i]]$wold_holdouts
  write.csv(wold_holdouts, paste("wold_holdouts_rank_",r,".csv",sep=""), row.names=F, quote=F)
}

write.csv(rank_res_data,"rank_res_data_permuted_nnmf_WH_2to20.csv", row.names=F, quote=F)

#===============================================================================
# Identify best rank from output
#===============================================================================
library(NMF)
library(NNLM)
library(pracma)
library(dplyr)
library(neurobase)

# Look at holdout CV in more detail

#==================
# original data
#==================
setwd(output_dir)
holdout_files <- dir(pattern="wold_holdouts_rank")
wold_holdouts <- NULL

for (f in holdout_files){
  rank <- as.numeric(substr(f, 20,21))
  wold_holdouts_f <- read.csv(f, header=T)
  wold_holdouts_f$rank <- rank
  wold_holdouts <- rbind(wold_holdouts_f,wold_holdouts)
}
wold_holdouts$rank <- factor(wold_holdouts$rank)
wold_holdouts <- wold_holdouts[order(wold_holdouts$rank),]
wold_holdouts$rank <- as.numeric(as.character(wold_holdouts$rank))
wold_holdouts$grad_frob <- NA
wold_holdouts$grad_rmse_test <- NA
wold_holdouts$grad_rmse_train <- NA

for (r in 1:50){
  wold_holdouts$grad_frob[wold_holdouts$rep==r] <- 
    gradient(wold_holdouts$recon_error_frob[wold_holdouts$rep==r])
  
  wold_holdouts$grad_rmse_test[wold_holdouts$rep==r] <- 
    gradient(wold_holdouts$rmse_test[wold_holdouts$rep==r])
  
  wold_holdouts$grad_rmse_train[wold_holdouts$rep==r] <- 
    gradient(wold_holdouts$rmse_train[wold_holdouts$rep==r])
}

wold_holdouts_aggr <- aggregate(recon_error_frob ~ rank, wold_holdouts, mean)
wold_holdouts_aggr$grad_frob <- gradient(wold_holdouts_aggr$recon_error_frob)
wold_holdouts_aggr$rmse_train <- aggregate(rmse_train ~ rank, wold_holdouts, mean)[,2]
wold_holdouts_aggr$rmse_test <- aggregate(rmse_test ~ rank, wold_holdouts, mean)[,2]
wold_holdouts_aggr$grad_rmse_train <- gradient(wold_holdouts_aggr$rmse_train)
wold_holdouts_aggr$grad_rmse_test <- gradient(wold_holdouts_aggr$rmse_test)

#=================
# permuted data
#=================
setwd(permuted_output_dir)
holdout_files_permut <- dir(pattern="wold_holdouts_rank")
wold_holdouts_permut <- NULL

for (f in holdout_files_permut){
  rank <- as.numeric(substr(f, 20,21))
  wold_holdouts_f <- read.csv(f, header=T)
  wold_holdouts_f$rank <- rank
  wold_holdouts_permut <- rbind(wold_holdouts_f,wold_holdouts_permut)
}
wold_holdouts_permut$rank <- factor(wold_holdouts_permut$rank)
wold_holdouts_permut <- wold_holdouts_permut[order(wold_holdouts_permut$rank),]
wold_holdouts_permut$rank <- as.numeric(as.character(wold_holdouts_permut$rank))
wold_holdouts_permut$grad_frob <- NA
wold_holdouts_permut$grad_rmse_test <- NA
wold_holdouts_permut$grad_rmse_train <- NA

for (r in 1:50){
  wold_holdouts_permut$grad_frob[wold_holdouts_permut$rep==r] <- 
    gradient(wold_holdouts_permut$recon_error_frob[wold_holdouts_permut$rep==r])
  
  wold_holdouts_permut$grad_rmse_test[wold_holdouts_permut$rep==r] <- 
    gradient(wold_holdouts_permut$rmse_test[wold_holdouts_permut$rep==r])
  
  wold_holdouts_permut$grad_rmse_train[wold_holdouts_permut$rep==r] <- 
    gradient(wold_holdouts_permut$rmse_train[wold_holdouts_permut$rep==r])
}

wold_holdouts_aggr_permut <- aggregate(recon_error_frob ~ rank, wold_holdouts_permut, mean)
wold_holdouts_aggr_permut$grad_frob <- gradient(wold_holdouts_aggr_permut$recon_error_frob)
wold_holdouts_aggr_permut$rmse_train <- aggregate(rmse_train ~ rank, wold_holdouts_permut, mean)[,2]
wold_holdouts_aggr_permut$rmse_test <- aggregate(rmse_test ~ rank, wold_holdouts_permut, mean)[,2]
wold_holdouts_aggr_permut$grad_rmse_train <- gradient(wold_holdouts_aggr_permut$rmse_train)
wold_holdouts_aggr_permut$grad_rmse_test <- gradient(wold_holdouts_aggr_permut$rmse_test)

#================================================
#                     Plots
#================================================
library(ggplot2)

wold_holdouts_aggr$Input <- "original"
wold_holdouts_aggr_permut$Input <- "permuted"
wold_holdouts_aggr_all <- rbind(wold_holdouts_aggr, wold_holdouts_aggr_permut)
wold_holdouts_aggr_all$Input <- factor(wold_holdouts_aggr_all$Input, levels=c("permuted","original"))

wold_holdouts$Input <- "original"
wold_holdouts_permut$Input <- "permuted"
wold_holdouts_all <- rbind(wold_holdouts, wold_holdouts_permut)
wold_holdouts_all$Input <- factor(wold_holdouts_all$Input, levels=c("permuted","original"))

mytheme <- theme(legend.background=element_rect(), 
      panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
      panel.background = element_rect(fill = "#d9d9d9", colour="black"))

# background #bdbdbd f0f0f0 #d9d9d9
# legend.position=c(0.138,0.08)

# mean
ggplot(data=wold_holdouts_aggr_all, aes(x=rank, y=recon_error_frob, color=Input)) + 
  geom_point(size=2) + scale_color_manual(values=c("#e41a1c","black")) +
  ylab("Mean reconstruction error") + xlab("Rank")  + 
  mytheme  + theme(legend.position=c(0.138,0.08))

# all points
ggplot(data=wold_holdouts_all, aes(x=rank, y=recon_error_frob, color=Input)) + 
  geom_point(pch=21) + scale_color_manual(values=c("#e41a1c","black")) +
  ylab("Reconstruction error") + xlab("Rank")  + 
  mytheme + theme(legend.position=c(0.12,0.11))


# Plot gradient
#======================================

library(dplyr)
ggplot(data=wold_holdouts_all %>%
         arrange(Input), 
       aes(x=rank, y=grad_frob, color=Input)) + 
  geom_point(pch=21) + scale_color_manual(values=c("#e41a1c","black")) +
  geom_point(data=wold_holdouts_aggr, aes(x=rank, y=grad_frob), color="white", size=2) +
  geom_point(data=wold_holdouts_aggr_permut, aes(x=rank, y=grad_frob), color="white", size=2) +
  ylab("Gradient of reconstruction error") + xlab("Rank")  + 
  geom_vline(xintercept = 15, col="blue", alpha=0.4) +
  mytheme + theme(legend.position=c(0.88, 0.11))


# TEST
#===============================================================================
for (r in unique(wold_holdouts$rank)){
  print(paste("rank",r,sep=" "))
  print("================", quote=F)
  print(t.test(wold_holdouts$grad_frob[wold_holdouts$rank==r], 
               wold_holdouts_permut$grad_frob[wold_holdouts_permut$rank==r]))
}

# --> at which rank is the difference between original and permuted reconstruction 
# error gradient no longer significant? Here it is 16, so I choose 15 as my 
# optimal rank

# Plot RMSE
#======================

plot_long <- reshape(wold_holdouts_all, direction="long",
                                  idvar=c("rep","rank","Input"),
                                  varying=list(c("grad_rmse_test","grad_rmse_train"),c("rmse_test","rmse_train")),
                                  v.names=c("grad_rmse","rmse"), timevar="train_or_test",
                                  time=c("test","train"),
                                  drop=c("reson_error_frob",
                                         "mse_train","mse_test","mae_train","mae_test",
                                         "grad_frob"))
plot_long$Data <- NA
plot_long$Data[plot_long$Input=="original" & plot_long$train_or_test=="train"] <- "Original train"
plot_long$Data[plot_long$Input=="original" & plot_long$train_or_test=="test"] <- "Original test"
plot_long$Data[plot_long$Input=="permuted" & plot_long$train_or_test=="train"] <- "Permuted train"
plot_long$Data[plot_long$Input=="permuted" & plot_long$train_or_test=="test"] <- "Permuted test"
plot_long$Data <- factor(plot_long$Data, levels=c("Permuted test","Permuted train","Original test","Original train"))


ggplot(data=plot_long %>%
         arrange(Data), 
       aes(x=rank, y=rmse, color=Data)) + 
  geom_point(pch=21) + scale_color_manual(values=c("#b2182b","#d6604d","#2166ac","#92c5de")) +
  ylab("RMSE") + xlab("Rank")  + 
  mytheme 

ggplot(data=plot_long %>%
         arrange(Data), 
       aes(x=rank, y=grad_rmse, color=Data)) + 
  geom_point(pch=21) + scale_color_manual(values=c("#b2182b","#d6604d","#2166ac","#92c5de")) +
  ylab("Gradient of RMSE") + xlab("Rank")  + 
  geom_vline(xintercept = 15, col="blue", alpha=0.4) +
  mytheme 