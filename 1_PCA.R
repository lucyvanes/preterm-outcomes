
library(caret)


setwd("C:/Users/vanes/Dropbox/KCL_BMEIS/ePrime/data/outputs/")
dat <- read.csv("ePrime_ALL_OUTCOMES_data.csv", header=T)
dat$mother_edu <- factor(dat$mother_edu)

psych_vars <- c("z_ecbq_neg_affect","z_ecbq_surgency","z_ecbq_effort_control",
                "z_emque_emo_contagion","z_emque_attent_others","z_emque_prosocial",
                "z_sdq_emo","z_sdq_conduct","z_sdq_adhd","z_sdq_peer","z_sdq_prosocial",
                "z_adhd_inattention","z_adhd_hyper",
                "z_srs_soc_awareness","z_srs_soc_cog","z_srs_soc_comm","z_srs_soc_mot","z_srs_rrb")


cog_vars <- c("z_brief_inhibit","z_brief_shift","z_brief_emo_control","z_brief_wm",
                "z_brief_plan",
                "z_wppsi_verb_comp","z_wppsi_visuospatial","z_wppsi_fluid_res","z_wppsi_wm","z_wppsi_proc_speed",
                "z_wppsi_vocab","z_wppsi_nonverbal","z_wppsi_gen_abilities","z_wppsi_cog_prof")

control_vars <- c("age4","ses")
all_vars <- c("id",psych_vars, cog_vars, control_vars)

dat_prep <- dat[,all_vars]
dat_prep <- dat_prep[complete.cases(dat_prep),]

# correct certain variables for age at assessment and ses (IMD)
#=================================================================================
vars_reg1 <- c("z_ecbq_neg_affect","z_ecbq_surgency","z_ecbq_effort_control",
          "z_emque_emo_contagion","z_emque_attent_others","z_emque_prosocial",
          "z_sdq_emo","z_sdq_conduct","z_sdq_adhd","z_sdq_peer","z_sdq_prosocial",
          "z_adhd_inattention","z_adhd_hyper",
          "z_srs_soc_awareness","z_srs_soc_cog","z_srs_soc_comm","z_srs_soc_mot","z_srs_rrb",
          "z_brief_inhibit","z_brief_shift","z_brief_emo_control","z_brief_wm",
          "z_brief_plan")
vars_reg2 <- c("z_wppsi_verb_comp","z_wppsi_visuospatial","z_wppsi_fluid_res","z_wppsi_wm","z_wppsi_proc_speed",
               "z_wppsi_vocab","z_wppsi_nonverbal","z_wppsi_gen_abilities","z_wppsi_cog_prof")

for (v in vars_reg1){
  lm1 <- lm(paste(v, "~ age4 + ses"), dat_prep)
  dat_prep[v] <- lm1$resid
}
for (v in vars_reg2){
  lm1 <- lm(paste(v, "~ ses"), dat_prep)
  dat_prep[v] <- lm1$resid
}

dat_prep$age4 <- NULL
dat_prep$ses <- NULL



# Run PCA
#================
pcadat <- dat_prep


pcadat$id <- NULL

pcadat_scaled <- as.data.frame(scale(pcadat, center=T, scale=T))


# permutation testing:
# for each permutation, shuffle rows in each column, re-compute PCA
#==============================================================================

# http://bioops.info/2015/01/permutation-pca/

# the fuction to assess the significance of the principal components.
sign.pc<-function(x,R=5000,s=10, cor=T,...){
  # run PCA
  pc.out<-princomp(x,cor=cor,...)
  # the proportion of variance of each PC
  pve=(pc.out$sdev^2/sum(pc.out$sdev^2))[1:s]
  
  # a matrix with R rows and s columns that contains
  # the proportion of variance explained by each pc
  # for each randomization replicate.
  pve.perm<-matrix(NA,ncol=s,nrow=R)
  for(i in 1:R){
    # permutation each column
    x.perm<-apply(x,2,sample)
    # run PCA
    pc.perm.out<-princomp(x.perm,cor=cor,...)
    # the proportion of variance of each PC.perm
    pve.perm[i,]=(pc.perm.out$sdev^2/sum(pc.perm.out$sdev^2))[1:s]
  }
  # calcalute the p-values
  pval<-apply(t(pve.perm)>pve,1,sum)/R
  return(list(pve=pve,pval=pval))
}


pca_sign <- sign.pc(pcadat_scaled,cor=T)
pca_sign
# PC 1-4 are significant in full sample

#==============================================================================
# Original PCA

out <- princomp(pcadat_scaled)
summary(out)

PCs_orig <- as.data.frame(out$scores)
PC_nos <- which(pca_sign$pval < 0.05)
PC_nos

PCs_orig_sig <- PCs_orig[,PC_nos]
names(PCs_orig_sig) <- paste("orig_PC", PC_nos, sep="")

#=============================================================================
# k-fold cross-validation

k <- 10

pred_obs_cor <- vector(mode = "list", length = k)

for (k in 2:10){
  myfolds <- createFolds(pcadat_scaled[,1], k = k)
  
  pred_obs_cor[[k]] <- data.frame(matrix(nrow=k, ncol=length(PC_nos)))
  
  for (f in 1:k){
  
    pcadat_train <- pcadat_scaled[-myfolds[[f]], ]
    pcadat_test <- pcadat_scaled[myfolds[[f]], ]
  
    out_f <- princomp(pcadat_train)
    pred <- predict(out_f, newdata=pcadat_test)

    pred_obs_cor[[k]][f,] <- abs(diag(cor(PCs_orig_sig[myfolds[[f]],], pred[,PC_nos])))
  }
  
  # output plots here  
}
pred_obs_cor


# Split-half
#==================================================================

n_half <- round(dim(pcadat_scaled)[1] * 0.5)

H <- 1000 # how many split half replications

loadings_cor <- data.frame(matrix(nrow=H, ncol=length(PC_nos)))
loadings_sig <- data.frame(matrix(nrow=H, ncol=length(PC_nos)*2))
names(loadings_sig) <- c("half1_PC1","half1_PC2","half1_PC3","half1_PC4",
                         "half2_PC1","half2_PC2","half2_PC3","half2_PC4")
for (n in names(loadings_sig)){
  loadings_sig[n] <- 0
}

for (h in 1:H){
  print(paste("Split half ", h, sep=""), quote=F)
  print("============================")
  half1 <- sort(sample(1:dim(pcadat_scaled)[1], size=n_half))
  half2 <- setdiff(1:dim(pcadat_scaled)[1], half1)

  pcadat1 <- pcadat_scaled[half1,]
  pcadat2 <- pcadat_scaled[half2,]

  out1 <- princomp(pcadat1)
  out2 <- princomp(pcadat2)

  loadings1 <- out1$loadings[,PC_nos]
  loadings2 <- out2$loadings[,PC_nos]

  loadings_cor[h,] <-  diag(cor(loadings1, loadings2)) # take abs?

  pca_sign1 <- sign.pc(pcadat1,cor=T)
  pca_sign2 <- sign.pc(pcadat2,cor=T)
  
  for (p in PC_nos){
    if (pca_sign1$pval[p] < 0.05){
      loadings_sig[h,p] <- 1
    }
    if (pca_sign2$pval[p] < 0.05){
      loadings_sig[h,p+4] <- 1
    }
  }
  loadings_sig[h,loadings_sig$haf1_PC1] <- 

  print(which(pca_sign1$pval < 0.05))
  print(which(pca_sign2$pval < 0.05))
  print(loadings_cor[h,])
}

colMeans(loadings_cor, na.rm=T)
mean(loadings_sig$half1_PC1)
mean(loadings_sig$half2_PC1)

# % both halves significant PC
dim(loadings_sig[loadings_sig$half1_PC1==1 & loadings_sig$half2_PC1==1,])[1] / 1000
dim(loadings_sig[loadings_sig$half1_PC2==1 & loadings_sig$half2_PC2==1,])[1] / 1000
dim(loadings_sig[loadings_sig$half1_PC3==1 & loadings_sig$half2_PC3==1,])[1] / 1000
dim(loadings_sig[loadings_sig$half1_PC4==1 & loadings_sig$half2_PC4==1,])[1] / 1000

# % at least one half significant PC
dim(loadings_sig[loadings_sig$half1_PC1==1 | loadings_sig$half2_PC1==1,])[1] / 1000
dim(loadings_sig[loadings_sig$half1_PC2==1 | loadings_sig$half2_PC2==1,])[1] / 1000
dim(loadings_sig[loadings_sig$half1_PC3==1 | loadings_sig$half2_PC3==1,])[1] / 1000
dim(loadings_sig[loadings_sig$half1_PC4==1 | loadings_sig$half2_PC4==1,])[1] / 1000

# % exactly one half significant PC
(dim(loadings_sig[loadings_sig$half1_PC1==1 & loadings_sig$half2_PC1==0,])[1] + 
    dim(loadings_sig[loadings_sig$half1_PC1==0 & loadings_sig$half2_PC1==1,])[1]) / 1000

(dim(loadings_sig[loadings_sig$half1_PC2==1 & loadings_sig$half2_PC2==0,])[1] + 
    dim(loadings_sig[loadings_sig$half1_PC2==0 & loadings_sig$half2_PC2==1,])[1]) / 1000

(dim(loadings_sig[loadings_sig$half1_PC3==1 & loadings_sig$half2_PC3==0,])[1] + 
  dim(loadings_sig[loadings_sig$half1_PC3==0 & loadings_sig$half2_PC3==1,])[1]) / 1000

(dim(loadings_sig[loadings_sig$half1_PC4==1 & loadings_sig$half2_PC4==0,])[1] + 
    dim(loadings_sig[loadings_sig$half1_PC4==0 & loadings_sig$half2_PC4==1,])[1]) / 1000


#===============================================================================
# inspect loadings
#===============================================================================
loadings <- out$loadings 
loadings <- as.data.frame.matrix(loadings)

# http://strata.uga.edu/8370/lecturenotes/principalComponents.html
threshold <- sqrt(1/ncol(pcadat))  # cutoff for 'important' loadings

loadings[abs(loadings) < threshold] <- NA
loadings

# will use PC1, PC2, and PC3 (latter with caution)

dat_prep$PC1_all <- out$scores[,1]
dat_prep$PC2_all <- out$scores[,2]
dat_prep$PC3_all <- out$scores[,3]

dat$id <- factor(dat$id)

dat_prep$id <- factor(dat_prep$id)

dat$PC1_all <- NA
dat$PC2_all <- NA
dat$PC3_all <- NA

for (i in levels(dat_prep$id)){
  dat$PC1_all[dat$id==i] <- dat_prep$PC1_all[dat_prep$id==i]
  dat$PC2_all[dat$id==i] <- dat_prep$PC2_all[dat_prep$id==i]
  dat$PC3_all[dat$id==i] <- dat_prep$PC3_all[dat_prep$id==i]
  
}


#=============================================================
# Plot heatmaps
#=============================================================

library(corrplot)
library(RColorBrewer)
library(dplyr)

# prepare correlations
#========================
pcs_indx <- which(!is.na(match(names(dat_prep), c("PC1_all","PC2_all","PC3_all"))))
vars_indx <- which(!is.na(match(names(dat_prep), c(psych_vars, cog_vars))))

n_p <- length(pcs_indx)
n_v <- length(vars_indx)

cor_data <- dat_prep[,c(pcs_indx,vars_indx)]
M <- cor(cor_data)
M <- M[(n_p+1):(n_v+n_p), 1:n_p]

dimnames(M)[[1]] <- c("Negative affect", "Surgency","Effortful control", # ECBQ
                      "Emotion contagion","Attention to others' feelings","Prosocial actions", # EmQue
                      "Emotional symptoms","Conduct problems","Hyperactivity-inattention","Peer problems", # SDQ
                      "Prosocial behaviour", # SDQ
                      "Inattention","Hyperactivity", # ADHD
                      "Social awareness","Social cognition","Social communication",
                      "Social motivation","Repetitive behaviours", # SRS
                      "Inhibit","Shift","Emotional control", "Working memory","Planning/Organisation", # BRIEF
                      "Verbal comprehension","Visuospatial skills",
                      "Fluid reasoning","Working memory","Processing speed","Vocabulary",
                      "Nonverbal","General abilities","Cognitive proficiency") # WPPSI
dimnames(M)[[2]] <- c("PC1","PC2","PC3")

res1 <- cor.mtest(cor_data, conf.level = .95)
res1$p <- res1$p[(n_p+1):(n_v+n_p), 1:n_p]
res1$lowCI <- res1$lowCI[(n_p+1):(n_v+n_p), 1:n_p]
res1$uppCI <- res1$uppCI[(n_p+1):(n_v+n_p), 1:n_p]


# prepare loadings
#====================

loadings <- as.matrix(loadings[,1:3])
dimnames(loadings)[[1]] <- c("Negative affect", "Surgency","Effortful control", # ECBQ
                             "Emotion contagion","Attention to others' feelings","Prosocial actions", # EmQue
                             "Emotional symptoms","Conduct problems","Hyperactivity-inattention","Peer problems", # SDQ
                             "Prosocial behaviour", # SDQ
                             "Inattention","Hyperactivity", # ADHD
                             "Social awareness","Social cognition","Social communication",
                             "Social motivation","Repetitive behaviours", # SRS
                             "Inhibit","Shift","Emotional control", "Working memory","Planning/Organisation", # BRIEF
                             "Verbal comprehension","Visuospatial skills",
                             "Fluid reasoning","Working memory","Processing speed","Vocabulary",
                             "Nonverbal","General abilities","Cognitive proficiency") # WPPSI
dimnames(loadings)[[2]] <- c("PC1","PC2","PC3")



# plots
#===================================================

setwd("C:/Users/vanes/Dropbox/KCL_BMEIS/Papers/NNMF_outcomes/Figures")

col=rev(colorRampPalette(brewer.pal(n=11, name="RdBu"))(100))

# jpeg("corrplot_R.jpeg")
corrplot(M, p.mat = res1$p,method="color", col=col,
         tl.col="black",tl.cex=1, # tl.offset=0.5,
         cl.align="l",cl.ratio = 1, 
         insig = "blank",
         addgrid.col="grey")
# dev.off() 

# jpeg("corrplot_loadings_thresh.jpeg")
corrplot(loadings, is.corr=F, method="color", col=col,
         tl.col="black",tl.cex=1,  
         cl.align="l", cl.ratio = 1,
         na.label="square",na.label.col="white",
         addgrid.col="grey")
# dev.off()

