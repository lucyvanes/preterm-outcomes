# 08/01/2021
# Lucy Vanes
# PCA code for: The effects of neurobiology and environment on childhood outcomes following very preterm birth
# This code takes raw questionnaire subscales as input
# regresses out age at assessment and index of multiple deprivation
# runs PCA with permutation testing and split-half reliability analysis
# outputs PC scores
# generates heatmaps

library(caret)

setwd("C:/Users/vanes/OneDrive/Documents/GitHub/preterm-outcomes/data")
dat <- read.csv("1_PCA_data.csv", header=T)

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

#===========================================
#    Run PCA with permutation testing
#===========================================
# for each permutation, shuffle rows in each column, re-compute PCA
# function for the permutation testing taken from:
# http://bioops.info/2015/01/permutation-pca/
pcadat <- dat_prep
pcadat$id <- NULL
pcadat_scaled <- as.data.frame(scale(pcadat, center=T, scale=T))

sign.pc<-function(x,R=5000,s=10, cor=T,...){
  pc.out<-princomp(x,cor=cor,...)  # run PCA
  pve=(pc.out$sdev^2/sum(pc.out$sdev^2))[1:s] # proportion of variance for each PC
  pve.perm<-matrix(NA,ncol=s,nrow=R)
  for(i in 1:R){
    x.perm<-apply(x,2,sample)
    pc.perm.out<-princomp(x.perm,cor=cor,...)
    pve.perm[i,]=(pc.perm.out$sdev^2/sum(pc.perm.out$sdev^2))[1:s]
  }
  pval<-apply(t(pve.perm)>pve,1,sum)/R # calculate p-values
  return(list(pve=pve,pval=pval))
}

pca_sign <- sign.pc(pcadat_scaled,cor=T)
pca_sign
# PC 1-4 are significant in full sample

#===========================================
#            Original PCA
#===========================================
out <- princomp(pcadat_scaled)
summary(out)

PCs_orig <- as.data.frame(out$scores)
PC_nos <- which(pca_sign$pval < 0.05)
PC_nos

PCs_orig_sig <- PCs_orig[,PC_nos]
names(PCs_orig_sig) <- paste("orig_PC", PC_nos, sep="")

#===========================================
#        Split-half analysis
#===========================================
n_half <- round(dim(pcadat_scaled)[1] * 0.5)
H <- 1000 # how many split half replications

loadings_cor <- data.frame(matrix(nrow=H, ncol=length(PC_nos)))
loadings_sig <- data.frame(matrix(nrow=H, ncol=length(PC_nos)*2))
names(loadings_sig) <- c(paste("half1_PC",PC_nos,sep=""),paste("half2_PC",PC_nos,sep=""))

for (n in names(loadings_sig)){
  loadings_sig[n] <- 0
}

# this loop takes a while - come back in the morning
for (h in 1:H){
  print(paste("Split half ", h, sep=""), quote=F)
  print("============================")
  half1 <- sort(sample(1:dim(pcadat_scaled)[1], size=n_half))
  half2 <- setdiff(1:dim(pcadat_scaled)[1], half1)
  pcadat1 <- pcadat_scaled[half1,]
  pcadat2 <- pcadat_scaled[half2,]
  out1 <- princomp(pcadat1)                   # PCA for half 1
  out2 <- princomp(pcadat2)                   # PCA for half 2
  loadings1 <- out1$loadings[,PC_nos]
  loadings2 <- out2$loadings[,PC_nos]
  loadings_cor[h,] <-  diag(cor(loadings1, loadings2)) 

  pca_sign1 <- sign.pc(pcadat1,cor=T)         # permutation testing on half 1
  pca_sign2 <- sign.pc(pcadat2,cor=T)         # permutation testing on half 2
  for (p in PC_nos){
    if (pca_sign1$pval[p] < 0.05){
      loadings_sig[h,p] <- 1
    }
    if (pca_sign2$pval[p] < 0.05){
      loadings_sig[h,p+4] <- 1
    }
  }
  print(which(pca_sign1$pval < 0.05))         # which PCs significant in half 1
  print(which(pca_sign2$pval < 0.05))         # which PCs significant in half 2
  print(loadings_cor[h,])                     # correlation between half 1 and half 2 PC loadings
}

# Get some stats (will only make sense if full loop was run)
#===========================================================
colMeans(loadings_cor, na.rm=T)
mean(loadings_sig$half1_PC1)
mean(loadings_sig$half2_PC1)

# % both halves significant PC
dim(loadings_sig[loadings_sig$half1_PC1==1 & loadings_sig$half2_PC1==1,])[1] / H
dim(loadings_sig[loadings_sig$half1_PC2==1 & loadings_sig$half2_PC2==1,])[1] / H
dim(loadings_sig[loadings_sig$half1_PC3==1 & loadings_sig$half2_PC3==1,])[1] / H
dim(loadings_sig[loadings_sig$half1_PC4==1 & loadings_sig$half2_PC4==1,])[1] / H

# % at least one half significant PC
dim(loadings_sig[loadings_sig$half1_PC1==1 | loadings_sig$half2_PC1==1,])[1] / H
dim(loadings_sig[loadings_sig$half1_PC2==1 | loadings_sig$half2_PC2==1,])[1] / H
dim(loadings_sig[loadings_sig$half1_PC3==1 | loadings_sig$half2_PC3==1,])[1] / H
dim(loadings_sig[loadings_sig$half1_PC4==1 | loadings_sig$half2_PC4==1,])[1] / H

# % exactly one half significant PC
(dim(loadings_sig[loadings_sig$half1_PC1==1 & loadings_sig$half2_PC1==0,])[1] + 
    dim(loadings_sig[loadings_sig$half1_PC1==0 & loadings_sig$half2_PC1==1,])[1]) / H

(dim(loadings_sig[loadings_sig$half1_PC2==1 & loadings_sig$half2_PC2==0,])[1] + 
    dim(loadings_sig[loadings_sig$half1_PC2==0 & loadings_sig$half2_PC2==1,])[1]) / H

(dim(loadings_sig[loadings_sig$half1_PC3==1 & loadings_sig$half2_PC3==0,])[1] + 
  dim(loadings_sig[loadings_sig$half1_PC3==0 & loadings_sig$half2_PC3==1,])[1]) / H

(dim(loadings_sig[loadings_sig$half1_PC4==1 & loadings_sig$half2_PC4==0,])[1] + 
    dim(loadings_sig[loadings_sig$half1_PC4==0 & loadings_sig$half2_PC4==1,])[1]) / H


#=================================================
# Inspect loadings and add PC scores to dataframe
#=================================================
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

#=================================================
#             Plot heatmaps
#=================================================
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
#====================
col <- rev(colorRampPalette(brewer.pal(n=11, name="RdBu"))(100))

corrplot(M, p.mat = res1$p,method="color", col=col,
         tl.col="black",tl.cex=1, # tl.offset=0.5,
         cl.align="l",cl.ratio = 1, 
         insig = "blank",
         addgrid.col="grey")

corrplot(loadings, is.corr=F, method="color", col=col,
         tl.col="black",tl.cex=1,  
         cl.align="l", cl.ratio = 1,
         na.label="square",na.label.col="white",
         addgrid.col="grey")
