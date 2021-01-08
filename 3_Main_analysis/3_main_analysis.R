

library(jtools)

setwd("C:/Users/vanes/Dropbox/KCL_BMEIS/ePrime/data/outputs")
dat <- read.csv("ePrime_ALL_OUTCOMES_data.csv", header=T)
dat$id <- factor(dat$id)
dim(dat)

anova(lm(cbind(PC1_all, PC2_all, PC3_all) ~ gest_weeks + sex + total_homes, dat))

summary(lm(cbind(PC1_all, PC2_all, PC3_all) ~ gest_weeks + sex + total_homes, dat))


# PC1_all influenced by total_homes

fit <- lm(PC1_all ~ gest_weeks + sex + total_homes, dat)

# tiff("Figure2_beh.tif", units="in", width=3.5, height=4, res=300)
effect_plot(fit, pred = total_homes, interval = TRUE, plot.points = TRUE, point.size=1.5,
            point.alpha=0.4, point.color="black",
            x.label="Cognitively stimulating parenting", y.label="PC1 partial residuals",
            partial.residuals = TRUE)
# dev.off()


# add NNMF networks
#=======================================================================================

dat <- dat[dat$use_Jacobians==1,]

setwd("C:/Users/vanes/Dropbox/KCL_BMEIS/ePrime/data/MRI/Jacobians/NNMF/NNMF15_AAL/")
jac <- read.csv("rank15_weighted_means.csv", header=T)
jac$id <- factor(jac$id)
dim(jac)

dat <- merge(dat, jac, by="id")
dim(dat)

fit <- lm(mean_NNMF15 ~ gest_weeks + pma + sex, dat)
summary(fit)
effect_plot(fit, pred = gest_weeks, interval = TRUE, plot.points = TRUE)
effect_plot(fit, pred = pma, interval = TRUE, plot.points = TRUE)

#======================================================================================

NNMF_networks <- paste("mean_NNMF",1:15, sep="")

# use multivariate multiple regression as omnibus test to protect from type I error

for (n in NNMF_networks){
  print("#=============================#")
  print(n)
  print("#=============================#")
  
  m <- lm(as.formula(paste("cbind(PC1_all, PC2_all, PC3_all) ~ gest_weeks + sex + pma + ",n)), dat) 
  m1 <- lm(as.formula(paste("PC1_all ~ gest_weeks + sex + pma + ",n)), dat)  
  m2 <- lm(as.formula(paste("PC2_all ~ gest_weeks + sex + pma + ",n)), dat)  
  m3 <- lm(as.formula(paste("PC3_all ~ gest_weeks + sex + pma + ",n)), dat)  
  print(anova(m))
  # print(summary(m3))
}

# significant effect of NNMF12 

summary(lm(cbind(PC1_all, PC2_all, PC3_all) ~ gest_weeks + sex + pma +  mean_NNMF12, dat))

# PLOTs
#===============================================================================

fit12_2 <- lm(PC2_all ~ gest_weeks + sex + pma + mean_NNMF12, dat)

# tiff("Figure4_brain.tif", units="in", width=3.5, height=4, res=300)
effect_plot(fit12_2, pred = mean_NNMF12, interval = TRUE, plot.points = TRUE,
            point.size=1.5,
            point.alpha=0.4, point.color="black",
            x.label="NNMF12 volume", y.label="PC2 partial residuals", partial.residuals = TRUE)
# dev.off()


# FULL MODEL
#===============================================================================

anova(lm(cbind(PC1_all, PC2_all, PC3_all) ~ gest_weeks + sex + pma + 
           total_homes + mean_NNMF12, dat))

summary(lm(cbind(PC1_all, PC2_all, PC3_all) ~ gest_weeks + sex + pma + 
             total_homes + mean_NNMF12, dat))



lm1 <- lm(PC1_all ~ gest_weeks + sex + pma + total_homes, dat)
lm2 <- lm(PC1_all ~ gest_weeks + sex + pma + total_homes + mean_NNMF12, dat)
anova(lm1, lm2)

lm1 <- lm(PC2_all ~ gest_weeks + sex + pma + total_homes, dat)
lm2 <- lm(PC2_all ~ gest_weeks + sex + pma + total_homes + mean_NNMF12, dat)
anova(lm1, lm2)

lm1 <- lm(PC2_all ~ gest_weeks + sex + pma + mean_NNMF12, dat)
lm2 <- lm(PC2_all ~ gest_weeks + sex + pma + mean_NNMF12 + total_homes, dat)
anova(lm1, lm2)

lm1 <- lm(PC3_all ~ gest_weeks + sex + pma + total_homes, dat)
lm2 <- lm(PC3_all ~ gest_weeks + sex + pma + total_homes + mean_NNMF12, dat)
anova(lm1, lm2)


# Test addition of variables
#===============================
for (n in NNMF_networks){
  print("#=============================#")
  print(n)
  print("#=============================#")
  
  m0 <- lm(PC1_all ~ gest_weeks + sex + pma + total_homes, dat) 
  m1 <- lm(as.formula(paste("PC1_all ~ gest_weeks + sex + pma + total_homes + ",n)), dat) 
 
  print(anova(m0, m1))
}

m0 <- lm(PC2_all ~ gest_weeks + sex + pma + mean_NNMF12, dat) 
m1 <- lm(PC2_all ~ gest_weeks + sex + pma + mean_NNMF12 + total_homes, dat) 
anova(m0, m1)





