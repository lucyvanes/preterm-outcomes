# 15/04/2021 
# Lucy Vanes
# Main analysis code for: The effects of neurobiology and environment on childhood outcomes following very preterm birth

library(jtools)

# 3.2.	Association between stimulating home environment and outcomes
#=======================================================================
setwd("C:/Users/vanes/OneDrive/Documents/GitHub/preterm-outcomes/data")
dat <- read.csv("3_Analysis_data.csv", header=T)
dat$id <- factor(dat$id)
dim(dat)

# omnibus test
anova(lm(cbind(PC1, PC2, PC3) ~ gest_weeks + sex + ses + cog_stim, dat))
# univariate tests
summary(lm(cbind(PC1, PC2, PC3) ~ gest_weeks + sex + ses + cog_stim, dat))

# plot effect
#=============
fit <- lm(PC1 ~ gest_weeks + sex + cog_stim, dat)

# tiff("PC1_scatter.tif", units="in", width=3.5, height=4, res=300)
effect_plot(fit, pred = cog_stim, interval = TRUE, plot.points = TRUE, point.size=1.5,
            point.alpha=0.4, point.color="black",
            x.label="Cognitively stimulating parenting", y.label="PC1 partial residuals",
            partial.residuals = TRUE)
# dev.off()

# 3.4.	Association between network volumes and outcomes
#=======================================================================
dat <- dat[dat$use_Jacobians==1,]

jac <- read.csv("3_NNMF15_weighted_means.csv", header=T)
jac$id <- factor(jac$id)
dim(jac)

dat <- merge(dat, jac, by="id")
dim(dat)

NNMF_networks <- paste("mean_NNMF",1:15, sep="")

# use multivariate multiple regression as omnibus test to protect from type I error
for (n in NNMF_networks){
  print("#=============================#")
  print(n)
  print("#=============================#")
  m <- lm(as.formula(paste("cbind(PC1, PC2, PC3) ~ gest_weeks + sex + pma + ses +",n)), dat) 
  print(anova(m))
}
# --> significant effect of NNMF12 (after correcting for 15 tests)
# investigate univariate effects
summary(lm(cbind(PC1, PC2, PC3) ~ gest_weeks + sex + pma + ses + mean_NNMF12, dat))

# Plot effect
#==============
fit12_2 <- lm(PC2 ~ gest_weeks + sex + pma + ses + mean_NNMF12, dat)

# tiff("PC2_scatter.tif", units="in", width=3.5, height=4, res=300)
effect_plot(fit12_2, pred = mean_NNMF12, interval = TRUE, plot.points = TRUE,
            point.size=1.5,
            point.alpha=0.4, point.color="black",
            x.label="NNMF12 volume", y.label="PC2 partial residuals", partial.residuals = TRUE)

# dev.off()
# 3.5.	Network volume, stimulating home environment, and outcomes
#====================================================================

anova(lm(cbind(PC1, PC2, PC3) ~ gest_weeks + sex + pma + ses + 
           cog_stim + mean_NNMF12, dat))

summary(lm(cbind(PC1, PC2, PC3) ~ gest_weeks + sex + pma + ses + 
             cog_stim + mean_NNMF12, dat))


# test whether adding interaction effect improves models
lm1 <- lm(PC1 ~ gest_weeks + sex + pma + ses + total_homes + mean_NNMF12, dat)
lm2 <- lm(PC1 ~ gest_weeks + sex + pma + ses + total_homes * mean_NNMF12, dat)
anova(lm1, lm2)

lm1 <- lm(PC2 ~ gest_weeks + sex + pma + ses + total_homes + mean_NNMF12, dat)
lm2 <- lm(PC2 ~ gest_weeks + sex + pma + ses + total_homes * mean_NNMF12, dat)
anova(lm1, lm2)

lm1 <- lm(PC3 ~ gest_weeks + sex + pma + total_homes + mean_NNMF12, dat)
lm2 <- lm(PC3 ~ gest_weeks + sex + pma + total_homes * mean_NNMF12, dat)
anova(lm1, lm2)

