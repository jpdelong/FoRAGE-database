###########################################################################################
### Extract Posterior Predictions of C and R, Matrix Calculations for wC and ahR
###########################################################################################

################################################# load libraries ################################################# 

library(dplyr); library(ggplot2); library(cowplot); library(brms);

### load Posteriod distributions of preda nd prey abundances

Abundance_prey_full <- as.matrix(read.csv("~Mass-Abundance Scaling/Abundance_prey_full.csv"))

Abundance_pred_full <- as.matrix(read.csv("~Mass-Abundance Scaling/Abundance_pred_full.csv"))

### load FoRAGE

forage <- read.csv('~Analysis/forage_modified_bolker_Abund_CI.csv', stringsAsFactors = TRUE)

################################################# Calculate C/R Matrix ################################################# 

C_R_predict <- Abundance_pred_full / Abundance_prey_full

################################################# Calculate wC Matrix ################################################# 

# create vector of fitted_w values

fitted_w_vec <- forage$fitted_w

length(fitted_w_vec)

#Create new matrix to represent posterior predictions of wC
# subtract 1/1e4 from pred density to account for when there is only 1 predator present
wC_predict <- (Abundance_pred_full - (1/1e4)) * fitted_w_vec 

################################################# Calculate ahR Matrix ################################################# 

#Create vector of scaled_a values

scaled_a_vec <- forage$scaled_a

# create vector of fitted_h values

fitted_h_vec <- forage$fitted_h

# create ah vector

ah_vec <- scaled_a_vec*fitted_h_vec

# create new matrix to represent posterior predictions of ahR

ahR_predict <- Abundance_prey_full * ah_vec

################################################# Calculate wC/ahR Matrix ################################################# 

wC_ahR_predict <- wC_predict/ahR_predict

################################################# Calculate Credible Intervals of Each Matrix ################################################# 

# C/R Credible intervals

C_R_ci <- t(apply(C_R_predict, 1, quantile, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE))

colnames(C_R_ci) <- c("C_R_5","C_R_25", "C_R_50", "C_R_75", "C_R_95")

# wC/ahR Credible intervals

wC_ahR_ci <- t(apply(wC_ahR_predict, 1, quantile, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE))

colnames(wC_ahR_ci) <- c("wC_ahR_5", "wC_ahR_25", "wC_ahR_50", "wC_ahR_75", "wC_ahR_95")

# Bind CIs to forage dataset

forage <- cbind(forage, C_R_ci, wC_ahR_ci)

################################################# Save FoRAGE file with CI's of Posteriors ################################################# 

write.csv(wC_ahR_predict, file = '~/Analysis/wC_ahR_predict.csv', row.names = FALSE)

write.csv(forage, file = '~/Analysis/forage_modified_bolker_Abund_CI.csv', row.names = FALSE)


