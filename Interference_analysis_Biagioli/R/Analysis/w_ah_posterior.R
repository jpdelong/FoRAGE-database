###########################################################################################
### Extract Posterior Predictions of a, h, w Matrix Calculations for w/ah
###########################################################################################

################################################# load libraries ################################################# 

library(dplyr); library(ggplot2); library(cowplot); library(brms);

### load Posteriod distributions of w, a and h

w_posterior <- as.matrix(read.csv("~/Julia/posterior_w_matrix.csv"))

h_posterior <- as.matrix(read.csv("~/Julia/posterior_h_matrix.csv"))

a_posterior <- as.matrix(read.csv("~/Julia/posterior_scaled_a_matrix.csv"))


### Calculate ah

ah_posterior <- a_posterior * h_posterior

## Calculate w/ah

w_ah_posterior <- w_posterior / ah_posterior

#Take a subsample of the w/ah posterior 
#set.seed(123)

#n_draws_samp <- sample(seq_len(ncol(w_ah_posterior)), 20000)

#w_ah_posterior_sample <- w_ah_posterior[, n_draws_samp]

#w_ah_posterior_sample <- cbind(Inter_ID = seq_len(nrow(w_ah_posterior)),
                               #as.data.frame(w_ah_posterior_sample))

#write_csv(w_ah_posterior_sample, "/Users/francisbiagioli/Documents/R/Interference/Bolker_Abund_CI/w_ah_post_samp.csv")

w_ah_posterior <- cbind(Inter_ID = seq_len(nrow(w_ah_posterior)),
                        as.data.frame(w_ah_posterior))

write_csv(w_ah_posterior, "~/w_ah_post.csv")


