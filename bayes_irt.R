# Bayesian analysis of open big5 data (Extraversion factor only)
library(rstan)
library(ggplot2)
library(tidyr)
library(dplyr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
setwd('~/Dropbox/R/srug_preso/')

# Data processing
# ---------------

# Get data
temp <- tempfile()
download.file('http://openpsychometrics.org/_rawdata/BIG5.zip',temp)
data <- read.csv(unz(temp, "BIG5/data.csv"), sep = '\t', stringsAsFactors = FALSE)
unlink(temp)

# Item columns (Extraversion only)
item_names <- paste0('E', 1:10)

# Summarize ratings
summary(data[, item_names])

# Drop 1 row with missing data
lapply(data[, item_names], function(x) which(x == 0))
data <- data[-19065,]

# English natives only
data <- data[data$engnat == 1, ]

# Items to reverse
# https://ipip.ori.org/newBigFive5broadKey.htm
items_reversed <- c('E2', 'E4', 'E6', 'E8', 'E10')

# Function to reverse item responses
reverse_answer <- function(ans, max_value=5, min_value=1) {
  (max_value + min_value) - ans
}

data[, items_reversed] <- lapply(data[, items_reversed], reverse_answer)

# Exploration
# -----------

# Look at distributions
data_long <- data[, item_names] %>%
  gather(key = 'question', value = 'response')

ggplot(data_long, aes(response)) +
  geom_histogram() +
  facet_wrap(~question) +
  labs(x = 'Response', y = 'Frequency')

ggsave('response_distributions.png', width=6, height=5)

# Classic analysis
# ----------------

# Item correlations
cor(data[, item_names])
image(cor(data[, item_names]))

# Item-total correlations
for (i in item_names) {
  other_items <- setdiff(item_names, i)
  totals <- rowSums(data[, other_items])
  cat(i, cor(data[, i], totals), '\n')
}

# Bayesian Modelling
# ------------------

# Specify data for Stan
ex_data <- list(
  K = 5,
  N = nrow(data),
  I = 10,
  X = data[, item_names]
)

# Load stan model Luo, Jiao 2018
ex_model <- stan_model('graded_response_model.stan', 'extraversion')

# Variational bayes to approximate posterior
ex_vb <- vb(ex_model, ex_data, seed=56327)

# Look at summary
summary(ex_vb, pars = c('alpha', 'kappa'), probs = c(0.25, 0.5, 0.75))

# Forest plot
plot(ex_vb, pars = 'alpha') 

# Extract posterior draws
trace <- rstan::extract(ex_vb)

alpha_means <- colMeans(trace$alpha)
alpha_lower <- apply(trace$alpha, 2, function(x) quantile(x, probs = 0.05))
alpha_upper <- apply(trace$alpha, 2, function(x) quantile(x, probs = 0.95))

kappa_means <- apply(trace$kappa, c(2,3), mean)
kappa_lower <- apply(trace$kappa, c(2,3), function(x) quantile(x, probs = 0.05))
kappa_upper <- apply(trace$kappa, c(2,3), function(x) quantile(x, probs = 0.95))

# Plot item response functions
# ----------------------------

# Set of thetas for x-axis
theta_space <- seq(-4, 4, 0.1)

# Sigmoid function
sigmoid <- function(z) {
  1 / (1 + exp(-z))
}

# Item Response Function
category_boundary <- function(thetas, alpha, kappa) {
  sigmoid(alpha * (thetas - kappa))
}

# Plot example of 2PL model
ggplot(data.frame(x=theta_space, y=category_boundary(theta_space, 1, 0)),
       aes(x, y)) + geom_line() + 
  labs(title = 'Example of binary question',
       x = 'theta', y = 'P(X = x)')

ggsave('first_example_2pl.png', width=3, height=3, units = 'in')

# Create plots for all response category boundaries using estimated parameters
plots <- list()
for (i in 1:10) {
  # Prepare data
  plot_list <- list() 
  for(k in 1:length(kappa_means[i,])) {
    boundary_name <- paste0('k=', k)
    plot_list[[boundary_name]] <- data.frame(
      x = theta_space,
      y_mean = category_boundary(theta_space, alpha_means[i], kappa_means[i, k]),
      y_lower = category_boundary(theta_space, alpha_lower[i], kappa_lower[i, k]),
      y_upper = category_boundary(theta_space, alpha_upper[i], kappa_upper[i, k])
    )
  }
  
  plot_df <- do.call(rbind, plot_list)
  plot_df$boundary <- substr(rownames(plot_df), 1, 3)
  
  plots[[paste0('E', i)]] <- ggplot(plot_df, aes(x, y_mean, colour = boundary, group = boundary)) + 
    geom_ribbon(aes(ymin = y_lower, ymax = y_upper, fill = boundary), colour = NA, alpha = 0.2) +
    geom_line() +
    labs(title = paste0('Question E', i), x = 'theta', y = 'P(X > k)')
}

# Plot E7, to be used as an example
plots$E7
alpha_means[7]
kappa_means[7,]
ggsave('second_example_grm.png', width=4, height=3, units='in')


# Plot E1 as example of a harder question
plots$E1
alpha_means[1]
kappa_means[1,]
ggsave('hard_question.png', width=4, height=3)

# Plot E6 as example of an easier question
plots$E6
alpha_means[6]
kappa_means[6,]
ggsave('easy_question.png', width=4, height=3)

# Bonus: MCMC sampling
# --------------------

# Subset data for shorter runtime 
N <- 500

ex_data_subset <- list(
  K = 5,
  N = N,
  I = 10,
  X = data[sample(1:nrow(data), N), item_names]
)

# Sample from posterior
ex_trace_subset <- sampling(
  ex_model, 
  data = ex_data_subset,
  chains = 4,
  iter = 1000,
  algorithm = 'NUTS',
  seed = 536278)

# Summarize posterior
summary(ex_trace_subset, pars = c('alpha', 'kappa'), probs = c(0.25, 0.5, 0.75))

# Forest plot
plot(ex_trace_subset, pars = 'alpha')

# Diagnostic plots
stan_diag(ex_trace_subset, information = 'sample')
stan_rhat(ex_trace_subset)
stan_ess(ex_trace_subset)
stan_mcse(ex_trace_subset)

# Trace plots
trace_subset <- rstan::extract(ex_trace_subset)
par(mfrow=c(5, 2), mar=c(2,2,1,1))
for (i in 1:10) {
  plot(trace_subset$alpha[,i], type='l')
}
par(mfrow=c(1, 1), mar=c(5,4,4,2))

# Posterior distributions
par(ask=TRUE)
for (i in 1:10) {
  hist(trace_subset$alpha[, i], main = i)
}
par(ask=FALSE)
