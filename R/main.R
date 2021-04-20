# ################### INFORMATION ##################################
#
# This is the main script for my time series forecast comparison for different
# parametric statistical models versus neural networks.
#
# All of the relevant work will be imported here so that the estimation process
# becomes comprehensible. Most of the code throughout the helper functions is
# commented to get an idea of the process and structure of my work.
#
#
#
#
############### PRE REQUISITS ###################################

# Load relevant packages and install them if not already done.
list_of_packages <- c("moments",
                      "tidyverse",
                      "lubridate",
                      "zoo",
                      "stats",
                      "aTSA",
                      "sn",
                      "forecast",
                      "parallel",
                      "foreach",
                      "keras",
                      "tfestimators",
                      "tfruns",
                      "rsample",
                      "recipes"
                      )

# Create a list of all missing packages
new_packages <- list_of_packages[!(list_of_packages %in%
                                   installed.packages()[, "Package"])]

# Install all missing packages at one
if (length(new_packages)) {
        install.packages(new_packages)
}
# Recursive load of packages
invisible(

          lapply(list_of_packages, library, character.only = T)
)

#### KERAS INFORMATION ########
#
# If you have not yet installed TensorFlow or Keras on your machine
# you can do so by executing the install_keras() function below.
# For details see the documentation at:
#                   'https://keras.rstudio.com/reference/install_keras.html'

# install_keras()

# Detect the number of cores on your system to enable multi-core processing.
# Last core is omitted to prevent system freeze during exec.
mc_cores <- detectCores() - 1

############ ONLY FOR R-STUDIO USERS ##################
#
# If you use RStudio you need to set the working dir to the file path
# to enable relative paths.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

################ SOURCING ################

# Create new environment to source the helper functions into.
# This avoids the .GlobalEnv to get too crowded.
funcEnv <- new.env()

# Vector with paths to the files to be sourced.
script_sources <- c("./import_pre-transforms.R",
                    "./descriptives.R",
                    "./auxf/helper.R",
                    "./exp_smoothing.R",
                    "./auxf/plotting.R",
                    "./ARIMA.R",
                    "./neural-nets.R")

# Source the functions and variables into the funcEnv
invisible(lapply(script_sources, source, local = funcEnv))

# Attach the environment to the search path
attach(funcEnv)

# Tidy up
rm(funcEnv, list_of_packages, script_sources)

######################################################################
#                       DESCRIPTIVE STATISTICS
######################################################################

# Plot the raw series
plot_Raw_Series()

# Plot ACFs and PACFs || can be found in ./img/ACFs
save_ACF(diff_list, n.lag = 90, type = "acf")
save_ACF(diff_list, n.lag = 90, type = "pacf")

# Tests for stationarity

# Augmented Dickey Fuller test
adf_test <- list_ADF(diff_list)

# KPSS test
kpss_test <- list_KPSS(diff_list)

# Tests for normality with both the raw data and the differenced one
jb_test_raw <- list_JB(ts_list)
jb_test_diff <- list_JB(lapply(ts_list, diff))

#######################################################################
#                           ESTIMATIONS 
#######################################################################

# Global forecast horizon
OUT_SAMPLE <- 7

########################################################################
#                       EXPONENTIAL SMOOTHING
########################################################################

# Create training and testing sets
train_test_ets <- make_training_and_testing_sets(ts_list, out.sample = OUT_SAMPLE)

# Fitting an Simple Exponential Smoothing model to the data
ses <- lapply(

              train_test_ets$train,
              forecast::ses,
              h = OUT_SAMPLE,
              ic = "aicc"
)

# Fitting Holt's linear trend method to the data
holt <- lapply(

               train_test_ets$train,
               forecast::holt,
               h = OUT_SAMPLE,
               ic = "aicc"
)

# Fitting Holt's linear damped-trend method to the data
damped_holt <- lapply(
  
  train_test_ets$train,
  forecast::holt,
  h = OUT_SAMPLE,
  ic = "aicc",
  damped = T
)

# Write ses parameters into data frame
ses_pars <- as.data.frame(matrix(ncol = 2, nrow = length(names_complete)))
ses_pars[,1] <- names_complete
names(ses_pars) <- c("Index", "alpha")
for (i in seq_len(length(names_complete))) {
  ses_pars[i,2] <- ses[[i]]$model$par[1]
}

# Write holt parameters into data frame
holt_pars <- as.data.frame(matrix(ncol = 3, nrow = length(names_complete)))
holt_pars[,1] <- names_complete
names(holt_pars) <- c("Index", "alpha", "gamma")
for (i in seq_len(length(names_complete))) {
  holt_pars[i,2] <- holt[[i]]$model$par[1]
  holt_pars[i,3] <- holt[[i]]$model$par[2]
}

# Write damped  holt parameters into data frame
damped_holt_pars <- as.data.frame(matrix(ncol = 4, nrow = length(names_complete)))
damped_holt_pars[,1] <- names_complete
names(damped_holt_pars) <- c("Index", "alpha", "gamma","phi")
for (i in seq_len(length(names_complete))) {
  damped_holt_pars[i,2] <- damped_holt[[i]]$model$par[1]
  damped_holt_pars[i,3] <- damped_holt[[i]]$model$par[2]
  damped_holt_pars[i,4] <- damped_holt[[i]]$model$par[3]
}

# Error Metrics
ses_error <- error_metrics_ets(ses, train_test_ets$train, train_test_ets$test)
holt_error <- error_metrics_ets(holt, train_test_ets$train, train_test_ets$test)
damped_holt_error <- error_metrics_ets(damped_holt, train_test_ets$train, train_test_ets$test)

# Naming the error metrics
ses_error <- bind_cols(

                       Index = names_complete,
                       ses_error
)

holt_error <- bind_cols(

                        Index = names_complete,
                        holt_error
)

damped_holt_error <- bind_cols(
  
  Index = names_complete,
  damped_holt_error
)

# Information criteria for ets fit
ses_info_crits <- info_crits_ets(ses)
holt_info_crits <- info_crits_ets(holt)
damped_holt_info_crits <- info_crits_ets(damped_holt)

# Naming the information criteria
ses_info_crits <- bind_cols(Index = names_complete, ses_info_crits)
holt_info_crits <- bind_cols(Index = names_complete, holt_info_crits)
damped_holt_info_crits <- bind_cols(Index = names_complete, damped_holt_info_crits)

# Ljung-Box test for the ses residuals
ses_resid <- lapply(ses, residuals)
ses_ljungBox <- lapply(ses_resid, Box.test, type = "L", lag = 2)
ses_ljungBox_summary <- data.frame(matrix(ncol = 3, nrow = length(names_complete)))
ses_ljungBox_summary[,1] <- names_complete
for (i in seq_len(length(names_complete))) {
  ses_ljungBox_summary[i,2] <- ses_ljungBox[[i]]$statistic
  ses_ljungBox_summary[i,3] <- ses_ljungBox[[i]]$p.value
}
names(ses_ljungBox_summary) <- c("Index", "Test statistic", "p.value")

# Plotting
plot_save_forecast(ses, train_test_ets$test, 20, "./img/exp_sm/ses")
plot_save_forecast(holt, train_test_ets$test, 20, "./img/exp_sm/holt")
plot_save_forecast(damped_holt, train_test_ets$test, 20, "./img/exp_sm/damped_holt")



######################################################################
#                         ARIMA assuming Normal
######################################################################

# Log_ transform for the indices
log_indices = lapply(ts_list,log)

# Training and testing sets for ARIMA estimation
train_test_arima <- make_training_and_testing_sets(

                                                   log_indices, 
                                                   OUT_SAMPLE
)

# Estimation of up to ARMA(14,14). Computationally intensive!
#
# Fitted object can directly be loaded via 
# auto_arima_list <- readRDS("./_objects/arima_fit.rds")
auto_arima_list <- list_Auto_Arima(

                                   log_indices, 
                                   crit = "aicc",
                                   parallel = T, 
                                   n.cores = mc_cores, 
                                   out.sample = OUT_SAMPLE

)

# Ljung-Box test for the arima residuals
arima_resid <- lapply(auto_arima_list, residuals)
arima_ljungBox <- lapply(arima_resid, Box.test, type = "L", lag = 1)
arima_ljungBox_summary <- data.frame(matrix(ncol = 3, nrow = length(names_complete)))
arima_ljungBox_summary[,1] <- names_complete
for (i in seq_len(length(names_complete))) {
  arima_ljungBox_summary[i,2] <- arima_ljungBox[[i]]$statistic
  arima_ljungBox_summary[i,3] <- arima_ljungBox[[i]]$p.value
}
names(arima_ljungBox_summary) <- c("Index", "Test statistic", "p.value")

# Extract the orders of the fitted arima models
arma_orders_aicc <- get_orders_from_arima_fit(auto_arima_list)

# Extract information criteria from the arima fits
arima_info_crits <- info_crits_arima(auto_arima_list)

# Forecast the ARIMA models based on their respective models
arima_forecast <- mclapply(

                           auto_arima_list, forecast, 
                           h = OUT_SAMPLE, 
                           mc.cores = mc_cores
)

# Error metrics 
arima_error <- error_metrics_arima(

                                   arima_forecast,
                                   train_test_arima$test,
                                   train_test_arima$train
)

# Naming the forecast table for better intuition 
arima_error <- bind_cols(

                         Index = names_complete, 
                         arima_error
)

# Plotting the forecasts and saving them
plot_save_forecast(

                   arima_forecast,
                   train_test_arima$test,
                   20,
                   "./img/arima_forecast/norm"
)

############################################################
#                 Arima assuming Skewed student t
#############################################################
cl <- makeCluster(mc_cores)

# Create training and testing sets
train_test_arima_std <- make_training_and_testing_sets(diff_list, OUT_SAMPLE)

arima_sstd <- lapply(train_test_arima_std$train, autoarfima, ar.max = 5, ma.max = 5, 
                     method = "full", distribution.model = "sstd", include.mean = F,
                     criterion = "BIC", cluster = cl)

# Extract residuals
sstd_resid <- lapply(arima_sstd, function(x) return(x$fit@fit$residuals))

# Ljung-Box test for the arima_sstd residuals
sstd_ljungBox <- lapply(sstd_resid, Box.test, type = "L", lag = 1)
sstd_ljungBox_summary <- data.frame(matrix(ncol = 3, nrow = length(names_complete)))
sstd_ljungBox_summary[,1] <- names_complete
for (i in seq_len(length(names_complete))) {
  sstd_ljungBox_summary[i,2] <- sstd_ljungBox[[i]]$statistic
  sstd_ljungBox_summary[i,3] <- sstd_ljungBox[[i]]$p.value
}
names(sstd_ljungBox_summary) <- c("Index", "Test statistic", "p.value")

# Forecast skew t arima
arima_sstd_forecast <-lapply(arima_sstd, 
                             function(x) return(arfimaforecast(x[[1]], 
                                                               n.ahead = OUT_SAMPLE)))

# Extract the forecasts, stitch them back to the series, and un-difference
sstd_retransformed <- sstd_retransform(arima_sstd_forecast, train_test_arima_std$train, const_int)

# Extract only the forecast values
sstd_plain_forecast <- lapply(sstd_retransformed, tail, n = OUT_SAMPLE)

# Write error metrics into data frame
arima_sstd_error <- error_metrics_arima(
  
  forecasts = sstd_plain_forecast,
  make_training_and_testing_sets(ts_list, OUT_SAMPLE)$test,
  make_training_and_testing_sets(ts_list, OUT_SAMPLE)$train
  
)

arima_sstd_error <- bind_cols(
  
  Index = names_complete,
  arima_sstd_error
  
)

# Re-create test-set for plotting
testSet_sstd <- make_training_and_testing_sets(ts_list, OUT_SAMPLE)$test

# Plotting and saving
plot_save_forecast(
  
  seriesWithForc = sstd_retransformed,
  testset = testSet_sstd,
  window = 20,
  path = "./img/arima_forecast/sstd",
  n.ahead = OUT_SAMPLE
)



##############################################################
#                ANNs
#############################################################

# Split the ts_list into two parts because of memory usage (16Gb on my machine 
# is not enough to train in one batch)
ts_list_1_7 <- list()
ts_list_8_14 <- list()
for (i in  seq_len(length(ts_list) / 2)) {
  ts_list_1_7[[i]] <- ts_list[[i]]
  ts_list_8_14[[i]] <- ts_list[[i + 7]]
}

# Create Vector of lags to train on the networks
lags_lstm <- c(2,4,7,14)
lags_mlp <- c(2,4,7,10,14)

# Create supervised learning problems (split the series into a labeled train 
# and test set), for all the lags from the lag vector
nn.gridsearch_1_7 <- lag_list_nn(ts_list_1_7, lags, OUT_SAMPLE)
nn.gridsearch_8_14 <- lag_list_nn(ts_list_8_14, lags, OUT_SAMPLE)

flags.mlp <-  list(dropout = c(.1),
                   dense.units = c(16, 32, 64, 128),
                   val.split = c(.1, .2),
                   neg.slope = c(0, 0.02))

flags.lstm <- list(dropout = c(0.1),
                   lstm.units1 = c(16,32,64),
                   lstm.units2 = c(16,32,64),
                   val.split = c(0.2),
                   neg.slope = c(0,0.02))

########## Actual Training #########################################

# This is the actual training process using the 'super_grid_search'
# function located in ./neural-nets.R
#
# Training takes around 4-6 hours for the MLP and 8-16 hours for the LSTM.
# If you run this yourself, please be aware that due to the nature of the 
# stochastic gradient descent that is used, your results are likely to be 
# different from the ones I saved from my training

# mlp_grid_search <- super_grid_search(nn.gridsearch,
#                                      flags = flags.mlp,
#                                      "./schedules/mlp_train.R")
# 
# lstm_grid_search_1_7 <- super_grid_search(nn.gridsearch_1_7,
#                                           flags = flags.lstm,
#                                           "./schedules/lstm_train.R")
# 
# lstm_grid_search_8_14 <- super_grid_search(nn.gridsearch_8_14,
#                                            flags = flags.lstm,
#                                            "./schedules/lstm_train.R")

#########################################################################

# Load in the finished training runs with all the hyperparameters,
# architectures and loss values.

#MLP
mlp_runs <- readRDS("./_objects/mlp_fit.rds")

# LSTM
lstm1 <- readRDS("./_objects/lstm_1_7_fit.rds")
lstm2 <- readRDS("./_objects/lstm_8_14_fit.rds")

lstm_runs <- c(lstm1, lstm2)
names(lstm_runs) <- names_complete
rm(lstm1,lstm2)

# Add a column with the lag at which training has been conducted
mlp_runs <- add_lag(mlp_runs, lags_mlp)
lstm_runs <- add_lag(lstm_runs, lags_lstm)

# Extract the best run from the tuning runs
mlp_best_runs <- get_best_runs(mlp_runs)
lstm_best_runs <- get_best_runs(lstm_runs)

# Extract the hyperparameters that yielded the best runs
mlp_best_flags <- extract_flags(mlp_best_runs)
lstm_best_flags <- extract_flags(lstm_best_runs)

# Train the final networks with the extracted flags
mlp_final <- train_and_predict(ts_list, mlp_best_flags, OUT_SAMPLE, "mlp")
lstm_final <- train_and_predict(ts_list, lstm_best_flags, OUT_SAMPLE, "lstm")

# Get error metrics for the neural network trainings
mlp_error <- error_metrics_nn(mlp_final)
lstm_error <- error_metrics_nn(lstm_final)

# Naming the errors for a better overview
mlp_error <- bind_cols(
  
  Index = names_complete, 
  mlp_error
)

lstm_error <- bind_cols(
  
  Index = names_complete, 
  lstm_error
)

# Plot the results of the neural network prediction
plot_save_nn(mlp_final, 20, "mlp")
plot_save_nn(lstm_final, 20, "lstm")
