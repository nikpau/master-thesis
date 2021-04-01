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
install_keras()

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
                    "./ARIMA.R")

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
              ic = "bic"
)

# Fitting Holt's linear trend method to the data
holt <- lapply(

               train_test_ets$train,
               forecast::holt,
               h = OUT_SAMPLE,
               ic = "bic"
)

# Error Metrics
ses_error <- error_metrics_ets(ses, train_test_ets$train, train_test_ets$test)
holt_error <- error_metrics_ets(holt, train_test_ets$train, train_test_ets$test)

# Naming the error metrics
ses_error <- bind_cols(

                       Index = names_complete,
                       ses_error
)

holt_error <- bind_cols(

                        Index = names_complete,
                        holt_error
)

# Information criteria for ets fit
ses_info_crits <- info_crits_ets(ses)
holt_info_crits <- info_crits_ets(holt)

# Naming the information criteria
ses_info_crits <- bind_cols(Index = names_complete, ses_info_crits)
holt_info_crits <- bind_cols(Index = names_complete, holt_info_crits)

# Plotting
plot_save_forecast(ses, train_test_ets$test, 20, "./img/exp_sm/ses")
plot_save_forecast(holt, train_test_ets$test, 20, "./img/exp_sm/holt")

######################################################################
#                         ARIMA 
######################################################################
 b = lapply(ts_list,log)
# Training and testing sets for ARIMA estimation
train_test_arima <- make_training_and_testing_sets(

                                                   b, 
                                                   OUT_SAMPLE
)

# Estimation of up to ARMA(7,7)
auto_arima_list_aicc <- list_Auto_Arima(

                                   b, 
                                   crit = "aicc",
                                   parallel = T, 
                                   n.cores = mc_cores, 
                                   out.sample = OUT_SAMPLE

)

# Extract the orders of the fitted arima models
arma_orders_bic <- get_orders_from_arima_fit(auto_arima_list)
arma_orders_aicc <- get_orders_from_arima_fit(auto_arima_list_aicc)

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
                   "./img/arima_forecast"
)
################ ANNs ##############################################
