# This script is for all outdated, paused or deprecated functions.
#
# 
############## NOT MAINTAINED ##################################
# 
#
#
# // origin file: ./ARIMA.R
#

# Source the rolling origins functions 
source("./aux/roll_or.R",
       echo = F)

# Create list with several splits each to enable backtesting later.
# (180 days train in-sample), 30 days test (out-of-sample), 30 days skip).
in.sample <- 6*30
out.sample <- 30
skip <- 30

# List with differenced values | Returns a "diff_split" class object
ts_list_splits_diff <- roll_or(diff_list,in.sample,out.sample,skip)

# Based on the ACF plots more than an AR(I)MA(1,0,1) seem unrealistic, however I will still implement
# a grid search algorithm based on the AIC of the model.

# Load package gtools for permutations
library(gtools)

# Source the regain function to regain numerical data from an rsplit object
source("~/Dropbox/Uni/Masterarbeit/Master_thesis/R/aux/regain.R",
       echo = F)

# Train an Arima model based on a 180 day training period and 60 days testing
# (Here I will use the last split from the roll_functions)
# The models will then be backtested for further evaluation
fit_arima <- function(rsplit_list, horizon = 60){
        if(!(class(rsplit_list) == "diff_split")){
                stop("Supplied list must be of type 'diff_split' \n 
                     Did you use the 'roll_or' function from /aux/roll_or.R ?")
        }

        res <- as.data.frame(matrix(nrow = length(rsplit_list),ncol = length(arima_names)-4))
        orders_res <- as.data.frame(matrix(nrow = length(rsplit_list), ncol = 4))
        orders_res[,1] <- names_complete
        res[,1] <- names_complete

        # Define Error metric vector for each index.
        RMSE <- vector()
        MASE <- vector()
        MAE <- vector()

        for(i in 1:length(rsplit_list)){

                single_index <- rsplit_list[[i]]

                train <- regain(single_index[["splits"]][[length(single_index[["splits"]])]])[[1]]
                test <- regain(single_index[["splits"]][[length(single_index[["splits"]])]])[[2]]

                # Prepare for gridsearching the order space of the arima model to return the
                # one with the lowest AIC (higher orders as 4 won't be searched to spped up the process)
                # As the serial data is already stationary the (I) order will be fixed to 1
                grid_AIC = vector()
                orders <- matrix(ncol = 3)
                order_permuts <- permutations(5, 2, v = 0:4, repeats.allowed = T)

                for(k in 1:length(order_permuts[,1])){

                        grid_arima <- stats::arima(train, 
                                                   order = c(order_permuts[k,1],0,order_permuts[k,2]),
                                                   method = "ML", 
                                                   include.mean = F,
                                                   optim.method = "BFGS",
                                                   optim.control = list(maxit = 1000))

                        grid_AIC <- append(grid_AIC,grid_arima[["aic"]])
                        orders <- rbind(orders,c(order_permuts[k,1],1,order_permuts[k,2]))

                }

                order_AIC <- cbind(grid_AIC,orders[-1,])
                order_AIC <- order_AIC %>% 
                        as.data.frame() %>% 
                        filter(grid_AIC == min(grid_AIC)) 

                # Print the current function status
                print(noquote(paste0("Working on Index ",names_complete[i])))

                # Final arima with the found order
                final_arima <- stats::arima(train,
                                            order = c(order_AIC[1,2],0,order_AIC[1,4]),
                                            method = "ML", 
                                            include.mean = F,
                                            optim.method = "BFGS",
                                            optim.control = list(maxit = 1000))

                # RMSE
                res[i,2]<- accuracy(forecast(final_arima, h = horizon)[["mean"]], test)[2]

                # MASE
                res[i,3] <- accuracy(forecast(final_arima, h = horizon), test)[2,6]

                # MAE
                res[i,4] <- accuracy(forecast(final_arima, h = horizon)[["mean"]], test)[3]

                # Order
                orders_res[i,2:4] <- order_AIC[,2:4]

                #AIC
                res[i,5] <- final_arima[["aic"]]


                # Save QQ-Plots of the respective residuslas for visual inspection
                pdf(file = paste0("./img/QQ/QQ-Norm"," ",
                                  names_complete[i],
                                  ".pdf"))
                qqnorm(final_arima$residuals, 
                       main = paste0("Normal Q-Q Plot for ",names_complete[i]), 
                       ylim = c(-.2,.2), lwd = 2)
                qqline(final_arima$residuals)

                dev.off()

                # Plot the forcasts and save the files for inspection

                # Vector with test-set at its end for plotting
                test_vec <- c(rep(NA,length(train)),test)

                pdf(file = paste0("./img/arima_forecast/",
                                  names_complete[i]," forecast",
                                  ".pdf"), height = 5)

                plot(forecast(final_arima, h = horizon), 
                     main = paste0(horizon,"-days-ahead forecast for ",names_complete[i]),
                     xlim = c(floor(2/3*length(train)),length(train)+length(test)))
                lines(test_vec, col = "green4")
                abline(h = 0, v = length(train), col = "red4",lty = 5)

                dev.off()



        }

        names(orders_res) <- arima_names[c(1,6,7,8)]
        names(res) <- arima_names[1:5]

        # Final return
        final_list <- list(ErrorMetrics = res, Orders = orders_res)
        class(final_list) <- "grid_arima"
        final_list
}
test <- fit_arima(ts_list_splits_diff, horizon = out.sample)

# Separate function for the backtesting procedure. 
backtesting_ARIMA  <- function(rsplit_list, orders, horizon = 60, gridSearch = F){

        if(!(class(rsplit_list) == "diff_split")){
                stop("Supplied list must be of type 'diff_split'. \n 
                     Did you use the 'roll_or' function from /aux/roll_or.R ?")
        }
        if(!(class(orders) == "grid_arima" && gridSearch == F)){
                stop(paste("Supplied orders must originate from a 'grid_arima' class. \n 
                           Did you use the 'fit_arima' function?"))
        }

        res <- as.data.frame(matrix(nrow = length(rsplit_list), ncol = length(arima_names)-5))
        res[,1] <- names_complete

        for(i in 1:length(rsplit_list)){

                # Extract the first split object
                single_index <- rsplit_list[[i]]

                RMSE <- vector()
                MASE <- vector()
                MAE <- vector()

                # Extract the orders from the fitting process one function above
                order_vec <- c(orders$Orders[i,c(2,4)])
                print(paste("Working on Index", names_complete[i]))

                for(k in 1:length(single_index$splits)){

                        train <- regain(single_index[["splits"]][[k]])[[1]]
                        test <- regain(single_index[["splits"]][[k]])[[2]]

                        fit <- arima(train, 
                                     order = c(order_vec[[1]],0,order_vec[[2]]),
                                     method = "ML", 
                                     include.mean = F,
                                     optim.method = "BFGS",
                                     optim.control = list(maxit = 1000))


                        RMSE <- append(RMSE,accuracy(forecast(fit, h = horizon)[["mean"]], test)[2])
                        MASE <- append(MASE,accuracy(forecast(fit, h = horizon), test)[2,6])
                        MAE <- append(MAE,accuracy(forecast(fit, h = horizon)[["mean"]], test)[3])


                }

                res[i,2] <- mean(RMSE)
                res[i,3] <- mean(MASE)
                res[i,4] <- mean(MAE)

        }

        names(res) <- arima_names[1:4]
        print(paste("Done! No errors occured"))
        list(ErrorMetrics = res)
}
test2 <- backtesting_ARIMA(ts_list_splits_diff, orders = test, gridSearch = F)
#################################################################################

############### BOOTSTRAPPING #################################

# Source the bootstrapping procedure after Pascual et al. (2004)
source("./bootstrap.R", 
       echo = F)

arima.residuals <- mclapply(auto.arima.list, FUN = residuals)

# Generate 1000 bootstrap repetition for each series

bootstrap_list <- list()

for(i in 1:length(auto.arima.list)){
        
        # Print current state.
        print(noquote(paste0("Bootstrapping ",names_complete[i])))
        
        bootstrap_list[[i]] <- make_bootsrapped_series_from_arima(arima_residuals = arima.residuals[[i]], 
                                                                  series = auto.arima.list[[i]][["x"]],
                                                                  coefs = coef(auto.arima.list[[i]]),
                                                                  order = c(auto.arima.list[[i]]$arma[1],
                                                                            auto.arima.list[[i]]$arma[2]),
                                                                  n = 1000, # 1000 bootstrapped series are re-sampled
                                                                  out.sample = out.sample,
                                                                  parallel = T)
        
}
# Name the list with the bootstrapped series
names(bootstrap_list) <- names_complete

# Forecast the bootstrapped series, with parallel support. Computationally demanding. 
boot.forecast <- lapply(bootstrap_list,Boots_Predict, parallel = T, n.ahead = out.sample)


