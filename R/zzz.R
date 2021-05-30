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

# In this script I implement the bootstrapping process from Pascual et al.(2004). 
# 
# Bootstrapped samples are generated by recalculating the series with randomly drawn 
# residuals with replacement from the empirical distribution function of the residuals.
#
# For the forecasting process from any ARMA(p,q) model, the last p obs. are fixed and the predictions 
# are calculated for every bootstrapped series.

# Multicore support for the bootstrapping process.
library(parallel)

# Rescale original residuals from fit after Stine (1987).
rescale_arima_residuals <- function(arima_residuals,order,series){

        p <- order[1]
        scale_factor <- sqrt((length(series) - p) / (length(series) - 2 * p))
        arima_residuals <- arima_residuals * scale_factor
        arima_residuals

}

# Prepare resampling of an ARIMA Series, based on coefs. Only for orders larger than (0,0,0).
calc_Series <- function(ar_coefs,ma_coefs,arima_residuals,series = NA, is.forecast = F,n.ahead = 30,init_vector = NULL){

        if(is.forecast == T){

                if(length(init_vector) != 7)
                        stop("Vector of starting values must be 7 digits long (for now).")

                end <- n.ahead
                y <- append(init_vector,rep(NA,n.ahead))
                new_series <- y

        }
        else{
                # Coerce series to double, and init intermediate series vector.
                series <- as.double(series)
                y <- vector()

                # Padding for max AR(7) is hard coded!
                # I use a random draw from the series added to a random residuum to initialize the resampling.
                y[1:7] <- sample(series,7) + sample(arima_residuals,7)
                y <- append(y,series)

                # For Bootstrapping the series, the whole length is used.
                end <- length(series)

                # Bootstrapped series will be created recursively.
                new_series <- vector()
        }

        for(i in 8:(end + 7)){

                # Naïve, isn't it? 
                # To be honest I had no better idea for implementing a 
                # recursive time series generator based on ARMA coefs.
                new_series[i] <- ar_coefs[1]*y[i-1] +
                        ar_coefs[2]*y[i-2] +
                        ar_coefs[3]*y[i-3] +
                        ar_coefs[4]*y[i-4] +
                        ar_coefs[5]*y[i-5] +
                        ar_coefs[6]*y[i-6] +
                        ar_coefs[7]*y[i-7] +
                        sample(arima_residuals,1)+
                        ma_coefs[1]* sample(arima_residuals,1)+
                        ma_coefs[2]* sample(arima_residuals,1)+
                        ma_coefs[3]* sample(arima_residuals,1)+
                        ma_coefs[4]* sample(arima_residuals,1)+
                        ma_coefs[5]* sample(arima_residuals,1)+
                        ma_coefs[6]* sample(arima_residuals,1)+
                        ma_coefs[7]* sample(arima_residuals,1)

                # Update vector with forecast value.
                if(is.forecast == T){y <- new_series}

        }

        #Remove padding and return series.
        return(new_series[-(1:7)])

}

# Derive Coefficients from list.
make_ar_ma_coefs <- function(coefs,order){


        # Initialize empty order vectors
        ar_coefs <- ma_coefs <- rep(0,7)

        # ARMA(0,0)
        if(order[1] == 0 && order[2] == 0){
                return(list(ar_coefs = ar_coefs, ma_coefs = ma_coefs))
        }
        # AMRA(0,q)
        else if(order[2] == 0){
                ar_coefs <- coefs
                ar_coefs <- ar_coefs %>% 
                        append(rep(0,7)) %>% 
                        head(7)

        }
        # ARMA(p,0)
        else if(order[1] == 0){
                ma_coefs <- coefs
                ma_coefs <- ma_coefs %>% 
                        append(rep(0,7)) %>% 
                        head(7)
        }
        else{
                # Fill vectors according to order
                ar_coefs[1:order[1]] <- head(coefs,order[1])
                ma_coefs[1:order[2]] <- tail(coefs,order[2])
        }

        list(ar_coefs = ar_coefs, ma_coefs = ma_coefs)
}

# Make recursive col names
recursive_col_Names <- function(n, String){

        String <- as.character(String)
        series_Names <- rep(NA,n)
        for(i in 1:n){
                series_Names[i] <- paste0(String," ", i)
        }
        series_Names
}

# Create Bootstrap replicates of an ARIMA series. I only consider ARIMA(p,1,q) models
# as all my seris are stationary already, and thus dont need diffrencing.
bootsrapp_from_arima <- function(arima_residuals, series, coefs, order = c(0,0), n = 100, 
                                               out.sample = 30, parallel = F){

        # Rescale arima_residuals
        arima_residuals <- rescale_arima_residuals(arima_residuals,order,series)

        # Check whether residuals and series are of smae length.
        if(length(arima_residuals) != length(series))
                stop("Residuals must have same length as series")

        # Check that only arma order are supplied
        if(length(order) != 2)
                stop("Please only provide ARMA(p,q) orders")

        # Init result vector
        res <- matrix(ncol = length(series),nrow = n)

        # Intercept only -> ARIMA(0,0,0)
        if(sum(order) == 0){

                # Set coefs to 0 to just sample with a random draw from the error term
                coefs <- 0
                # This is just c + ê_t (intercept + resampled residuum)
                res <- replicate(n,sample(arima_residuals,length(series)))

        }

        else if(sum(order) != 0){

                # Other order exepct (0,0,0)
                arma  <- make_ar_ma_coefs(coefs,order)

                #Multicore support
                if(parallel == T){

                        cl <- makeCluster(detectCores()-1)

                        # Get library support needed to run the code
                        clusterEvalQ(cl,library(MASS))

                        clusterExport(cl, c("arma","arima_residuals","series"), envir = environment())
                        clusterExport(cl,"calc_Series", envir = .GlobalEnv)

                        res <- parSapply(cl,1:n, 
                                         function(i, ...){
                                                 x <- calc_Series(arma[[1]],arma[[2]],arima_residuals, series = series)
                                         }
                        )

                        stopCluster(cl)

                }

                else if (parallel == F){

                        # Resample Series
                        res <- replicate(n,calc_Series(arma[[1]],arma[[2]],arima_residuals, series = series))

                }

        }
        # Name the cols
        series_names <- recursive_col_Names(n,"Series")

        #Name the orders
        names(order) <- c("AR(p)","MA(q)")

        # Return Bootstrapped series and the orders for re-checking
        final <- list(bootSeries =  res, 
                      orders = order, 
                      coefs = coefs, 
                      origSeries = series, 
                      residuals = arima_residuals)
        class(final) <- "bootstrap_Arima"
        return(final)
}

# NOT RUN:
# f <- make_bootsrapped_series_from_arima(arima_residuals = auto.arima.list$CAC40$residuals,
# 										series = diff_list[[1]],
# 										coefs = coef(auto.arima.list$CAC40),
# 										order = c(auto.arima.list$CAC40$arma[1], auto.arima.list$CAC40$arma[2]), parallel = T)

# Calculate n days ahead from the bootstrapped series.
Boots_Predict <- function(bootstraps, n.ahead = 30, parallel = F){

        if(!(class(bootstraps) == "bootstrap_Arima")){

                stop("Supplied input is not of class 'boots_Arima'.\n
                     Did you use the 'make_bootsrapped_series_from_arima' function? ")

        }

        # Initialize result vector
        raw <- as.data.frame(matrix(ncol = n.ahead, nrow = ncol(bootstraps$bootSeries)))

        # Naming for result vector
        names(raw) <- recursive_col_Names(n.ahead,"T +")

        # Disassemble class 
        arima_residuals <- bootstraps$residuals
        bootMatrix <- bootstraps$bootSeries
        series <- c(bootstraps$origSeries)

        # Open a multi-core cluster for parallel computing.
        if(parallel == T){

                print(noquote("Init..."))
                cl <- makeCluster(detectCores())
                clusterExport(cl,c("bootMatrix","bootstraps"),envir = environment())

                print(noquote("Calculating..."))

                res <- parApply(cl = cl ,bootMatrix,2,arima,
                                order = c(bootstraps$orders[1],0,bootstraps$orders[2]),
                                method = "ML", 
                                SSinit = c("Rossignol2011"),
                                include.mean = F,
                                seasonal = list(order = c(0, 0, 0), period = NULL),
                                optim.method = "Nelder-Mead",
                                optim.control = list(trace = 1, maxit = 1000))#, fnscale = .1))
        }
        else if(parallel == F){

                res <- apply(bootMatrix,2,arima,
                             order = c(bootstraps$orders[1],0,bootstraps$orders[2]),
                             method = "ML", 
                             include.mean = F,
                             seasonal = list(order = c(0, 0, 0), period = NULL),
                             optim.method = "Nelder-Mead",
                             optim.control = list(maxit = 1000))

        }

        print(noquote("Done."))

        # If the coefs are all 0 the coef list has to be filled manually,
        # because the arima returns coef(fit) = numeric(0).
        if(bootstraps$orders[1] == 0 && bootstraps$orders[2] == 0){
                new_arima_coefs <- matrix(0,ncol = ncol(bootMatrix), nrow = 2)
        }
        else{

                # Get new coefs from resampling
                new_arima_coefs <- sapply(res,coef)
        }
        # Fix the last p observations from the original series.
        last_p_obs <- tail(series, bootstraps$orders[1])

        # Init vector for prediction.
        init <- rep(0,7)
        if(bootstraps$orders[1] != 0){
                init[1:bootstraps$orders[1]] <- rev(last_p_obs)
                init <- rev(init)
        }

        # Convert the coef list back to matrices with only one type (AR or MA).
        coef_list <- apply(new_arima_coefs,2,make_ar_ma_coefs, order = bootstraps$orders)
        coef_list <- sapply(coef_list,unlist)
        ar_coefs <- coef_list[1:7,]
        ma_coefs <- coef_list[8:14,]

        for(i in 1:ncol(bootMatrix)){

                raw[i,] <- calc_Series(ar_coefs[,i], 
                                       ma_coefs[,i], 
                                       arima_residuals = arima_residuals,
                                       is.forecast = T,
                                       series = NA,
                                       n.ahead = n.ahead,
                                       init_vector = init)

        }

        quants <- matrix(ncol = 2, nrow = n.ahead)
        for(i in 1:ncol(raw)){
                q <- quantile(raw[,i], probs = c(.05, .95))
                quants[i,1] <- q[1]
                quants[i,2] <- q[2]
        }
        summary <- apply(raw,2,mean)

        summary <- as.data.frame(cbind(summary,quants))
        names(summary) <- c("Mean","Lo 95", "Hi 95")

        res_tot <- list(rawPredicts = raw, Summary = summary)
        class(res_tot) <- "arimaBootstrapSummary"

        # Close multi-core cluster
        if(parallel == T){stopCluster(cl)}

        res_tot
}

# NOT RUN:
# x <- Boots_Predict(f, parallel = T)

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

################## GARCH MODELS#############################################


# Load the rugarch package for GARCH model specs
library(rugarch)

# Source the import and transformations file
source("./import_pre-transforms.R", echo = F)

out.sample <- 30

################### GARCH SPECIFICATION ###############################
#
# Create a list of GARCH specification objects based on the orders of the
# fitted ARIMA models from ./ARIMA.R
get_GARCH_spec_list <- function(autoArimaList, garch.model = "gjrGARCH") {

        stopifnot(class(autoArimaList) == "autoArimaList")

        armaOrders <- get_orders_from_arima_fit(autoArimaList)[,2:3]

        spec.list <- list()

        for(i in 1:length(autoArimaList)){
                spec.list[[i]] <- ugarchspec(variance.model =
                                             list(model = garch.model,
                                                  garchOrder = c(1,1),
                                                  variance.targeting = T),
                                             mean.model = 
                                                     list(armaOrder = c(armaOrders[i,1],
                                                                        armaOrders[i,2]),
                                                          include.mean = F),
                                             distribution.model = "sstd")
        }
        names(spec.list) <- names_complete
        class(spec.list) <- "ugarchSpecList"
        return(spec.list)
}
garch.specs <- get_GARCH_spec_list(auto.arima.list)

################### GARCH FITTING ####################################

# Fit the specified models to the data 
fit_GARCH_List <- function(ugarchSpecList, data, out.sample = 10){

        stopifnot(class(ugarchSpecList) == "ugarchSpecList")

        fitted.list <- list()
        for(i in 1:length(ugarchSpecList)){
                print(noquote(paste0("Fitting ", names_complete[i], "...")))
                fitted.list[[i]] <- ugarchfit(ugarchSpecList[[1]],
                                              data = data[[i]],
                                              out.sample = out.sample)
        }
        names(fitted.list) <- names_complete
        return(fitted.list)
}

garch.fits <- fit_GARCH_List(garch.specs,data = diff.list, out.sample = out.sample)

################## GARCH FORECASTING ##################################
# #
# # Here I forecast the non-bootstrapped GARCH model in order to obtain
# # the in-sample MAE which is used to calculate the MASE
# forecast_garch_from_fit_list <- function(GARCHfitList, n.ahead = 10){
#   
#   forc <- list()
#   for(i in 1:length(GARCHfitList)){
#     print(noquote(paste0("Forecasting ", names_complete[i])))
#     forc[[i]] <- ugarchforecast(GARCHfitList[[i]],n.ahead = n.ahead)
#   }
#   names(forc) <- names_complete
#   return(forc)
# }
# 
# garch.forecast <- forecast_garch_from_fit_list(garch.fits,n.ahead = out.sample)

################## GARCH BOOTSTRAPPING #################################
#
# This function implements the bootstrapping process from Pascual et. al (2006).
fit_garch_bootstrap_list <- function(garchFITobjectLIST, n.ahead = 10){

        bootstrap <- list()

        for(i in 1:length(garchFITobjectLIST)){

                print(noquote("Initializing cluster..."))
                cl  <-  makeCluster(detectCores() - 1)

                print(noquote(paste0("Fitting ",names_complete[i],"...")))
                bootstrap[[i]] <- ugarchboot(garchFITobjectLIST[[i]], 
                                             method = "Partial",
                                             n.ahead = n.ahead,
                                             #n.bootpred = 0,
                                             n.bootfit = 1000, 
                                             cluster = cl)
                stopCluster(cl)
                print(noquote("--------------------------------------"))
        }
        names(bootstrap) <- names_complete
        print(noquote("Done!"))
        class(bootstrap) <- "garchBOOTlist"
        return(bootstrap)
}

garch.boots <- fit_garch_bootstrap_list(garch.fits, n.ahead = 30)

################ ERROR METRICS ################################################
#
# Get MASE scaling factor from in-sample one-step naive forecast MAE of series.
# (Hyndman and Koehler, 2005)
MASE_scaling_factor <- function(training.data) {
        naive.forecast <- vector()
        for(i in 1:(length(training.data)- 1)){
                naive.forecast[i] <- training.data[i+1]-training.data[i]
        }
        scaling.factor <- sum(abs(naive.forecast)) / (length(training.data) - 1)
        return(scaling.factor)
}

# This function retrieves error metrics from the bootstrapped GARCH forecasts
# compared to the testing sets.
error_metrics_garch_bootstrap <- function(garchBOOTlist, trainingset, testsets) {

        # Check if supplied class is correct
        stopifnot(class(garchBOOTlist) == "garchBOOTlist")

        # Initialize result matrix
        resdf <- data.frame(matrix(nrow = length(garchBOOTlist), ncol = 3))

        for(entry in 1:length(garchBOOTlist)){

                # Get mean of bootstrapped forecasts
                mean.boot.forecast <- apply(garchBOOTlist[[entry]]@fseries, 2, mean)
                training.data <- trainingset[[entry]]
                testing.data <- testsets[[entry]]
                # Get MASE scaling factor
                MASE.scaling.factor <- MASE_scaling_factor(training.data)

                # RMSE
                resdf[entry,1] <- accuracy(mean.boot.forecast,testing.data)[2]
                #MASE
                sum_abs_err <- sum(abs(testing.data - mean.boot.forecast))
                resdf[entry,2] <- (sum_abs_err / length(testing.data)) / MASE.scaling.factor
                # MAE
                resdf[entry,3] <- accuracy(mean.boot.forecast,testing.data)[3]
        }
        names(resdf) <- c("RMSE","MASE","MAE")
        return(resdf)
}

garch.errors <- error_metrics_garch_bootstrap(garch.boots,train.test.diff$train,
                                                       train.test.diff$test)










# Apply test Garch Model
padded.forecast.vector <- rep(NA,length(diff.list[[1]])-out.sample)
padded.forecast.vector <- append(padded.forecast.vector,garch.bootp@fseries[1,])

plot(diff.list$CAC40, type = "b", xlim = c(1150,1250), lwd = 2)
lines(padded.forecast.vector, type = "b", col = "green4", lwd = 2)
abline(v = 1220)
