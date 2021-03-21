
# Set working directoty and enable relative paths.
setwd(getwd())

# If you use RStudio please uncomment the following line
# to enable relative paths.
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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
