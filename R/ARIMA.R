# Function for finding the best ARIMA model based on some information
# criterion.
#
# Basically just a wrapper for the forecast::auto.arima function with added list capability.
#
# (Compilation Time: 1 minute for ~14*1200 obs. on 16 cores)
list_Auto_Arima <- function(list,out.sample = 30,crit = c("aicc","aic","bic"), parallel = F, n.cores = 2) {

        # Result list
        res <- list()

        for (entry in 1:length(list)) {

                cat(paste0("Calculating ", names_complete[entry], "\n"))

                # Define training sets
                train <- head(list[[entry]], -out.sample)

                res[[entry]] <- forecast::auto.arima(train, d = 1, D = 0,
                                                     max.p = 7,
                                                     max.q = 7,
                                                     max.order = 14,
                                                     max.d = 1, max.D = 0,
                                                     stationary = F,
                                                     seasonal = F,
                                                     stepwise = F,
                                                     approximation = F,
                                                     allowdrift = F,
                                                     allowmean = F,
                                                     parallel = parallel,
                                                     num.cores = n.cores)


        }
        names(res) <- names_complete
        class(res) <- "autoArimaList"
        res
}

info_crits_arima <- function(arimaFITorLIST, name = "Series") {

        resdf <- data.frame(matrix(ncol = 4, nrow = 1))

        if(class(arimaFITorLIST) == "autoArimaList"){
                resdf <- data.frame(matrix(ncol = 4, nrow = length(arimaFITorLIST)))
                for (entry in seq_len(length(arimaFITorLIST))) {
                        resdf[entry,2] <- arimaFITorLIST[[entry]]$aic
                        resdf[entry,3] <- arimaFITorLIST[[entry]]$aicc
                        resdf[entry,4] <- arimaFITorLIST[[entry]]$bic
                }
                resdf[,1] <- names(arimaFITorLIST)
        }
        else{
                resdf[1,2] <- arimaFITorLIST$aic
                resdf[1,3] <- arimaFITorLIST$aicc
                resdf[1,4] <- arimaFITorLIST$bic
                resdf[,1] <-  name
        }
        names(resdf) <- c("Index", "AIC","AICc","BIC")
        return(resdf)
}

# Write the ARIMA orders into a data frame.
get_orders_from_arima_fit <- function(arimaFITorList) {

        if(class(arimaFITorList) == "autoArimaList"){
                resdf <- data.frame(matrix(nrow = length(arimaFITorList), 
                                           ncol = 3))
                for(entry in 1:length(arimaFITorList)){
                        resdf[entry,c(2,3)] <- arimaFITorList[[entry]][["arma"]][1:2]
                }
                resdf[,1] <- names_complete
                names(resdf) <- c("Index","AR(p)","MA(q)")
        }
        else{
                resdf <- arimaFITorList[["arma"]][c(1,2)]
                names(resdf) <- c("AR(p)","MA(q)")
        }
        return(resdf)
}

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

########## ERROR METRICS #############################

# Get Error metrics for forecast values.
error_metrics_arima <- function(forecasts, testsets, trainsets) {

        resdf <- data.frame(matrix(ncol = 3, nrow = length(forecasts)))

        for (entry in seq_len(length(forecasts))) {
                
                test <- testsets[[entry]]
                train <- trainsets[[entry]]
                forc <- forecasts[[entry]]$mean
                MASE_scaling_factor <- MASE_scaling_factor(train)
                
                #MASE
                sum_abs_err <- sum(abs(test - forc))
                mase <- (sum_abs_err / length(test)) / MASE_scaling_factor
                resdf[entry, 2] <- mase
                
                #RMSSE
                sum_err <- sum(test - forc)
                rmsse <- sqrt(((sum_err / length(test)) / MASE_scaling_factor)^2)
                resdf[entry, 1] <- rmsse
                
                #MdASE
                med_abs_err <- median(abs(test - forc))
                mdase <- med_abs_err / MASE_scaling_factor
                resdf[entry, 3] <- mdase
        }
        names(resdf) <- c("MASE", "RMSSE", "MdASE")
        return(resdf)
}

save_resid_qq <- function(autoARIMAList) {

        if (class(autoARIMAList) != "autoArimaList")
                stop("Input class must be an auto Arima List. \n
                     Please use the list_Auto_Arima function")

                     for (i in seq_len(length(autoARIMAList))) {

                             pdf(file = paste0("./img/residual_qq/",
                                               names_complete[i],
                                               "_resid_qq_norm.pdf"))
                             qqnorm(autoARIMAList[[i]]$residuals,
                                    main = paste0("Residual QQ-Normal Plot for ",
                                                  names_complete[i]),
                                    lwd = 2)
                             qqline(autoARIMAList[[i]]$residuals)
                             dev.off()

                     }
}
