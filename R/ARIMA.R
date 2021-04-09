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
                                                     max.p = 12,
                                                     max.q = 12,
                                                     max.order = 24,
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
                resdf[entry, 1] <- calculate_mase(test, train, forc, MASE_scaling_factor)
                
                #RMSSE
                resdf[entry, 2] <- calculate_rmsse(test, train, forc, MASE_scaling_factor)
                
                #MdASE
                resdf[entry, 3] <- calculate_mdase(test, train, forc, MASE_scaling_factor)
        }
        names(resdf) <- c("MASE", "RMSSE", "MdASE")
        return(resdf)
}

# Save QQ-Normal plots for the residuals of the arima fits
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
