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

# Get orders from the rugarch fit
get_orders_rugarch_arima <- function(rugarchFitlList) {

        coefs <- lapply(rugarchFitlList, function(x) return(coef(x$fit)))

        resdf <- data.frame(matrix(nrow = length(rugarchFitlList), 
                                   ncol = 3))

        for (i in seq_len(length(rugarchFitlList))) {

                names <- names(coefs[[i]])
                names <- head(names, -3)
                ar_order <- sum(str_count(names, pattern = "ar"))
                ma_order <- sum(str_count(names, pattern = "ma"))

                resdf[i,c(2,3)] <- c(ar_order, ma_order)
        }
        resdf[,1] <- names_complete
        names(resdf) <- c("Index","AR(p)","MA(q)")
        return(resdf)

}

########## ERROR METRICS #############################

# Get Error metrics for forecast values.
# The forecasts coming from the arima fit with normal distribution 
# will be retransformed right here. 
error_metrics_arima <- function(forecasts, testsets, trainsets) {

        resdf <- data.frame(matrix(ncol = 3, nrow = length(forecasts)))

        for (entry in seq_len(length(forecasts))) {



                if (any(class(forecasts[[1]]) == "forecast")) {
                        forc <- exp(forecasts[[entry]]$mean)
                        test <- exp(testsets[[entry]])
                        train <- exp(trainsets[[entry]])
                }
                else {
                        forc <- forecasts[[entry]]
                        test <- testsets[[entry]]
                        train <- trainsets[[entry]]
                }
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

        op <- par(no.readonly = T)


        if (class(autoARIMAList) != "autoArimaList")
                stop("Input class must be an auto Arima List. \n
                     Please use the list_Auto_Arima function")

                     for (i in seq_len(length(autoARIMAList))) {

                             pdf(file = paste0("./img/residual_qq/",
                                               names_complete[i],
                                               "_resid_qq_norm.pdf"))
                             par(mar = c(4,6,4,4)+.1)
                             qqnorm(autoARIMAList[[i]]$residuals,
                                    main = paste0("Residual QQ-Normal Plot for ",
                                                  names_complete[i]),
                                    lwd = 2, cex = 2, cex.axis = 2, cex.lab = 2,
                                    cex.main = 2, pch = 18)
                             qqline(autoARIMAList[[i]]$residuals)
                             grid()
                             dev.off()

                     }
                     par(op)
}

# Function for extracting the forecasts from a rugarch forecast class, glue
# them to the original series and writing them into a list.
sstd_retransform <- function(rugarchForecast, trainingData, const_int) {

        result <- list()
        for (i in seq_len(length(rugarchForecast))) {
                stitched <- c(trainingData[[i]], rugarchForecast[[i]]@forecast$seriesFor)
                result[[i]] <- exp(diffinv(stitched, xi = const_int[i]))
        }
        return(result)
}

# Plot the p-values of JB tests for different lag-values
plot_JB_lags <- function(residualList, up.to = 20, names, filename) {

        p_vals <- matrix(ncol = up.to, nrow = length(residualList))

        for (k in seq_len(length(residualList))){
                for (i in seq_len(up.to)) {
                        p_vals[k,i] <- Box.test(residualList[[k]], lag = i, type = "L")$p.value

                }
        }

        pdf(file = paste0("./img/residual_JB/",filename,".pdf"), 
            height = 5, family = "Courier")

        plot(p_vals[1,], type = "b", pch = 12, xaxt = "n",
             ylim = c(0,1), lwd = 1.2, main = paste("p-values for JB-test up to Lag",up.to),
             ylab = "p.value",xlab = "Lag")
        axis(1, seq(1,20,2))
        grid(lty = 1)

        for (j in seq_len(nrow(p_vals[-1,]))) {

                lines(p_vals[j+1,], type = "b", pch = j+1,
                      ylim = c(0,1), lwd = 1.2)

        }

        abline(h = .05, lty = 4, lwd = 2, col = "#B58900")


        legend("topright", legend = names,
               pch = c(12, seq(2,nrow(p_vals[-1,]) + 1)),
               bg = "white",
               lwd = 1)

        dev.off()

}
























