# Retrieve information criteria from an ets fit.
error_metrics_ets <- function(etsLIST, training_set, testing_set) {

        resdf <- data.frame(matrix(ncol = 3, nrow = length(etsLIST)))

        for (i in seq_len(length(etsLIST))) {

                test <- testing_set[[i]]
                train <- training_set[[i]]
                forc <- etsLIST[[i]]$mean
                MASE_scaling_factor <- MASE_scaling_factor(train)

                #MASE
                sum_abs_err <- sum(abs(test - forc))
                mase <- (sum_abs_err / length(test)) / MASE_scaling_factor
                resdf[i, 2] <- mase

                #RMSSE
                sum_err <- sum(test - forc)
                rmsse <- sqrt(((sum_err / length(test)) / MASE_scaling_factor)^2)
                resdf[i, 1] <- rmsse

                #MdASE
                med_abs_err <- median(abs(test - forc))
                mdase <- med_abs_err / MASE_scaling_factor
                resdf[i, 3] <- mdase

        }
        names(resdf) <- c("MASE", "RMSSE", "MdASE")
        return(resdf)

}

# Function for receiving information criteria from an
# ets model.
info_crits_ets <- function(etsList) {
        resdf <- data.frame(matrix(ncol = 3, nrow = length(etsList)))

        for (i in seq_len(length(etsList))) {
                resdf[i,1] <- etsList[[i]][["model"]][["aic"]]
                resdf[i,2] <- etsList[[i]][["model"]][["aicc"]]
                resdf[i,3] <- etsList[[i]][["model"]][["bic"]]
        }
        names(resdf) <- c("AIC", "AICc", "BIC")
        return(resdf)
}
