# Retrieve information criteria from an ets fit.
error_metrics_ets <- function(etsLIST, training_set, testing_set) {

        resdf <- data.frame(matrix(ncol = 3, nrow = length(etsLIST)))

        for (i in seq_len(length(etsLIST))) {

                test <- testing_set[[i]]
                train <- training_set[[i]]
                forc <- etsLIST[[i]]$mean
                MASE_scaling_factor <- MASE_scaling_factor(train)

                #MASE
                resdf[i, 2] <- calculate_mase(test, train, forc, MASE_scaling_factor)
                
                #RMSSE
                resdf[i, 1] <- calculate_rmsse(test, train, forc, MASE_scaling_factor)
                
                #MdASE
                resdf[i, 3] <- calculate_mdase(test, train, forc, MASE_scaling_factor)

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
