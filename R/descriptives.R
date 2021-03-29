########################### RAW PLOTTING #############################
# Plot function for the raw series and the differenced series 
# for visual inspection.
plot_Raw_Series <- function(differencing = F ){

        ts <- date_ts
        ts_diff <- diff_list

        for (entry in seq_len(length(ts_diff))) {
                print(noquote(paste0("Plotting series ", names_complete[entry])))
                if(differencing == F){
                        pdf(paste0("./img/raw_series/undifferenced/",
                                   names_complete[entry],
                                   "_ts_raw.pdf"), height = 5)

                        plot(ts[[1]],ts[[entry+1]], type = "l", lwd = 1.5, 
                             main = paste0("Original Series for ", 
                                           names_complete[entry]), 
                             xlab = "Date", ylab = "Value")
                        grid()
                        dev.off()
                }
                else{
                        pdf(paste0("./img/raw_series/differenced/",
                                   names_complete[entry],
                                   "_ts_diff.pdf"), height = 5)

                        plot(ts_diff[[entry]], type = "l", lwd = 1, 
                             main = paste0("(0,1,0) Series for ", 
                                           names_complete[entry]), 
                             xlab = "Index", ylab = "Value")
                        grid()
                        dev.off()
                }

        }

}

# Create Boxplots for the time series
make_Boxplot <- function(list) {

        for (i in seq_len(length(list))) {

                pdf(file = paste0("./img/raw_series/boxplots/",
                                  names_complete[i],
                                  "_boxPlot.pdf"), width = 4)
                boxplot(list[[i]], main = paste0("Boxplot for ", names_complete[i]),
                        xlab = NULL, range = 2)
                abline(h = 0)
                dev.off()

        }

}

# Create autocorrelation function plots for each time series and save them into the /img folder
save_ACF <- function(list,n.lag = 30, type = "acf") {

        for (i in seq_len(length(list))) {

                if(type == "acf") {
                        print(noquote(paste0("Working on ACF for", 
                                             names(diff_list)[i])))
                        index_name <- names(list)[i]

                        pdf(file = paste0("./img/ACFs/ACF"," ",
                                          index_name,
                                          ".pdf"), height = 5)
                        plot(forecast::Acf(list[[i]],lag.max = n.lag, plot = F), 
                             main = paste("ACF for",index_name,
                                          "(n.lag = ",n.lag,")"), 
                             ylim = c(-.1,.1), lwd = 2)

                        dev.off()
                }
                else if (type  == "pacf") {
                        print(noquote(paste0("Working on PACF for", 
                                             names(diff_list)[i])))
                        index_name <- names(list)[i]

                        pdf(file = paste0("./img/PACFs/PACF"," ",
                                          index_name,
                                          ".pdf"), height = 5)
                        plot(forecast::Pacf(list[[i]],lag.max = n.lag, plot = F), 
                             main = paste("PACF for",index_name,
                                          "(n.lag = ",n.lag,")"), 
                             ylim = c(-.1,.1), lwd = 2)

                        dev.off()
                }
                else {
                        stop("Type can only be one of either 'acf' or 'pacf' ")
                }
        }

        print(noquote(paste0("Done!")))
}

########################### Stationarity #############################

# The following test all require the input-lists to only have one atomic structure per list entry
# inside. Since this is a tailored function just for the purpose of my master thesis I will not
# include any checks for valid input arguments.
#
# Perform ADF test and print results in table. Multicore support (4 cores as default)
list_ADF <- function(list, n.cores = 4) {

        adf_results <- mclapply(list, 
                                aTSA::adf.test,
                                output = F,
                                mc.preschedule = T,
                                mc.cores = n.cores)
        res  <- matrix(ncol = 3)
        index_name <- vector()

        for(i in 1:length(adf_results)){

                type2 <- rep("type1", length(adf_results[[i]][[1]][, 1]))
                type2 <- rep("type2", length(adf_results[[i]][[2]][, 1]))
                type3 <- rep("type3", length(adf_results[[i]][[3]][, 1]))

                temp1  <- cbind(type1, adf_results[[i]][[1]][, c(1, 3)])
                temp2  <- cbind(type2, adf_results[[i]][[2]][, c(1, 3)])
                temp3  <- cbind(type3, adf_results[[i]][[3]][, c(1, 3)])

                res <- rbind(res,temp1,temp2,temp3)
                index_name <- append(index_name, rep(names(list)[i],
                                                    (length(type1) +
                                                     length(type2) +
                                                     length(type3))))


        }

        res <- cbind(index_name, res[-1, ])
        names(res) <- c("Index", "test_type", "lag", "p.val")
        as.data.frame(res)

}

# KPSS-Test
list_KPSS <- function(list) {

        KPSS_results <- lapply(list, aTSA::kpss.test, output = F)
        res  <- matrix(ncol = 2)
        index_name <- vector()
        type_name <- vector()

        for(i in 1:length(KPSS_results)){

                index_name <- append(index_name, rep(names(list)[i], length(KPSS_results[[i]][,1])))
                type_name <- append(type_name, row.names(KPSS_results[[i]]))
                res  <- rbind(res, KPSS_results[[i]][,c(1,3)])
        }

        res <- cbind(index_name,type_name,res[-1, ])
        row.names(res) <- NULL
        as.data.frame(res)

}

###################### Tests for Normality ###################################

# Jarque-Bera test
list_JB <- function(list) {

        result <- as.data.frame(matrix(ncol = 2))
        index_name <- vector()

        for(i in seq_len(length(list))) {

                test <- moments::jarque.test(list[[i]])
                temp_res <- c(test[[1]], test[[2]])
                index_name <- append(index_name, names(list)[i])

                result <- result %>%
                        rbind(temp_res)


        }

        result <- cbind(index_name,result[-1, ])
        names(result) <- c("Index", "Test Statistc", "p.val")
        row.names(result) <- NULL
        result
}

################### Distribution Tests / Distr. params estimation ############

# Define likelihood functions for three different univariate distributions.
# (Skew-normal, Skew-t, Skew-Chauchy)
# I use the 'sn' package for convenience as the distr. functions are implemented in C 
# for computational effciency.
#
# Likelihood functions | NOTE: the functions will be implemented as negative log-lik funcs.
#
# Skew-normal
sn_llik <- function(par,data){

        -sum(log(sn::dsn(data,xi = par[1],omega = par[2],alpha = par[3])))

}
# Skew-t
st_llik <- function(par,data){

        -sum(log(sn::dst(data,xi = par[1],omega = par[2],alpha = par[3],nu = par[4])))

}
# Skew-Cauchy
sc_llik <- function(par,data){

        -sum(log(sn::dst(data,xi = par[1],omega = par[2],alpha = par[3], log = F)))

}

# Fit data to distribution using the Nelder-Mead algorithm
# Returns a data frame with all indices and their respective estimated parameters. 

# Starting values (arbitrary)
p0_1 <- c(.02,.03,1)
p0_2 <- c(.02,.03,1,4)

# Function
est_Par <- function(list){

        res <- as.data.frame(matrix(ncol = 4 + 3 + 3, nrow = length(list)))
        index_name <- vector()

        # Information Criteria result matrices
        AIC_res <- as.data.frame(matrix(ncol = 3))
        BIC_res <- as.data.frame(matrix(ncol = 3))
        HQC_res <- as.data.frame(matrix(ncol = 3))
        ic_names <- c("Index","Skew-Normal","Skew-t","Skew-Cauchy")

        # Information criteria functions
        AIC <- function(optim) {2 * optim$value + 2 * length(optim$par)}
        BIC <- function(optim,ts) { 2 * optim$value + length(optim$par) * log(length(ts))}
        HQC <- function(optim,ts) {2 * optim$value + 2 * length(optim$par) * log(log(length(ts)))}



        for(i in  seq_len(length(list))) {


                index_name <- append(index_name,names(list)[i])

                sno <- optim(p0_1,
                             sn_llik, 
                             data = list[[i]], 
                             method = "N", 
                             control = list(maxit = 1000))

                sto <- optim(p0_2,
                             st_llik, 
                             data = list[[i]],
                             method = "N", 
                             control = list(maxit = 1000))

                sco <- optim(p0_1,
                             sc_llik, 
                             data = list[[i]], 
                             method = "N", 
                             control = list(maxit = 1000))

                res[i,] <- c(sno$par,
                             sto$par,
                             sco$par)

                # Information Criteria
                # AIC
                AIC_res <- rbind(AIC_res, c(AIC(sno), AIC(sto), AIC(sco)))
                BIC_res <- rbind(BIC_res, c(BIC(sno,list[[i]]), BIC(sto,list[[i]]), BIC(sco, list[[i]])))
                HQC_res <- rbind(HQC_res, c(HQC(sno,list[[i]]), HQC(sto, list[[i]]), HQC(sco,list[[i]])))



        }

        # Bind names vector to the Indices
        AIC_res <- cbind(index_name,AIC_res[-1, ])
        BIC_res <- cbind(index_name,BIC_res[-1, ])
        HQC_res <- cbind(index_name,HQC_res[-1, ])

        # Name columns
        names(AIC_res) <- ic_names
        names(BIC_res) <- ic_names
        names(HQC_res) <- ic_names

        res <- cbind(index_name,res)
        names(res) <- c("Index", "Location (Norm)","Scale (Norm)","Skew (Norm)",
                        "Location (t)","Scale (t)","Skew (t)","df(t)",
                        "Location (Cauchy)","Scale (Cauchy)","Skew (Cauchy)")
        res_list <- list(Estimates = res, 
                         InformationCriteria = list(AIC = AIC_res,
                                                    BIC = BIC_res,
                                                    HQC = HQC_res))
        res_list

}
