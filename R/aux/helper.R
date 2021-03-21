# Get MASE scaling factor from in-sample one-step naive forecast MAE of series.
# (Hyndman and Koehler, 2005)
MASE_scaling_factor <- function(training.data) {
        naive.forecast <- vector()
        for (i in 1:(length(training.data) - 1)) {
                naive.forecast[i] <- training.data[i+1]-training.data[i]
        }
        scaling.factor <- sum(abs(naive.forecast)) / (length(training.data) - 1)
        return(scaling.factor)
}

# Function for preparing a list of training and
# testing sets with different lags to feed into the
# neural network.
lag_list_nn <- function(data, lag.vector, out.sample) {

        lag.list <- list()
        for (i in seq_len(length(lag.vector))) {
                cat(paste0("Preparing Lag ", lag.vector[i], "\n"))
                lag.list[[i]] <- pre_process_mlp(data, lag.vector[[i]], out.sample)

        }
        class(lag.list) <- "nnLagGrid"
        return(lag.list)

}

# Pre-process training and testing data in order to to be passed
# into a neural network. This function has list capabilities
# to process all univariate time series at once.
# holdoud = length of test set
# n.lag = numer of past obervations used for one timestep
pre_process_mlp <- function(data, n.lag, holdout) {

        # Data transformations
        transforms <- list()
        history <- list()

        for (i in seq_len(length(data))) {
                df <- tibble(raw.ts = data[[i]])
                rec <- recipe(~., df) %>%
                        step_log(raw.ts) %>% # log transform
                        step_center(raw.ts) %>% # center data to μ = 0
                        step_scale(raw.ts) %>% # scale to variance = 1
                        prep()
                ts.transform <- bake(rec, df) %>%
                        unlist()
                transforms[[i]] <- ts.transform

                # Center and scale history for later back-transformation
                center.history <- rec[["steps"]][[2]][["means"]][["raw.ts"]]
                scale.history <- rec[["steps"]][[3]][["sds"]][["raw.ts"]]
                history[[i]] <- list(center.history = center.history,
                                     scale.history = scale.history)
        }

        # Create lagged x.input and x.test matrices to train and test
        # the neural net.
        if (length(n.lag) == 1) {
                lag.matrices <- lapply(transforms, split_seq, n.in = n.lag)
                lag.matrices <- lapply(lag.matrices, head, n = -1)

                # Split matrices into train and test sets
                x.train <- lapply(lag.matrices, head, n = -holdout)
                x.test <- lapply(lag.matrices, tail, n = holdout)
        }
        else {
                stopifnot(length(n.lag) == length(transforms))
                lag.matrices <- list()
                x.train <- list()
                x.test <- list()
                for (j in seq_len(length(transforms))) {
                        lag.matrices[[j]] <- split_seq(transforms[[j]], n.lag[[j]]) %>%
                                head(-1)
                        x.train[[j]] <- head(lag.matrices[[j]], -n.lag[[j]])
                        x.test[[j]] <- tail(lag.matrices[[j]], n.lag[[j]])
                }

        }

        # Create y.input list based on the lagged values from above
        if (length(n.lag) == 1) {
                y.train <- lapply(transforms, tail, n = -n.lag) %>%
                        lapply(head, n = -holdout)
        }
        else {
                y.train <- list()
                for (k in seq_len(length(transforms))) {
                        y.train[[k]] <- tail(transforms[[k]], -n.lag[[k]]) %>%
                                head(-holdout)
                }
        }

        # Create y.test labels for error checking
        y.test <- lapply(transforms, tail, n = holdout)

        return(list(train = list(x.train = x.train, y.train = y.train),
                    test = list(x.test = x.test, y.test = y.test),
                    history = history))
}
# Retrieve summarized error metrics from the neural net training process
error_metrics_nn <- function(pred, test, train) {

        res <- matrix(ncol = 3)

        #RMSE
        res[1,1] <- accuracy(pred, test)[2]
        MASE.scaling.factor <- MASE_scaling_factor(train)
        sum_abs_err <- sum(abs(test - pred))
        res[1,2] <- (sum_abs_err / length(test)) / MASE.scaling.factor
        res[1,3] <- accuracy(pred, test)[3]

        return(res)

}

# Workhorse of the neural net data transformation
# Function for splitting sequences into a supervised learning
# problem that can be fed into the neural net.
# n.in = # of time-steps fed into the network as input per output.
split_seq <- function(sequence, n.in) {

        in.vector <- sequence[seq_len(n.in)]
        i <- 2
        while (i <= length(sequence) - (n.in - 1)) {
                in.vector <- append(in.vector, sequence[i:(i + n.in - 1)])
                i <- i + 1
        }
        mat.in <- matrix(in.vector, ncol = n.in, byrow = T)
        return(mat.in)

}
