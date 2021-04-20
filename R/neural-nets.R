# In this function I gridsearch through different lags
# for the neural network in order to later
# optimize hyperparameters, based on the best training runs.

# Test different values for the number of 
# units in the hidden layer and the dropout rate for every 
# lag returned from the lag_list_nn function and every index.
super_grid_search <- function(nnLagGrid, flags, train.schedule) {
        
        if (class(nnLagGrid) != "nnLagGrid")
                stop("Only nnLagGrids allowed as inputs.")
        
        ind <- list()
        naming <- vector()
        
        # Iterate over time series
        for (i in seq_len(length(nnLagGrid[[1]]$train$x.train))) {
                lags <- list()
                
                # Iterate over lags
                for (j in seq_len(length(nnLagGrid))) {
                        
                        naming[j] <- paste0("Lag ",
                                            dim(nnLagGrid[[j]][["train"]][["x.train"]][[1]])[2])
                        
                        # Extract training and testing sets from the nnLagGrid
                        x.train <- nnLagGrid[[j]]$train$x.train[[i]]
                        y.train <- nnLagGrid[[j]]$train$y.train[[i]]
                        x.test <- nnLagGrid[[j]]$test$x.test[[i]]
                        y.test <- nnLagGrid[[j]]$test$y.test[[i]]
                        
                        lags[[j]] <- tuning_run(train.schedule,
                                                flags = flags,
                                                confirm = F,
                                                envir = environment(),
                                                echo = F)
                        
                }
                names(lags) <- naming
                ind[[i]] <- lags
        }
        return(ind)
}

# Add a column to the tf_runs df with the lag at which the run was done
add_lag <- function(runsList, lag.vector) {
        
        if (length(runsList[[1]]) != length(lag.vector))
                stop("wrong length for lag vector")
        
        for (i in seq_len(length(runsList))) {
                
                for (j in seq_len(length(runsList[[i]]))) {
                        
                        runsList[[i]][[j]][["Lag"]] <- rep(lag.vector[j],
                                                           nrow(runsList[[i]][[j]]))
                }
        }
        
        return(runsList)
}

# Extract the first and second best run
# (second one for comparison)
get_best_runs <- function(runsList) {
        
        reduced <- lapply(runsList, do.call, what = rbind.data.frame)
        best <- list()
        for (i in seq_len(length(reduced))) {
                
                ordered <- reduced[[i]][order(reduced[[i]]$metric_val_loss), ]
                best[[i]] <- ordered[1, ]
                
        }
        names(best) <- names_complete
        return(best)
        
}

# Extract flags from best run for training
extract_flags <- function(runsList) {
        
        flags <- list()
        for (i in seq_len(length(runsList))) {
                flags[[i]] <- runsList[[i]] %>%
                        select(tidyselect::contains(c("flag", "Lag")))
        }
        return(flags)
}

train_and_predict <- function(series, flags, out.sample, network.type) {
        
        pred_out <- list()
        
        for (i in seq_len(length(series))) {
                
                # Preprocess the data
                processed <- pre_process(series[[i]],flags[[i]]$Lag, out.sample)
                
                if (network.type == "lstm") {
                        
                        #Train the network 
                        trained <- do_train_lstm(processed, flags[[i]])
                        
                        # Reshape test array
                        x.test.arr <- array(processed$test$x.test[[1]], 
                                            dim = c(nrow(processed$test$x.test[[1]]), 
                                                    ncol(processed$test$x.test[[1]]), 1))
                        
                        # Predict 
                        predicted <- trained$model %>% 
                                predict(x.test.arr) %>% 
                                .[,1]
                        
                }
                else if (network.type == "mlp") {
                        trained <- do_train_mlp(processed, flags[[i]])
                        
                        # Predict 
                        predicted <- trained$model %>% 
                                predict(processed$test$x.test[[1]]) %>% 
                                .[,1]
                        
                }
                else {
                        stop("wrong network type")
                }
                
                
                scale_history <- processed$history$scale.history
                center_history <- processed$history$center.history
                
                pred_out[[i]] <- list(
                        forc = exp(predicted * scale_history + center_history),
                        train = exp(processed$train$y.train[[1]] * scale_history + center_history),
                        test = exp(processed$test$y.test[[1]] * scale_history + center_history)
                )
                
        }
        return(pred_out)
        
}

# Literal training procedure for the LSTM network
do_train_lstm <- function(data, flags) {
        
        # Reshaping inputs to arrays for LSTM model
        x.train.arr <- array(data$train$x.train[[1]], 
                             dim = c(nrow(data$train$x.train[[1]]), 
                                     ncol(data$train$x.train[[1]]), 1))
        y.train.arr <- array(data$train$y.train[[1]], 
                             dim = c(length(data$train$y.train[[1]]), 1))
        
        # Define model ----------------------------------
        
        lstm_model <- keras_model_sequential()
        lstm_model %>% 
                layer_lstm(units = flags$flag_lstm.units1,
                           input_shape = c(dim(data$train$x.train[[1]])[2], 1),
                           return_sequences = T,
                           stateful = F) %>% 
                layer_activation_relu(negative_slope = flags$flag_neg.slope) %>%
                layer_dropout(rate = flags$flag_dropout) %>%
                layer_lstm(units = flags$flag_lstm.units2, 
                           return_sequences = F,
                           stateful = F) %>%
                layer_activation_relu(negative_slope = flags$flag_neg.slope) %>%
                layer_dense(units = 1)
        
        lstm_model %>%
                compile(
                        loss = "mse",
                        optimizer = optimizer_adam(lr = 0.001)
                )
        # Training & Eval ------------------------------
        
        history <- lstm_model %>% fit(
                x.train.arr, y.train.arr,
                epochs = 100,
                verbose = 0,
                callback = callback_early_stopping(monitor = "val_loss",
                                                   patience = 5),
                validation_split = flags$flag_val.split
                #batch_size = batch_size
        )
        return(list(model = lstm_model, history = history))
}

do_train_mlp <- function(data, flags) {
        
        mlp_model <- keras_model_sequential()
        mlp_model %>%
                layer_dense(units = flags$flag_dense.units,
                            input_shape = dim(data$train$x.train[[1]])[2]) %>%
                layer_activation_relu(negative_slope = flags$flag_neg.slope) %>%
                layer_dense(units = flags$flag_dense.units) %>%
                layer_activation_relu(negative_slope = flags$flag_neg.slope) %>%
                layer_dense(units = 1)
        
        mlp_model %>%
                compile(
                        loss = "mse",
                        optimizer = optimizer_adam(lr = 0.001)
                )
        
        # Training & Eval ------------------------------
        
        history <- mlp_model %>% fit(
                data$train$x.train[[1]], data$train$y.train[[1]],
                epochs = 100,
                verbose = 0,
                callback = callback_early_stopping(monitor = "val_loss",
                                                   patience = 10),
                validation_split = flags$flag_val.split
        )
        return(list(model = mlp_model, history = history))
        
}

# Retrieve summarized error metrics from the neural net training process
error_metrics_nn <- function(nnPredictLIst) {
        
        resdf <- data.frame(matrix(ncol = 3, nrow = length(nnPredictLIst)))
        
        for (i in seq_len(length(nnPredictLIst))) {
                
                test <- nnPredictLIst[[i]]$test
                train <- nnPredictLIst[[i]]$train
                forc <- nnPredictLIst[[i]]$forc
                MASE_scaling_factor <- MASE_scaling_factor(train)
                
                #MASE
                resdf[i, 1] <- calculate_mase(test, train, forc, MASE_scaling_factor)
                
                #RMSSE
                resdf[i, 2] <- calculate_rmsse(test, train, forc, MASE_scaling_factor)
                
                #MdASE
                resdf[i, 3] <- calculate_mdase(test, train, forc, MASE_scaling_factor)
                
        }
        names(resdf) <- c("MASE", "RMSSE", "MdASE")
        return(resdf)
        
}

