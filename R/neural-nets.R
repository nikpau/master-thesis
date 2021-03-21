setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Source the import and transformations file
source("./import_pre-transforms.R", echo = F)

# Source helper functions
source("./aux/helper.R")

out.sample <- 7

# Split the ts.list into two parts because of memory usage (16Gb on my machine 
# is not enough to train in one batch)
ts.list_1_7 <- list()
ts.list_8_14 <- list()
for (i in  seq_len(length(ts.list) / 2)) {
        ts.list_1_7[[i]] <- ts.list[[i]]
        ts.list_8_14[[i]] <- ts.list[[i + 7]]
}

nn.gridsearch_1_7 <- lag_list_nn(ts.list_1_7, c(2,4,7,14), out.sample)
nn.gridsearch_8_14 <- lag_list_nn(ts.list_8_14, c(2,4,7,14), out.sample)

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

        # Iterate over indices
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
        #names(ind) <- names_complete
        return(ind)
}
flags.mlp <-  list(dropout = c(.1),
                   dense.units = c(16, 32, 64, 128),
                   val.split = c(.1, .2),
                   neg.slope = c(0, 0.02))

flags.lstm <- list(dropout = c(0.1),
                   lstm.units1 = c(16,32,64),
                   lstm.units2 = c(16,32,64),
                   val.split = c(0.2),
                   neg.slope = c(0,0.02))
mlp_grid_search <- super_grid_search(nn.gridsearch,
                                     flags = flags.mlp,
                                     "./schedules/mlp_train.R")

lstm_grid_search_1_7 <- super_grid_search(nn.gridsearch_1_7,
                                      flags = flags.lstm,
                                      "./schedules/lstm_train.R")

lstm_grid_search_8_14 <- super_grid_search(nn.gridsearch_8_14,
                                      flags = flags.lstm,
                                      "./schedules/lstm_train.R")

lstm1 <- readRDS("./_objects/lstm_1_7_fit.rds")
lstm2 <- readRDS("./_objects/lstm_8_14_fit.rds")

lstm_gridsearch <- c(lstm1, lstm2)
names(lstm_gridsearch) <- names_complete

# Extract the first and second best run
# (second one for comparison)
get_best_runs <- function(runsList) {

        reduced <- lapply(runsList, do.call, what = rbind.data.frame)
        best <- list()
        for (i in seq_len(length(reduced))) {

                ordered <- reduced[[i]][order(reduced[[i]]$metric_val_loss), ]
                best[[i]] <- ordered[c(1,2), ]

        }
        names(best) <- names_complete
        return(best)

}
ppp <- get_best_runs(lstm_gridsearch)
