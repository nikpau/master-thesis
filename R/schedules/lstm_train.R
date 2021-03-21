library(keras)

# Hyperparameter flags ----------------------------

FLAGS <- flags(
  flag_numeric("dropout", 0.2),
  flag_numeric("lstm.units1", 16),
  flag_numeric("lstm.units2", 16),
  flag_numeric("val.split", 0.1),
  flag_numeric("neg.slope", 0)
)

# Data preparation -------------------------------

# Data preparation takes place in in neural-nets.R.
# In here I always use one x.train and y.train for training
# which are fed sequentially into this script via the super_grid_search
# function.

# Reshaping inputs to arrays for LSTM model
x.train.arr <- array(x.train, dim = c(nrow(x.train), ncol(x.train), 1))
y.train.arr <- array(y.train, dim = c(length(y.train), 1))

# Define model ----------------------------------

lstm_model <- keras_model_sequential()
lstm_model %>%
  layer_lstm(units = FLAGS$lstm.units1, 
             input_shape = c(dim(x.train)[2], 1),
             return_sequences = T,
             stateful = F) %>%
  layer_dropout(rate = FLAGS$dropout) %>%
  layer_lstm(units = FLAGS$lstm.units2, 
             return_sequences = F,
             stateful = F) %>%
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
  validation_split = FLAGS$val.split
  #batch_size = batch_size
)
