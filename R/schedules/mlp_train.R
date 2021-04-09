library(keras)

# Hyperparameter flags ----------------------------

FLAGS <- flags(
               flag_numeric("dropout", 0.2),
               flag_numeric("dense.units", 16),
               flag_numeric("val.split", 0.1),
               flag_numeric("neg.slope", 0)
)

# Data preparation -------------------------------

# Data preparation takes place in in the nn_gridsearch.R and
# neural-nets.R.
# In here I always just one x.train and y.train for training

# Define model ----------------------------------

model <- keras_model_sequential()
model %>%
        layer_dense(units = FLAGS$dense.units,
                    input_shape = dim(x.train)[2]) %>%
        layer_activation_relu(negative_slope = FLAGS$neg.slope) %>%
        layer_dense(units = FLAGS$dense.units) %>%
        layer_activation_relu(negative_slope = FLAGS$neg.slope) %>%
        layer_dense(units = 1)

model %>%
        compile(
                loss = "mse",
                optimizer = optimizer_adam(lr = 0.001)
        )

# Training & Eval ------------------------------

history <- model %>% fit(
                         x.train, y.train,
                         epochs = 100,
                         verbose = 2,
                         callback = callback_early_stopping(monitor = "val_loss",
                                                            patience = 10),
                         validation_split = FLAGS$val.split
)