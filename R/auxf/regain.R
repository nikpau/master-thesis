# Regain numeric data from split
regain <- function(single_split){
  train <- rsample::training(single_split) %>% 
    unlist()
  test <- rsample::testing(single_split) %>% 
    unlist()
  list(train = train,test = test)
}
