# import necessary library
library(rsample)
#
# List functionality for creating a rolling origin with arbitrary train, 
# test and skip lengths (creates "splits").
roll_or <- function(list, train, test, skip){
        roll_ls <- list()
        for (i in 1:length(list)) {
                roll_ls[[i]] <- rsample::rolling_origin(as.data.frame(list[[i]]),
                                                        initial = train,
                                                        assess = test,
                                                        skip = skip,
                                                        cumulative = F)
        }
        class(roll_ls) <- "diff_split"
        roll_ls
}
