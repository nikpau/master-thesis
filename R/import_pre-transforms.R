# Load data | Origin: Yahoo Finance and investing.com
dir <- "./data"

indicesYahoo <- list(
                     CAC40 <- read_csv(paste0(dir, "/CAC40.csv")),
                     DAX <- read_csv(paste0(dir, "/DAX.csv")),
                     DJI <- read_csv(paste0(dir, "/DJI.csv")),
                     EUROSTOXX50 <- read_csv(paste0(dir, "/EUROSTOXX50.csv")),
                     KOSPI <- read_csv(paste0(dir, "/KOSPI.csv")),
                     NIFTY50 <- read_csv(paste0(dir, "/NIFTY50.csv")),
                     N225 <- read_csv(paste0(dir, "/N225.csv")),
                     SP500 <- read_csv(paste0(dir, "/SP500.csv")),
                     SPASX200 <- read_csv(paste0(dir, "/SPASX200.csv")),
                     SPTSX <- read_csv(paste0(dir, "/SP_TSX-Composite_Canada.csv")),
                     SSE <- read_csv(paste0(dir, "/SSE_Composite.csv")),
                     TA35 <- read_csv(paste0(dir, "/TA35.csv"))
)

indicesInvesting <- list(
                         BOVESPA <- read_csv(paste0(dir, "/Bovespa.csv")),
                         BIST100 <- read_csv(paste0(dir, "/BIST100.csv"))
)

# Define Names
namesYahoo <- c("CAC40","DAX","DJI","EUROSTOXX50","KOSPI","NIFTY50",
                "N225","SP500","SPASX200","SPTSX","SSE","TA35")
namesInvesting <- c("BOVESPA","BIST100")

names_complete <- c(namesYahoo, namesInvesting)

# Name the Lists
names(indicesInvesting) <- namesInvesting
names(indicesYahoo) <- namesYahoo

# Drop unnecessary columns
for (i in seq_len(length(indicesYahoo))) {
        indicesYahoo[[i]] <- indicesYahoo[[i]][,c(1,6)]
}
for (i in seq_len(length(indicesInvesting))) {
        indicesInvesting[[i]] <- indicesInvesting[[i]][,c(1,2)]
}

# Check whether dates are correctly specified
for (i in seq_len(length(indicesInvesting))) {
        indicesInvesting[[i]][,1] <- unlist(indicesInvesting[[i]][,1]) %>%
                as_date(format = "%b %d, %Y")
}

# Join lists
indices <- append(indicesYahoo, indicesInvesting)

# Rename indicies and define the 'value' column to be numeric
for (i in seq_len(length(indices))) {
        names(indices[[i]]) <- c("date", "value")
}

suppressWarnings(
                 for (i in seq_len(length(indices))) {
                         indices[[i]]["value"] <- unlist(indices[[i]]["value"]) %>%
                                 as.double()
                 })

# Delete residue dataframes
rm(CAC40,DAX,DJI,EUROSTOXX50,KOSPI,NIFTY50,N225,
   SP500,SPASX200,SPTSX,SSE,TA35,BOVESPA,BIST100)

#Remove NA
indices <- mclapply(indices, na.omit)

# Order dataframes in list by Date
for (i in seq_len(length(indices))) {
        indices[[i]] <- indices[[i]][order(indices[[i]]$date), ]
}

# Create dataframe with all ts values and a continuous date to inspect data structure
date_ts <- tibble(date = seq(as.Date("2015-01-01"),
                             as.Date("2020-01-01"),
                             by = "days"))
for (i in seq_len(length(indices))) {
        date_ts <- left_join(date_ts, indices[[i]],
                             by = "date")
}
names(date_ts) <- c("date", namesYahoo, namesInvesting)

# Extract only numeric values without dates into a list for further processing
ts_list <- list()
for (i in seq_len(length(indices))) {
        ts_list <- append(ts_list, indices[[i]]["value"])
}
names(ts_list) <- c(namesYahoo, namesInvesting)

# First Difference and log transform
diff_list <- mclapply(ts_list, log) %>%
        mclapply(diff)

# Save constants of integration to transform the forecasts
# back to their original scale after estimation
const_int <- sapply(ts_list, head, n = 1)
const_int <- log(const_int)

# Create Boxplots for the time series
make_Boxplot <- function(list) {

        for (i in seq_len(length(list))) {

                pdf(file = paste0("./img/raw_series/boxplots/",
                                  names_complete[i],
                                  "_boxPlot"), width = 4)
                boxplot(list[[i]], main = paste0("Boxplot for ", names_complete[i]),
                        xlab = NULL, range = 2)
                abline(h = 0)
                dev.off()

        }

}

# Remove the extreme outliers which lie outside the .01 and.99 quantile
for (i in seq_len(length(diff_list))) {

        range  <- c(0.01, 0.99)
        q <- quantile(ts_list[[i]], probs = range)
        q_diff <- quantile(diff_list[[i]], probs = range)
        diff_list[[i]] <- diff_list[[i]][diff_list[[i]] > q_diff[1] & diff_list[[i]] < q_diff[2]]
        ts_list[[i]] <- ts_list[[i]][ts_list[[i]] > q[1] & ts_list[[i]] < q[2]]
}
rm(q, q_diff, range)

# Create a training and testing set, based on a user specified holdout peroid
make_training_and_testing_sets <- function(tsORListOfTs, out.sample = 30) {

        length_of_series <- sapply(tsORListOfTs, length)

        training.set <- list()
        for (i in seq_len(length(tsORListOfTs))) {

                training.set[[i]] <- tsORListOfTs[[i]][1:(length_of_series[i] - out.sample)]

        }
        names(training.set) <- names_complete

        testing.set <- mclapply(tsORListOfTs, tail, n = out.sample)

        return(list(train = training.set, test = testing.set))
}


# Remove unused Indices and variables
rm(i, namesInvesting, namesYahoo, indicesInvesting, indicesYahoo)
