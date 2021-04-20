# Plot and save function for forecast classes
plot_save_forecast <- function(forecast = NULL, testset, window, path, seriesWithForc = NULL, n.ahead = NULL) {
  
  
  for (i in seq_len(length(testset))){
    
    if (is.null(seriesWithForc)) {
      single_fit <- forecast[[i]]
      
      series <- tail(forecast[[i]]$x, window)
      forc <- forecast[[i]]$mean
      test <- testset[[i]]
      test <- ts(test, start = length(forecast[[i]]$x)+1, 
                 end = length(forecast[[i]]$x) + length(test))
      trans_true <- ts(c(tail(series,1), head(test,1)), 
                       start = tail(index(series),1),
                       end = head(index(test),1))
      trans_forc <- ts(c(tail(series,1), head(forc,1)), 
                       start = tail(index(series),1),
                       end = head(index(forc),1))
      ui <- forecast[[i]]$upper[,2]
      li <- forecast[[i]]$lower[,2]
      
      length_series <- length(forecast[[i]]$x)
      
      ending <- "_forecast.pdf"
    }
    else {
      
      series_raw <- series <- ts(head(seriesWithForc[[i]], -(n.ahead)), start = 1)
      forc <- tail(seriesWithForc[[i]], n.ahead)
      test <- testset[[i]]
      
      test <- ts(test, start = length(series)+1, 
                 end = length(series) + length(test))
      forc <- ts(forc, start = length(series)+1, 
                 end = length(series) + length(forc))
      trans_true <- ts(c(tail(series,1), head(test,1)), 
                       start = tail(index(series),1),
                       end = head(index(test),1))
      trans_forc <- ts(c(tail(series,1), head(forc,1)), 
                       start = tail(index(series),1),
                       end = head(index(forc),1))
      
      series <- tail(series, window)
      
      length_series <- length(series_raw)
      
      ending <- "_forecast.pdf"
    }
    
    pdf(file = paste0(path, "/", names_complete[i],ending), 
        height = 5, family = "Courier")
    # Plot series
    
    plot(
      series, 
      type = "b",
      lty = 1,
      pch = 18,
      frame = F,
      lwd = 3,
      xlab = "Index",
      ylab = "Series",
      col = "#002B36",
      ylim = c(miny(series,test), maxy(series,test)),
      xlim = c(length_series - window, 
               length_series + length(test)),
      main = paste0(names_complete[i], " ", length(test), " days ahead")
    )
    abline(v = length_series, lwd  = 2)
    
    add_lines(forc, trans_forc, test, trans_true)
    
    #CI
    if (is.null(seriesWithForc)) {
    lines(ui, type = "b", col = "#B58900", pch = 2)
    lines(li, type = "b", col = "#B58900", pch = 2)
    }
    grid()
    plot_legend()
    dev.off()
  }
  cat("Done \n")
  
}
# Find minimum y value of test set and series
miny <- function(series, test) {
  if (min(series, na.rm = T) < min(test, na.rm = T))
    min(series, na.rm = T)
  else
    min(test, na.rm = T)
}

# Find minimum y value of test set and series
maxy <- function(series, test) {
  if (max(series, na.rm = T) > max(test, na.rm = T))
    max(series, na.rm = T)
  else
    max(test, na.rm = T)
}

# Add forecasted values as lines to the plot
add_lines <- function(forc, trans_forc, true, trans_true) {
  # Forecast
  lines(forc, type = "b", lwd = 3, col = "#859900", pch = 20)
  lines(trans_forc, type = "l", lty = 4, lwd = 3, col = "#859900")
  # True model
  lines(true, type = "b", lwd = 3, col = "#268BD2", pch = 18)
  lines(trans_true, type = "l", lty = 4, lwd = 3, col = "#268BD2")
}

# Legend
plot_legend <- function() {
  legend("topleft", legend = c("Series", "True", "Predict"),
         col = c("#002B36","#268BD2","#859900"),
         lty = c(1,4,4),
         pch = c(18,18,20),
         bg = "white",
         lwd = 3)
}

# Plotting function for the neural networks
plot_save_nn <- function(nnPrediction, window, type) {
  
  for (i in seq_len(length(nnPrediction))) {
    
    series <- nnPrediction[[i]]$train %>% 
      ts(start = 1, end = length(nnPrediction[[i]]$train)) %>% 
      tail(window)
    
    forc <-  nnPrediction[[i]]$forc %>% 
      ts(start = length(nnPrediction[[i]]$train) + 1, 
         end = length(nnPrediction[[i]]$train) + 1 + length(nnPrediction[[i]]$forc))
    
    test <-  nnPrediction[[i]]$test %>% 
      ts(start = length(nnPrediction[[i]]$train) + 1, 
         end = length(nnPrediction[[i]]$train) + 1 + length(nnPrediction[[i]]$forc))
    
    trans_true <- ts(c(tail(series,1), head(test,1)), 
                     start = length(nnPrediction[[i]]$train),
                     end = length(nnPrediction[[i]]$train) +1)
    trans_forc <- ts(c(tail(series,1), head(forc,1)), 
                     start = length(nnPrediction[[i]]$train),
                     end = length(nnPrediction[[i]]$train) + 1)
    
    if (type == "lstm") {
      path <- "./img/ann_forecast/lstm/"
      ending <- paste0(names_complete[i], "_lstm_forecast.pdf")
    }
    else {
      path <- "./img/ann_forecast/mlp/"
      ending <- paste0(names_complete[i], "_mlp_forecast.pdf")
    }
    
    pdf(file = paste0(path,ending), 
        height = 5, family = "Courier")
    
    plot(
      series, 
      type = "b",
      lty = 1,
      pch = 18,
      frame = F,
      lwd = 3,
      xlab = "Index",
      ylab = "Series",
      col = "#002B36",
      ylim = c(miny(series,test), maxy(series,test)),
      xlim = c(length(nnPrediction[[i]]$train) - window, 
               length(nnPrediction[[i]]$train) + length(nnPrediction[[i]]$forc)),
      main = paste0(names_complete[i], " ", length(nnPrediction[[i]]$forc), " days ahead")
    )
    
    abline(v = length(nnPrediction[[i]]$train), lwd  = 3)
    
    add_lines(forc, trans_forc, test, trans_true)
    
    grid()
    plot_legend()
    
    dev.off()
    
  }
  
}
