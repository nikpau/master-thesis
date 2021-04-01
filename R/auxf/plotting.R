# Plot and save function for forecast classes
plot_save_forecast <- function(forecast, testset, window, path) {
  
  
  for (i in seq_len(length(forecast))){
    
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
    
    pdf(file = paste0(path, "/", names_complete[i],"_exp_sm.pdf"), 
        height = 5)
    # Plot series
    
    plot(
      series, 
      type = "b",
      lty = 1,
      pch = 18,
      frame = F,
      lwd = 2,
      xlab = "Index",
      ylab = "Series",
      col = "#002B36",
      ylim = c(miny(series,test), maxy(series,test)),
      xlim = c(length(forecast[[i]]$x) - window, 
               length(forecast[[i]]$x) + length(test)),
      main = paste0(names_complete[i], " ", length(test), " days ahead")
    )
    abline(v = length(forecast[[i]]$x), lwd  = 2)
    
    add_lines(forc, trans_forc, test, trans_true)
    
    #CI
    lines(ui, type = "b", col = "#B58900", pch = 2)
    lines(li, type = "b", col = "#B58900", pch = 2)
    grid()
    plot_legend()
    dev.off()
  }
  cat("Done \n")
  
}
# Find minimum y value of test set and series
miny <- function(series, test) {
  if (min(series) < min(test))
    min(series)
  else
    min(test)
}

# Find minimum y value of test set and series
maxy <- function(series, test) {
  if (max(series) > max(test))
    max(series)
  else
    max(test)
}

# Add forecasted values as lines to the plot
add_lines <- function(forc, trans_forc, true, trans_true) {
  # Forecast
  lines(forc, type = "b", lwd = 2, col = "#859900", pch = 4)
  lines(trans_forc, type = "l", lty = 4, lwd = 2, col = "#859900")
  # True model
  lines(true, type = "b", lwd = 2, col = "#268BD2", pch = 18)
  lines(trans_true, type = "l", lty = 4, lwd = 2, col = "#268BD2")
}

# Legend
plot_legend <- function() {
  legend("topleft", legend = c("Series", "True", "Predict"),
         col = c("#002B36","#268BD2","#859900"),
         lty = c(1,4,4),
         pch = c(18,18,4),
         bg = "white",
         lwd = 2)
}