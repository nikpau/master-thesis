# Plot and save function for forecast classes
plot_save_forecast <- function(forecast, test, window, path) {
  
  
  for (i in seq_len(length(forecast))){
    
    single_fit <- forecast[[i]]
    
    series <- tail(forecast[[i]]$x, window)
    forc <- forecast[[i]]$mean
    test <- test[[i]]
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
    miny <- function() {if(min(series) < min(test)) min(series) else min(test)}
    maxy <- function() {if(max(series) > max(test)) max(series) else max(test)}
    
    pdf(file = paste0(path, "/", names_complete[i],"_exp_sm.pdf"), 
        height = 5,
        family = "Times")
    
    plot(
      series, 
      type = "l", 
      lwd = 2,
      xlab = "Index",
      ylab = "Series",
      ylim = c(floor(miny()), ceiling(maxy())),
      xlim = c(length(forecast[[i]]$x) - window, 
               length(forecast[[i]]$x) + length(test)),
      main = paste0(names_complete[i], " ", length(test), " days ahead")
    )
    abline(v = length(forecast[[i]]$x))
    lines(forc, type = "l", lwd = 2, col = "blue2")
    lines(trans_forc, type = "l", lty = 2, lwd = 2, col = "blue2")
    lines(test, type = "l", lwd = 2, col = "green4")
    lines(trans_true, type = "l", lty = 2, lwd = 2, col = "green4")
    #CI
    lines(ui, type = "b", col = "red4")
    lines(li, type = "b", col = "red4")
    grid()
    dev.off()
  }
  cat("Done \n")
  
}