plot_distribution <- function(X, k) {
  
  # Calculate the average of X
  avg_X <- mean(X, na.rm = TRUE)
  
  # Create a histogram of X
  hist(X, breaks = 50, main = "Distribution of X with Vertical Lines", xlab = "Values", col = "gray", border = "black")
  
  # Add a vertical dashed line at the average of X
  abline(v = avg_X, col = "red", lty = "dashed", lwd = 2)
  
  # Add a vertical dashed line at k
  abline(v = k, col = "blue", lty = "dashed", lwd = 2)
  
  # Add a legend
  legend("topright", legend = c("Average of X", "Value k"), col = c("red", "blue"), lty = "dashed")
}
