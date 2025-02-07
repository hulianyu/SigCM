source('jackstraw_cluster.R')
# Dataset filenames and row names
filename <- c('lenses', 'lung-cancer', 'soybean-small', 'zoo', 'dna-promoter',
              'hayes-roth', 'lymphography', 'heart-disease', 'solar-flare', 'primary-tumor',
              'dermatology', 'house-votes', 'balance-scale', 'credit-approval', 'breast-cancer-wisconsin',
              'mammographic-mass', 'tic-tac-toe', 'car')

rowNames <- c('Ls', 'Lc', 'So', 'Zo', 'Ps', 'Hr', 'Ly', 'Hd', 'Sf',
              'Pt', 'De', 'Hv', 'Bs', 'Ca', 'Bc', 'Mm', 'Tt', 'Ce')

# Loop through each dataset
for (I in 1:18) {
  
  # Load the current dataset
  dataset_name <- filename[I]
  X_data <- read.table(paste0(dataset_name, '-onehot.txt'))  # Load the txt file
  X <- X_data
  N <- nrow(X)  # Number of samples
  
  # Convert X to matrix if it's not already a matrix
  X <- as.matrix(X)  # Ensures X is a matrix
  
  # Get the current GT partition
  all_pi <- read.table(paste0(dataset_name, '-GT.txt'))  # Load the txt file
  
  # Initialize a vector to store all p-values
  all_pvs <- numeric(0)  # Start with an empty vector
  
  total_time <- 0  # Start with zero total time
  
  # Run the loop for the specified number of iterations
  for (run in 1:1) {
    pi <- all_pi[, run]
    K <- max(pi)  # Number of clusters
    
    # Compute the cluster centers (mean of each cluster)
    my_centers <- matrix(NA, nrow = K, ncol = ncol(X))  # Create an empty matrix for cluster centers
    for (i in 1:K) {
      my_centers[i, ] <- colMeans(X[pi == i, , drop = FALSE])  # Compute the mean for each cluster
    }
    
    # Convert my_centers to matrix (in case it's not already a matrix)
    my_centers <- as.matrix(my_centers)
    
    start_time <- Sys.time()  # Capture the start time
    
    # Call the jackstraw_cluster function
    result <- jackstraw_cluster(
      dat = X,         # Data matrix
      k = K,           # Number of clusters
      cluster = pi,    # Cluster labels
      centers = my_centers  # Cluster centers
    )
    
    # Capture the end time and calculate the duration
    end_time <- Sys.time()  # Capture the end time
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))  # Time difference in seconds
    
    # Get the p-values and append to the all_pvs vector
    all_pvs <- c(all_pvs, result$p.F)
    # Accumulate the total execution time for all runs
    total_time <- total_time + execution_time
    cat("Processing dataset", rowNames[I], "- Run", run, "of", RT, "\n")
  }
  
  # Save the all_pvs to a file with the corresponding name
  filename_pvs <- paste0('GT_Jackstraw_all_pvs_', rowNames[I], '.txt')
  # write.table(all_pvs, file = filename_pvs, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  
  filename_time <- paste0('GT_Jackstraw_runtime_', rowNames[I], '.txt')
  # write.table(total_time, file = filename_time, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Output progress
  cat("Finished processing dataset: ", rowNames[I], "\n")
  cat("Total time for all runs:", total_time, "seconds\n")
}
