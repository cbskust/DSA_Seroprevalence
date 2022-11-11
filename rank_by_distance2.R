
rank_by_distance2 <- function(frame, metric='mode'){
  
  # Define
  parameter.draws <- frame
  
  # Determine metric to define 'center'
  if (metric == 'mean'){
    # Mean of parameters
    pt <- apply(parameter.draws, 2, mean)
    
  } else if (metric == 'median'){
    # Median of parameters
    pt <- apply(parameter.draws, 2, median)
    
  } else if (metric == 'mode'){
    # Mode of parameters
    pt <- apply(parameter.draws, 2, mlv, method = "lientz", bw = 0.2)
    
  } else {
    # Throw error
    stop("Metric must be 'mean', 'median', or 'mode'.")
    
  }
  
  # vector of weights 
  weights = 1/apply(parameter.draws, 2, var)
  # Vector of distances
  distances <- apply(parameter.draws, 1, function(x) weights%*%(x - pt) ^ 2)
  
  # Append to data and sort
  parameter.draws$distance <- distances
  
  # Rank and re-index
  parameter.draws.ranked <- parameter.draws[order(parameter.draws$distance), ]
  row.names(parameter.draws.ranked) <- NULL
  
  # Return sorted results and the point used ranking
  return(list(parameter.draws.ranked, pt))
}

################################################################################