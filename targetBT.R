#calculate the target function to minimize in the B-T model
#as well as the gradient for each observation (in Jacobian-like style)

targetBT = function(score, data, mu, sigma){
  #browser()
  
  #let m be the number of varieties,
  #let n be the number of farmers.
  #data is an n*m matrix,
  #data(i, j) represents the rank of variety i by farmer j
  #the entry where varieties are not included is 0
  
  #(mu, sigma) are parameters of the normal prior
  
  nobs = nrow(data)
  nvar = ncol(data)
  colnames(data) = 1:nvar #assign labels to varieties
  
  #initialize, add the normal term
  inv_sigma = solve(sigma)
  target_value = as.numeric(0.5 * (t(score - mu) %*% inv_sigma %*% (score - mu)))
  #J-matrix, each row is the gradient
  J = matrix(0, nobs, nvar)
  
  
  #loop over all observations
  for(i in 1:nobs){
    
    #calculate the ranking of the form A>B>C...
    ranks = data[i, ][data[i, ] != 0]
    ranking = as.numeric(names(sort(ranks)))
    
    #the length of i-th observation
    nrank = length(ranks)
    
    #initialize the gradient
    gradient = 1 / nobs * inv_sigma %*% (score - mu)
    
    
    for(j in 1:(nrank - 1)){
      
      for(k in (j + 1):nrank){
        
        #ranking[j] wins over ranking[k]
        win = ranking[j]
        lose = ranking[k]
        exp_term = exp(-score[win] + score[lose])
        
        #update the value of the target function
        target_value = target_value + log(1 + exp_term)
        
        #update gradient
        gradient[win] = gradient[win] - exp_term / (1 + exp_term)
        gradient[lose] = gradient[lose] + exp_term / (1 + exp_term)
          
      }
      
    }
    
    J[i, ] = gradient
    
    
  }
  
  return(list(value = target_value, Jacobian = J))
}

