#Plackett-Luce model

#calculate the target function to minimize in the Plackett-Luce model
#as well as the gradient for each observation (in Jacobian-like style)

targetPL = function(score, data, mu, sigma){
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
  
  #J-matrix, each row is the gradient
  J = matrix(0, nobs, nvar)
  #initialize
  inv_sigma = solve(sigma)
  target_value = as.numeric(0.5 * (t(score - mu) %*% inv_sigma %*% (score - mu)))
  
  #loop over all observations
  for(i in 1:nobs){
    
    #calculate the ranking of the form A>B>C...
    ranks = data[i, ][data[i, ] != 0]
    ranking = as.numeric(names(sort(ranks)))
    
    #the length of i-th observation
    nrank = length(ranks)
    
    #initialize the gradient
    gradient = J[i, ]
    gradient[1:nvar] = 1 / nobs * inv_sigma %*% (score - mu)
    
    #loop over all pairwise comparisons
    for(j in 1:(nrank - 1)){
      
      sum_temp = 0
      sum_temp2 = 0
      
      for(k in (j + 1):nrank){
        
        #ranking[j] wins over ranking[k]
        win = ranking[j]
        lose = ranking[k]
        
        sum_temp = sum_temp + exp(-(score[win] - score[lose]))
        sum_temp2 = sum_temp2 + exp(-(score[win] - score[lose])) *
          (-score[win] + score[lose])
        
      }
      
      #update the value of the target function
      target_value = target_value + log(1 + sum_temp)
      
      #update the gradient w.r.t. score
      after_j = ranking[(j + 1):nrank] #varieties ranked after variety j
      gradient[after_j] = gradient[after_j] + (1 / (1 + sum_temp)) * 
        exp(-(score[win] - score[after_j]))
      gradient[win] = gradient[win] - (1/ (1 + sum_temp)) * sum_temp
      
      
    }
    
    J[i, ] = gradient
    
  }
  
  return(list(value = target_value, Jacobian = J))
  
  
}

