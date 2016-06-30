#Thurstone model modified to accommodate the 'adherence parameter'

#calculate the target function to minimize in the B-T model
#as well as the gradient for each observation (in Jacobian-like style)

targetThur2 = function(score, adherence, data, mu, sigma){
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
  #in each row, the first nvar element is the gradient for score
  #the last nobs element is the gradient for adherence
  J = matrix(0, nobs, (nvar + nobs))
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
      
      for(k in (j + 1):nrank){
        
        #ranking[j] wins over ranking[k]
        win = ranking[j]
        lose = ranking[k]
        
        #calculate the term involved with cdf of std. normal
        quant = sqrt(adherence[i]) * (score[win] - score[lose]) / sqrt(2)
        cdf_term = pnorm(quant)
        
        #update the value of the target function
        target_value = target_value - log(cdf_term)
        
        #update the gradient w.r.t. score
        grad_change = (1 / cdf_term) * (1 / sqrt(2 * pi)) * 
          exp(-0.5 * quant^2) * (sqrt(adherence[i]) / sqrt(2))
        gradient[win] = gradient[win] - grad_change
        gradient[lose] = gradient[lose] + grad_change
        
        #update the gradient w.r.t. adherence
        gradient[nvar + i] = gradient[nvar + i] -
          (1 / cdf_term) * (1 / sqrt(2 * pi)) * exp(-0.5 * quant^2) *
          (score[win] - score[lose]) / (sqrt(2) * 2 * sqrt(adherence[i]))
        
        
      }
    }
    
    
    J[i, ] = gradient
    
  }
  
  return(list(value = target_value, Jacobian = J))
  
}
