sgdBT = function(data, mu, sigma, rate, maxiter = 1000, tol = 1e-9, start){
  source('targetBT.R')
  #let m be the number of varieties,
  #let n be the number of farmers.
  #data is an n*m matrix,
  #data(i, j) represents the rank of variety i by farmer j
  #the entry where varieties are not included is 0
  
  #rate is the step size/learning rate
  
  nobs = nrow(data)
  nvar = ncol(data)
  colnames(data) = 1:nvar #assign labels to varieties
  inv_sigma = solve(sigma)
  
  #initialize
  niter = 0
  score = start
  target = rep(0, niter)
  
  flag = TRUE
  #loop until meet the convergence criteria
  while(flag){
    
    #shuffle the data
    #perm = sample(nobs)
    for(i in 1:nobs){ 
      
      niter = niter + 1
      #evaluate the likelihood as well as the gradient
      #only used for small dataset (where we want to decide the learning rate)
      #if used for big dataset, where we don't want to
      #evaluate likelihood everytime, the function should be modified
      res_temp = targetBT(score, data, mu, sigma)
      #store the value of the target function
      target[niter] = res_temp[[1]]
      
      #extract the gradient
      gradient = res_temp[[2]][i, ]
      
      #update score
      score = score - rate * gradient
      
      #check the convergence criteria: square of the change of target values
      if(niter > 1){
        if((target[niter] - target[niter - 1]) ^ 2 < tol | niter > maxiter){
          flag = FALSE
          break
        }
      }
    
    }
  }
  
  return(list(value = target, niter = niter, score = score))
}
