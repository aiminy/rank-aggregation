sgdMPM = function(data, mu, sigma, rate, maxiter = 1000, tol = 1e-9, start){
  source('targetMPM.R')
  #let m be the number of varieties,
  #let n be the number of farmers.
  #data is an n*m matrix,
  #data(i, j) represents the rank of variety i by farmer j
  #the entry where varieties are not included is 0
  
  #(mu, sigma) is the parameters of the normal prior
  
  #rate is the starting step size/learning rate
  #start is the starting value of the parameter, ordered as 
  #c(score, uncertainty, adherence)
  
  nobs = nrow(data)
  nvar = ncol(data)
  colnames(data) = 1:nvar #assign labels to varieties
  inv_sigma = solve(sigma)
  
  #initialize
  niter = 0
  #the first nvar element is the score
  #the next nvar element is the uncertainty
  #the last nobs element is the adherence
  param = start
  target = rep(0, niter)
  
  flag = TRUE
  #loop until the convergence criteria are met
  while(flag){
    
    for(i in 1:nobs){
      
      niter = niter + 1
      score_temp = param[1:nvar]
      uncertainty_temp = param[(nvar + 1):(2*nvar)]
      adherence_temp = param[(2*nvar + 1):(2*nvar + nobs)]
      #evaluate the log-posterior as well as the gradient
      #only used for small dataset (where we want to decide the learning rate)
      #if used for big dataset, where we don't want to
      #evaluate log-posterior everytime, the function should be modified
      res_temp = targetMPM(i, score_temp, uncertainty_temp, adherence_temp, data, mu, sigma)
      
      #store the value of the target function
      target[niter] = res_temp[[1]]
      
      #extract the gradient
      gradient = res_temp[[2]]
      
      #update the parameters
      param = param - rate * gradient
      
      #check the convergence criteria: square of the change of target values
      if(niter > 1){
        if((target[niter] - target[niter - 1]) ^ 2 < tol | niter > maxiter){
          flag = FALSE
          break
        }
        
        #update learning rate if the target value don't decrease
        if((target[niter - 1] - target[niter]) / target[niter - 1] < 0){
          rate = rate / 1.05
        }
      }
      
      
    }
    
  }
  
  return(list(value = target, niter = niter, score = param[1:nvar],
              uncertainty = param[(nvar + 1):(2 * nvar)],
              adherence = param[(2 * nvar + 1):(2 * nvar + nobs)]))
  
}