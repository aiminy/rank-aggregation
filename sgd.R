#transform the ranking data into a count matrix
rank2count = function(data){
  
  #inpute data is a vector,
  #where i-th element in the vector is the rank assigned to the i-th item
  #the entry where item not ranked is replaced by 0
  
  #assume there are no ties
  
  nitem = length(data)
  names(data) = 1:length(data) #assign labels
  
  rank = data[data != 0]
  nitem_ranked = length(rank)
  C = matrix(0, nitem, nitem)
  
  #loop over all pairwise comparisons
  for(i in 1:(nitem_ranked - 1)){
    
    for(j in (i + 1):nitem_ranked){
      
      entry = rank[i] - rank[j]
      item_i = as.numeric(names(rank)[i])
      item_j = as.numeric(names(rank)[j])
      
      if(entry > 0){ #item_j wins over item_i
        C[item_j, item_i] = entry
      } else{  #item_i wins over item_j
        C[item_i, item_j] = -entry
      }
      
    }
  }
  
  return(C)
  
}




######TARGET FUNCTIONS
#return the value of the function to be minimized
#as well as the gradient w.r.t. index-th observation
targetThurs = function(index, score, adherence, data, mu, sigma){
  #let m be the number of varieties,
  #let n be the number of farmers.
  #data is an n*m matrix,
  #data(i, j) represents the rank of variety i by farmer j
  #the entry where varieties are not included is 0
  
  #(mu, sigma) are parameters of the normal prior
  
  nobs = nrow(data)
  nvar = ncol(data)
  colnames(data) = 1:nvar #assign labels to varieties
  
  #the first nvar element is the gradient for score
  #the last nobs element is the gradient for adherence
  gradient = rep(0, nvar + nobs)
  #initialize
  inv_sigma = solve(sigma)
  target_value = as.numeric(0.5 * (t(score - mu) %*% inv_sigma %*% (score - mu)))
  gradient[1:nvar] = 1 / nobs * inv_sigma %*% (score - mu)
  
  #loop over all observations
  for(i in 1:nobs){
    
    #calculate the ranking of the form A>B>C...
    ranks = data[i, ][data[i, ] != 0]
    ranking = as.numeric(names(sort(ranks)))
    
    #the length of i-th observation
    nrank = length(ranks)
    
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
        
        if(i == index){
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
    }
    
  }
  
  return(list(value = target_value, gradient = gradient))
  
}






#return the value of the function to be minimized
#as well as the gradient w.r.t. index-th observation
targetBT = function(index, score, adherence, data, mu, sigma){
  
  #let m be the number of varieties,
  #let n be the number of farmers.
  #data is an n*m matrix,
  #data(i, j) represents the rank of variety i by farmer j
  #the entry where varieties are not included is 0
  
  #(mu, sigma) are parameters of the normal prior
  
  nobs = nrow(data)
  nvar = ncol(data)
  colnames(data) = 1:nvar #assign labels to varieties
  
  #the first nvar element is the gradient for score
  #the last nobs element is the gradient for adherence
  gradient = rep(0, nvar + nobs)
  
  #initialize
  inv_sigma = solve(sigma)
  target_value = as.numeric(0.5 * (t(score - mu) %*% inv_sigma %*% (score - mu)))
  gradient[1:nvar] = 1 / nobs * inv_sigma %*% (score - mu)
  
  #loop over all observations
  for(i in 1:nobs){
    
    #calculate the ranking of the form A>B>C...
    ranks = data[i, ][data[i, ] != 0]
    ranking = as.numeric(names(sort(ranks)))
    
    #the length of i-th observation
    nrank = length(ranks)
    
    #loop over all pairwise comparisons
    for(j in 1:(nrank - 1)){
      
      for(k in (j + 1):nrank){
        
        #ranking[j] wins over ranking[k]
        win = ranking[j]
        lose = ranking[k]
        
        exp_term = exp(-adherence[i] * (score[win] - score[lose]))
        
        #update the value of the target function
        target_value = target_value + log(1 + exp_term)
        
        if(index == i){
          #update the gradient w.r.t. score
          gradient[win] = gradient[win] - adherence[i] * exp_term / (1 + exp_term)
          gradient[lose] = gradient[lose] + adherence[i] * exp_term / (1 + exp_term)
          
          #update the gradient w.r.t. adherence
          gradient[nvar + i] = gradient[nvar + i] +
            (-score[win] + score[lose]) * exp_term / (1 + exp_term)
          
        }
        
      }
    }
    
    
  }
  
  return(list(value = target_value, gradient = gradient))
}



#return the value of the function to be minimized
#as well as the gradient w.r.t. index-th observation
targetPL = function(index, score, adherence, data, mu, sigma){
  #let m be the number of varieties,
  #let n be the number of farmers.
  #data is an n*m matrix,
  #data(i, j) represents the rank of variety i by farmer j
  #the entry where varieties are not included is 0
  
  #(mu, sigma) are parameters of the normal prior
  
  nobs = nrow(data)
  nvar = ncol(data)
  colnames(data) = 1:nvar #assign labels to varieties
  
  #the first nvar element is the gradient for score
  #the last nobs element is the gradient for adherence
  gradient = rep(0, nvar + nobs)
  #initialize
  inv_sigma = solve(sigma)
  target_value = as.numeric(0.5 * (t(score - mu) %*% inv_sigma %*% (score - mu)))
  gradient[1:nvar] = 1 / nobs * inv_sigma %*% (score - mu)
  
  #loop over all observations
  for(i in 1:nobs){
    
    #calculate the ranking of the form A>B>C...
    ranks = data[i, ][data[i, ] != 0]
    ranking = as.numeric(names(sort(ranks)))
    
    #the length of i-th observation
    nrank = length(ranks)
    
    
    #loop over all pairwise comparisons
    for(j in 1:(nrank - 1)){
      
      sum_temp = 0
      sum_temp2 = 0
      
      for(k in (j + 1):nrank){
        
        #ranking[j] wins over ranking[k]
        win = ranking[j]
        lose = ranking[k]
        
        sum_temp = sum_temp + exp(-adherence[i] * (score[win] - score[lose]))
        sum_temp2 = sum_temp2 + exp(-adherence[i] * (score[win] - score[lose])) *
          (-score[win] + score[lose])
        
      }
      
      #update the value of the target function
      target_value = target_value + log(1 + sum_temp)
      
      if(index == i){
        #update the gradient w.r.t. score
        after_j = ranking[(j + 1):nrank] #varieties ranked after variety j
        gradient[after_j] = gradient[after_j] + (adherence[i] / (1 + sum_temp)) * 
          exp(-adherence[i] * (score[win] - score[after_j]))
        gradient[win] = gradient[win] - (adherence[i] / (1 + sum_temp)) * sum_temp
        
        #update the gradient w.r.t. adherence
        gradient[nvar + i] = gradient[nvar + i] + (1 / (1 + sum_temp)) * sum_temp2
        
      }
      
    }
    
    
  }
  
  return(list(value = target_value, gradient = gradient))
  
  
}





#return the value of the function to be minimized
#as well as the gradient w.r.t. index-th observation
targetMPM = function(index, score, uncertainty, adherence, data, mu, sigma){
  
  #let m be the number of varieties,
  #let n be the number of farmers.
  #data is an n*m matrix,
  #data(i, j) represents the rank of variety i by farmer j
  #the entry where varieties are not included is 0
  
  #(mu, sigma) are parameters of the normal prior
  
  nobs = nrow(data)
  nvar = ncol(data)
  colnames(data) = 1:nvar #assign labels to varieties
  
  #the first nvar elements are the gradient for score,
  #the next nvar elements are the gradient for uncertainty,
  #the last nobs elements are the gradient for adherence
  gradient = rep(0, (nvar + nvar + nobs))
  
  #initialize
  inv_sigma = solve(sigma)
  target = as.numeric(0.5 * (t(score - mu) %*% inv_sigma %*% (score - mu)))
  gradient[1:nvar] = 1 / nobs * inv_sigma %*% (score - mu)
  ranks = data[index, ][data[index, ]!= 0] #the index-th observation
  ranking = as.numeric(names(sort(ranks))) #the index-th ranking
  nrank = length(ranking)
  
  #pre-compute numbers (pairwise) for later computations
  pair_score = outer(score, score, '-')
  pair_uncer = outer(uncertainty, uncertainty, '+')
  frac_term = pair_score / pair_uncer
  
  ##calculate the target value
  vec1 = rep(0, nobs)
  sum_count = rep(0, nobs)
  for(i in 1:nobs){ #loop over all observations
    count = rank2count(data[i, ]) #count matrix
    sum_count[i] = sum(count)
    
    temp_res = adherence[i] * frac_term
    
    #update target value
    target = target - sum(count * temp_res)
    
    temp_res = exp(temp_res)
    diag(temp_res) = 0
    
    #store a specific temp_res for later computations
    if(i == index){
      exp_term = temp_res
      mycount = count
    }
    
    vec1[i] = sum(temp_res)
  }
  #update target value
  target = target + sum(sum_count * log(vec1))
  
  
  ##compute gradient
  vec2 = rep(0, nvar)
  vec3 = rep(0, nvar)
  for(i in 1:nvar){
    temp_res = (exp_term[i, ] - exp_term[, i]) * adherence[index] 
    vec2[i] = sum(temp_res / pair_uncer[i, ])
    vec3[i] = sum(temp_res * pair_score[i, ] / (pair_uncer[i, ])^2)
    
    #update gradient w.r.t. score
    gradient[i] = gradient[i] + sum(mycount[, i] * adherence[index] / pair_uncer[i, ])
    gradient[i] = gradient[i] - sum(mycount[i, ] * adherence[index] / pair_uncer[i, ])
    
    #update gradient w.r.t. uncertainty
    gradient[nvar + i] = gradient[nvar + i] + sum(mycount[, i]
                                                  * adherence[index] * pair_score[, i]
                                                  / (pair_uncer[i, ])^2)
    gradient[nvar + i] = gradient[nvar + i] + sum(mycount[i, ]
                                                  * adherence[index] * pair_score[i, ] 
                                                  / (pair_uncer[i, ])^2)
    
  }
  #update gradient w.r.t. score
  gradient[1:nvar] = gradient[1:nvar] + (1 / vec1[index]) * vec2 * sum_count[index]
  
  #update gradient w.r.t. uncertainty
  gradient[(nvar + 1):(2 * nvar)] = gradient[(nvar + 1):(2 * nvar)] - 
    (1 / vec1[index]) * vec3 * sum_count[index]
  
  #update gradient w.r.t. adherence
  gradient[2 * nvar + index] = gradient[2 * nvar + index] - sum(mycount * frac_term)
  temp_num = sum(exp_term * frac_term)
  gradient[2 * nvar + index] = gradient[2 * nvar + index] + 
    (1 / vec1[index]) * temp_num * sum_count[index]
  
  return(list(value = target, gradient = gradient))
  
}






#####SGD FUNCTIONS
sgdThurs = function(data, mu, sigma, rate, maxiter = 1000, tol = 1e-9, start, decay){
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
  #the first nvar element is the score
  #the last nobs element is the adherence
  param = start
  target = rep(0, niter)
  
  flag = TRUE
  #loop until the convergence criteria are met
  while(flag){
    
    for(i in 1:nobs){
      
      niter = niter + 1
      
      score_temp = param[1:nvar]
      adherence_temp = param[(nvar + 1):(nvar + nobs)]
      #evaluate the log-posterior as well as the gradient
      #only used for small dataset (where we want to decide the learning rate)
      #if used for big dataset, where we don't want to
      #evaluate log-posterior everytime, the function should be modified
      res_temp = targetThurs(i, score_temp, adherence_temp, data, mu, sigma)
      
      #store the value of the target function
      target[niter] = res_temp[[1]]
      
      #extract the gradient
      gradient = res_temp[[2]]
      
      #update the parameters
      param = param - rate * gradient
      
      #if the updated adherence is negative, change it to 0.5
      param[(nvar + 1):(nvar + nobs)][param[(nvar + 1):(nvar + nobs)] < 0] = 0.5
      
      
      #check the convergence criteria: square of the change of target values
      if(niter > 1){
        if((target[niter] - target[niter - 1]) ^ 2 < tol | niter > maxiter){
          flag = FALSE
          break
        }
        
        #update learning rate if the target value don't decrease
        if((target[niter - 1] - target[niter]) / target[niter - 1] < 0){
          rate = rate / decay
        }
      }
      
      
    }
    
  }
  
  return(list(value = target, niter = niter, score = param[1:nvar],
              adherence = param[(nvar + 1):(nvar + nobs)]))
  
}




sgdBT = function(data, mu, sigma, rate, maxiter = 1000, tol = 1e-9, start, decay){
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
  #the first nvar element is the score
  #the last nobs element is the adherence
  param = start
  target = rep(0, niter)
  
  flag = TRUE
  #loop until the convergence criteria are met
  while(flag){
    
    for(i in 1:nobs){
      niter = niter + 1
      score_temp = param[1:nvar]
      adherence_temp = param[(nvar + 1):(nvar + nobs)]
      #evaluate the log-posterior as well as the gradient
      #only used for small dataset (where we want to decide the learning rate)
      #if used for big dataset, where we don't want to
      #evaluate log-posterior everytime, the function should be modified
      res_temp = targetBT(i, score_temp, adherence_temp, data, mu, sigma)
      
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
          rate = rate / decay
        }
      }
      
      
    }
    
  }
  
  return(list(value = target, niter = niter, score = param[1:nvar],
              adherence = param[(nvar + 1):(nvar + nobs)]))
}




sgdPL = function(data, mu, sigma, rate, maxiter = 1000, tol = 1e-9, start, decay){
  
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
  #the first nvar element is the score
  #the last nobs element is the adherence
  param = start
  target = rep(0, niter)
  
  flag = TRUE
  #loop until the convergence criteria are met
  while(flag){
    
    for(i in 1:nobs){
      niter = niter + 1
      score_temp = param[1:nvar]
      adherence_temp = param[(nvar + 1):(nvar + nobs)]
      #evaluate the log-posterior as well as the gradient
      #only used for small dataset (where we want to decide the learning rate)
      #if used for big dataset, where we don't want to
      #evaluate log-posterior everytime, the function should be modified
      res_temp = targetPL(i, score_temp, adherence_temp, data, mu, sigma)
      
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
          rate = rate / decay
        }
      }
      
      
    }
    
  }
  
  return(list(value = target, niter = niter, score = param[1:nvar],
              adherence = param[(nvar + 1):(nvar + nobs)]))
  
}




sgdMPM = function(data, mu, sigma, rate, maxiter = 1000, tol = 1e-9, start, decay){

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
          rate = rate / decay
        }
      }
      
      
    }
    
  }
  
  return(list(value = target, niter = niter, score = param[1:nvar],
              uncertainty = param[(nvar + 1):(2 * nvar)],
              adherence = param[(2 * nvar + 1):(2 * nvar + nobs)]))
  
}
