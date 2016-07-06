#return the value of the function to be minimized
#as well as the gradient w.r.t. index-th observation
targetMPM = function(index, score, uncertainty, adherence, data, mu, sigma){
  source('rank2count.R')
  
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