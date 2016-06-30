#Plackett-Luce model modified to accommodate the 'adherence parameter'

#calculate the target function to minimize in the Plackett-Luce model
#as well as the gradient for each observation (in Jacobian-like style)

targetPL2 = function(score, adherence, data, mu, sigma){
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
      
      #update the gradient w.r.t. score
      after_j = ranking[(j + 1):nrank] #varieties ranked after variety j
      gradient[after_j] = gradient[after_j] + (adherence[i] / (1 + sum_temp)) * 
               exp(-adherence[i] * (score[win] - score[after_j]))
      gradient[win] = gradient[win] - (adherence[i] / (1 + sum_temp)) * sum_temp
      
      #update the gradient w.r.t. adherence
      gradient[nvar + i] = gradient[nvar + i] + (1 / (1 + sum_temp)) * sum_temp2
      
    }
    
    J[i, ] = gradient
    
  }
  
  return(list(value = target_value, Jacobian = J))
  
  
}



#test
score = c(1,2,3)
adherence = 0.8
data = t(c(3,2,1)); data
mu = c(1,1,1)
sigma = diag(1, 3)
res1 = targetPL2(score, adherence, data, mu, sigma)

#test w.r.t. score
#change the score a little bit
score2 = c(1.00, 2.001, 3.00)
res2 = targetPL2(score2, adherence, data, mu, sigma)
#test using finite difference
(res2[[1]] - res1[[1]]) / 0.001
res1[[2]][2]
res2[[2]][2]

#test w.r.t. adherence
#change the adherence a little bit
adherence2 = 0.801
res3 = targetPL2(score, adherence2, data, mu, sigma)
(res3[[1]] - res1[[1]]) / 0.001
res1[[2]][[4]]
res3[[2]][[4]]













