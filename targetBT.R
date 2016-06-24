#calculate the target function to minimize in the B-T model
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
  res = as.numeric(0.5 * (t(score - mu) %*% solve(sigma) %*% (score - mu)))
  
  #loop over all observations
  for(i in 1:nobs){
    
    #calculate the ranking of the form A>B>C...
    ranks = data[i, ][data[i, ] != 0]
    ranking = as.numeric(names(sort(ranks)))
    
    #the length of i-th observation
    nrank = length(ranks)
    
    #calculate likelihood
    for(j in 1:(nrank - 1)){
      
      for(k in (j + 1):nrank){
        
        res = res + log(1 + exp(-score[ranking[j]] + score[ranking[k]]))
        
      }
      
    }
  }
  
  return(res)
}



