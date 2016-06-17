plMM = function(Dmatrix, nitem, tol = 1e-9, maxiter = 1000){
  #Each row in Dmatrix is the ranking data, 
  #the length of each row may vary,
  #entry with no ranking is replaced by 0.
  
  #nitem is the total number of items to rank,
  #originally the item is labeled from 1 to nitem
  
  #note: input data must be proper, items labeled from 1 to nitem
  #the item violating assumption 1 in the paper should be deleted
  
  #initialize
  gamma_old = rep(1/nitem, nitem) #assign equal probabilities
  gamma_new = rep(1/nitem, nitem)
  niter = 0
  nobs = nrow(Dmatrix) #total number of rankings observed
  
  while(1){ #do until meeting the convergence criteria
    niter = niter + 1
    
    for(t in 1:nitem){ #loop over all items
      #number of rankings in which t-th item is ranked higher than last
      #also the numerator in the update equation
      win_t = 0
      #denominator in the update equation
      denom = 0
      
      for(j in 1:nobs){ #loop over all observations
        
        ranking = Dmatrix[j, ][Dmatrix[j, ] != 0] #non-zero entry
        nitem_each = length(ranking) #number of items in each observation
        #test if item t is ranked higher than last
        win_t = win_t + !is.na(match(t ,ranking[1:(nitem_each-1)]))
        
        #test if item t is in the j-th observation
        if(!is.na(match(t, ranking))){
          for(i in 1:(nitem_each - 1)){ #loop over all items in j-th obs
            
            #test if item t is ranked no better than i-th place
            rank_t = match(t, ranking)
            if(!rank_t < i){
              
              partial_sum = 0
              
              for(s in i:nitem_each){
                #which item is ranked in the s-th place in j-th obs
                item_ranked_s = ranking[s]
                partial_sum = partial_sum + gamma_old[item_ranked_s]
                
              }
              
              denom = denom + 1 / partial_sum
              
            }
            
            
          }
          
        }
        
      }
      
      gamma_new[t] = win_t / denom #update t-th parameter
      
      
    }
    #normalize
    gamma_new = gamma_new / sum(gamma_new)
    
    cat('Now at iteration ', niter, '\n')
    
    #test the convergence criteria
    if(sqrt(sum((gamma_new - gamma_old)^2)) < tol | niter > maxiter){
      break
    } else{
      gamma_old = gamma_new
    }
    
  }
  
  if(niter > maxiter){
    cat('MM doesn\'t converge \n')
  } else{
    cat('At iteration ', niter, ', MM converges \n')
    
  }
  
  return(gamma_new)
}

