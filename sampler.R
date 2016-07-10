library(stats) #sampler from multinomial
rMallow = function(sigma, phi){
  #sigma is the reference ranking in the form (sig1, sig2,..., sigm )
  #phi is the dispersion parameter
  if(phi <= 0 | phi > 1){
    warning('phi should be in (0, 1]')
  }
  nitem = length(sigma)
  #initialize
  ranking = c(sigma[1])
  
  #repeated insersion
  for(i in 2:nitem){
    j = 1:i #the place to insert
    if(phi == 1){
      pvec = rep(1/i, i) #insersion probabilities
    } else{
      pvec = phi^(i-j) * (1-phi) / (1-phi^i) #insersion probabilities
    }
    ind = which(rmultinom(1, 1, pvec) == 1) - 1 #insert after which index
    ranking = append(ranking, sigma[i], ind)
  }
  
  #calculate ranks of each item
  ranks = match(1:nitem, ranking)
  
  return(list(ranks = ranks, ranking = ranking))
  
}



library(FAdist) #for sampling from Gumbel distribution
#this sampler is based on the equivalence
#of Thurstonian model and Luce's Choice
#when utilities are drawn independently
rPL = function(S, beta = 2){
  #S: score vector 
  #beta: arbitrary positive scale parameter
  nitem = length(S)
  
  #compute the location parameter of Gumbel distribution
  mu = beta * S
  
  #sample from Gumbell distribution
  utility = rgumbel(nitem, scale = beta, location = mu)
  names(utility) = 1:nitem #assign labels to items
  #get ranking
  ranking = as.numeric(names(sort(utility, decreasing = TRUE)))
  
  #calculate ranks of each item
  ranks = match(1:nitem, ranking)
  
  return(list(ranks = ranks, ranking = ranking))
}


rThurstone = function(S, Svar){
  #S: Score vector
  #Svar: variance vector
  nitem = length(S)
  utility = rnorm(nitem, S, sqrt(Svar))
  names(utility) = 1:nitem #assign labels to varieties
  #ranking is the items listed in the order A succeeds B succeeds C, etc.
  ranking = as.numeric(names(sort(utility, decreasing = TRUE)))
  
  #calculate ranks of each item
  ranks = match(1:nitem, ranking)
  
  return(list(ranks = ranks, ranking = ranking))
}


rBT = function(S){
  #S: score vector
  nitem = length(S)
  #generate count matrix C
  C = matrix(0, nitem, nitem)
  
  for(i in 1:(nitem - 1)){
    for(j in (i + 1):nitem){
      #generate indicator of whether v_i defeats v_j
      count = rbinom(1, 1, 1 / (1 + exp(-S[i] + S[j])))
      if(count == 1){
        C[i, j] = 1
      } else{
        C[j, i] = 1
      }
    }
  }
  
  #number of wins
  wins = apply(C, 1, sum)
  
  #assign labels to farmers
  names(wins) = 1:nitem
  
  ranking = as.numeric(names(sort(wins, decreasing = TRUE)))
  
  #calculate ranks of each item
  ranks = match(1:nitem, ranking)
  
  return(list(ranks = ranks, ranking = ranking))
}




