library(FAdist) #for sampling from Gumbel distribution
#this sampler is based on the equivalence
#of Thurstonian model and Luce's Choice
#when utilities are drawn independently
plSampler = function(ability, beta = 2){
  #ability is the vector of parameters in Plackett-Luce model, 
  #and is normalized if not sum to 1;
  #beta is the scale parameter of Gumbel distribution.
  nitem = length(ability)
  #normalize
  ability = ability / sum(ability)
  #compute the location parameter of Gumbel distribution
  mu = beta * log(ability)
  #sample from Gumbell distribution
  utility = rgumbel(nitem, scale = beta, location = mu)
  names(utility) = 1:nitem #assign labels to items
  #get rank
  ranking = as.numeric(names(sort(utility, decreasing = TRUE)))
  
  return(ranking)
}

