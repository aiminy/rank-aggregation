thurstoneSampler = function(mu, theta = c(0, 1)){
  #mu: vector of expected utility (like phenotypic traits)
  #theta = (mean, sd): the parameters passed to normal distributions
  nitem = length(mu)
  utility = mu + rnorm(nitem, theta[1], theta[2]) #add normal noise
  names(utility) = 1:nitem #assign labels to items
  #ranking is the items listed in the order A succeeds B succeeds C, etc.
  ranking = as.numeric(names(sort(utility, decreasing = TRUE)))
  return(ranking)
}
