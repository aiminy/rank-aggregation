library(lme4)
library(EMMREML)


###this function can only be used to analysi tricot comparisons

#require: each variety at least be compared once, i.e. nobs * 3 >= nvar
rankLM = function(data, K = NULL){
  #browser()
  #let m be the number of varieties,
  #let n be the number of farmers.
  #data is an n*m matrix,
  #data(i, j) represents the rank of variety i by farmer j
  #the entry where varieties are not included is 0
  
  #if K is provided, then the analysis is done using additive relationship matrix: K
  #otherwise, then the marker information isn't taken into analysis
  
  nvar = ncol(data)
  nobs = nrow(data)
  
  data_raw = as.vector(data)
  
  data_linear = matrix(0, nvar * nobs, 3)
  # the second column is the label of each variety
  data_linear[, 2] = rep(1:nvar, each = nobs)
  # the third column is the label for each farmer
  data_linear[, 3] = rep(1:nobs, times = nvar)
  # the first column is the rank
  for(i in 1:(nvar * nobs)){
    
    data_linear[i, 1] = data_raw[i]
    
  }
  #delete all the entry where the rank is 0 (not ranked)
  data_linear = data_linear[data_linear[ ,1] != 0, ]
  
  colnames(data_linear) = c('rank', 'variety', 'farmer')
  data_linear = as.data.frame(data_linear)
  data_linear$variety = factor(data_linear$variety)
  data_linear$farmer = factor(data_linear$farmer)
  
  
  if(is.null(K)){
    fit = lmer(rank ~ 1 + (1 | variety), data = data_linear)
    blup_var = coef(fit)$variety #larger value means less competitive
    
    ranking = order(blup_var)
    ranks = match(1:nvar, ranking)
    
  } else{
    
    
    y = data_linear$rank #the response
    X = as.matrix(rep(1, length(y))) #the intercept
    #Z: design matrix for random effects
    Z = matrix(0, length(y), nvar)
    for(i in 1:length(y)){
      
      Z[i, as.numeric(data_linear$variety[i])] = 1
      
    }
    
    fit_m = emmreml(y, X, Z, K)
    
    ranking = order(as.vector(fit_m$uhat))
    ranks = match(1:nvar, ranking)
    
    
  }
  
  return(list(ranks = ranks, ranking = ranking))
  
}


