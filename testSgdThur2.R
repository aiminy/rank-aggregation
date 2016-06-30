#only used for small dataset
#should be optimized for large dataset

sgdThur2 = function(data, mu, sigma, rate, maxiter = 1000, tol = 1e-9, start){
  #let m be the number of varieties,
  #let n be the number of farmers.
  #data is an n*m matrix,
  #data(i, j) represents the rank of variety i by farmer j
  #the entry where varieties are not included is 0
  
  #rate is the step size/learning rate
  source('targetThur2.R')
  
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
      res_temp = targetThur2(score_temp, adherence_temp, data, mu, sigma)
      
      #store the value of the target function
      target[niter] = res_temp[[1]]
      
      #extract the gradient
      gradient = res_temp[[2]][i, ]
      
      #update the parameters
      param = param - rate * gradient
      
      #check the convergence criteria: square of the change of target values
      if(niter > 1){
        if((target[niter] - target[niter - 1]) ^ 2 < tol | niter > maxiter){
          flag = FALSE
          break
        }
      }
      
      
    }
    
  }
  
  return(list(value = target, niter = niter, score = param[1:nvar],
              adherence = param[(nvar + 1):(nvar + nobs)]))
  
}


#test
set.seed(1)
score = runif(99, 1, 100)
names(score) = 1:99 #assign labels
#99 varieties, 99 farmers
data = matrix(0, 99, 99)
source('../sampler/rThurstone.R')
for(i in 1:99){
  data[i, ] = rThurstone(S = score, Svar = 1)$ranks
}
#View(data)
#use alpha design to assign varieties to farmers
set.seed(10)
library(agricolae)
design = design.alpha(1:99, 3, 3)$sketch
#allocation of varieties
alloc = matrix(0, 99, 3)
for(i in 1:3){
  alloc[((i - 1) * 33 + 1):(i * 33), ] = matrix(as.numeric(design[[i]]), 33, 3)
}
#select the varieties in data
for(i in 1:99){
  data[i, ][-alloc[i, ]] = 0
}
#View(data)

res = sgdThur2(data, sigma = diag(rep(1, 99)), rate = 0.1, 
             mu = rep(0, 99), start = rep(1, 198), maxiter = 50000, tol = 1e-9)

#plot to see if the function value decreases
plot(1:length(res[[1]]), res[[1]], type = 'l')

#estimated scores
est_score = res$score
names(est_score) = 1:99 #assign labels
#calculate ranking
est_ranking = as.numeric(names(sort(score, decreasing = T)))
real_ranking = as.numeric(names(sort(est_score, decreasing = T)))
est_ranking
real_ranking
#calculate ranks
real_rank = match(1:99, real_ranking)
est_rank = match(1:99, est_ranking)
#calculate kendall's correlation coefficient
cor(cbind(est_rank, real_rank), method = 'kendall')
#the adherence
res$adherence

