#only used for small dataset
#should be optimized for large dataset
sgdBT = function(data, mu, sigma, rate, maxiter = 1000, tol = 1e-9, start){
  source('targetBT.R')
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
  score = start
  target = rep(0, niter)
  
  flag = TRUE
  #loop until meet the convergence criteria
  while(flag){
    
    #shuffle the data
    #perm = sample(nobs)
    for(i in 1:nobs){ 
      
      niter = niter + 1
      #evaluate the log-posterior as well as the gradient
      #only used for small dataset (where we want to decide the learning rate)
      #if used for big dataset, where we don't want to
      #evaluate likelihood everytime, the function should be modified
      res_temp = targetBT(score, data, mu, sigma)
      #store the value of the target function
      target[niter] = res_temp[[1]]
      
      #extract the gradient
      gradient = res_temp[[2]][i, ]
      
      #update score
      score = score - rate * gradient
      
      #check the convergence criteria: square of the change of target values
      if(niter > 1){
        if((target[niter] - target[niter - 1]) ^ 2 < tol | niter > maxiter){
          flag = FALSE
          break
        }
      }
    
    }
  }
  
  return(list(value = target, niter = niter, score = score))
}



#test
set.seed(1)
score = runif(100, 1, 100)
names(score) = 1:100 #assign labels
data = matrix(0, 100, 100)
for(i in 1:100){
  data[i, ] = rThurstone(S = score, Svar = 1)$ranks
}


res = sgdBT(data, sigma = diag(rep(1, 100)), rate = 0.1, 
             mu = rep(1, 100), start = rep(1, 100), maxiter = 20, tol = 1e-9)

#plot to see if the function value decreases
plot(1:length(res[[1]]), res[[1]], type = 'l')

#estimated scores after 20 iterations
est_score = res$score
names(est_score) = 1:100 #assign labels
#calculate ranking
est_ranking = as.numeric(names(sort(score, decreasing = T)))
real_ranking = as.numeric(names(sort(est_score, decreasing = T)))
est_ranking
real_ranking
#calculate ranks
real_rank = match(1:100, real_ranking)
est_rank = match(1:100, est_ranking)
#20 steps in sgd produces a fairly good estimation of ranking
cor(cbind(est_rank, real_rank), method = 'kendall')

