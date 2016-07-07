#test the best learning rate

#the asymptotic SGD convergence rates are independent of the training data size
#so we decide the best learning rate by using the data
#with 33 farmers and 33 varieties

#number of iterations required and the value of the target function achived are computed

#learning rates seq(0.05, 0.5, 0.05) are tested


#main test function
testRate = function(rate_vec, nsim){
  #browser()
  source('../sampler/rThurstone.R')
  library(agricolae)
  nrate = length(rate_vec)
  #set up parameters
  score = runif(33, 1, 100)
  names(score) = 1:33 #assign labels
  #33 varieties, 33 farmers
  data = matrix(0, 33, 33)
  
  #store the results
  niter = list(thurs = rep(0, nrate), bt = rep(0, nrate)
               , pl = rep(0, nrate), mpm = rep(0, nrate))
  
  target = list(thurs = rep(0, nrate), bt = rep(0, nrate)
                , pl = rep(0, nrate), mpm = rep(0, nrate))
  
  for(j in 1:nsim){
    
    #generate data
    for(i in 1:33){
      data[i, ] = rThurstone(S = score, Svar = 1)$ranks
    }
    #use alpha design to assign varieties to farmers
    design = design.alpha(1:33, 3, 3)$sketch
    #allocation of varieties
    alloc = matrix(0, 33, 3)
    for(i in 1:3){
      alloc[((i - 1) * 11 + 1):(i * 11), ] = matrix(as.numeric(design[[i]]), 11, 3)
    }
    #select the varieties in data
    for(i in 1:33){
      data[i, ][-alloc[i, ]] = 0
    }
    
    #change the rank so that the order is reserved but the ranks only contain 1,2,3
    colnames(data) = 1:33
    for(i in 1:33){
      ranks_temp = data[i, ][data[i, ] != 0]
      sorted_ranks = sort(ranks_temp)
      label = as.numeric(names(sorted_ranks))
      data[i, ][label] = c(1,2,3)
    }
    
    mu = rep(1, 33)
    sigma = diag(1, 33)
    maxiter = 1000
    tol = 1e-6
    start = rep(1, 66)
    decay = 1.05
    
    #fit the data with different learning rates
    for(i in 1:nrate){
      
      temp_res1 = sgdThurs(data, mu, sigma, rate_vec[i], maxiter, tol, start, decay)
      temp_res2 = sgdBT(data, mu, sigma, rate_vec[i], maxiter, tol, start, decay)
      temp_res3 = sgdPL(data, mu, sigma, rate_vec[i], maxiter, tol, start, decay)
      #temp_res4 = sgdMPM(data, mu, sigma, rate_vec[i], maxiter, tol, c(start, rep(1, 33)), decay)
      
      niter[[1]][i] = niter[[1]][i] + temp_res1[[2]]
      niter[[2]][i] = niter[[2]][i] + temp_res2[[2]]
      niter[[3]][i] = niter[[3]][i] + temp_res3[[2]]
      #niter[[4]][i] = niter[[4]][i] + temp_res4[[2]]
      
      target[[1]][i] = target[[1]][i] + temp_res1[[1]][length(temp_res1[[1]])]
      target[[2]][i] = target[[2]][i] + temp_res2[[1]][length(temp_res2[[1]])]
      target[[3]][i] = target[[3]][i] + temp_res3[[1]][length(temp_res3[[1]])]
      #target[[4]][i] = target[[4]][i] + temp_res4[[1]][length(temp_res4[[1]])]
      
      
    }
  }
  
  #average
  for(i in 1:4){
    niter[[i]] = niter[[i]] / nsim
    target[[i]] = target[[i]] / nsim
  }
  
  #plot
  par(mfrow = c(2, 2))
  model_name = c('thurs', 'bt', 'pl', 'mpm')
  #plot for number of iterations
  for(i in 1:4){
    plot(x = rate_vec, y = niter[[i]], type = 'b', xlab = 'rate', ylab = 'niter', 
         main = model_name[i])
  }
  
  #plot for target values
  for(i in 1:4){
    plot(x = rate_vec, y = target[[i]], type = 'b', xlab = 'rate', ylab = 'target', 
         main = model_name[i])
    
  }
  
  return(list(niter = niter, target = target))
  
}

#run the function
res = testRate(seq(0.1, 0.5, 0.05), nsim = 50)
res$rate =  seq(0.5, 1, 0.05)
save(res, file = 'rate_choice.RData')

