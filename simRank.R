source('sampler.R')
source('sgd.R')

#main simulation function
simRank = function(h2, nvar_vec, nobs_vec, marker, marker_mu, marker_Sig, nsim, seed){
  library(MASS)
  library(agricolae) #generating alpha designs
  library(rrBLUP)
  #browser()
  #h2: heritability;
  
  #nvar_vec: vector containing number of varieties;
  
  #nobs_vec: vector containing number of farmers/observers;
  
  #marker: marker matrix coded in {0, 1, 2}, with nrow = number of clones
  #and ncol = number of loci, we also assume nrow > max(nobs_vec)
  
  #marker_mu: the mean vector for causal loci effects
  #marker_Sig: the covariance matrix for causal loci effects
  #length of marker_mu is the size of the causal loci
  
  #nsim: number of iterations in the simulation
  
  #prop: the proportion of loci chosen to be the causal loci
  
  #seed: the seed for randomization
  
  set.seed(seed)
  
  #set up the file to be written
  write.table(t(c('Sampler', 'Model', 'nvar', 'nobs', 'Kendall_marker', 'Kendall_nomarker')), 
            file = 'simRes.csv', col.names = F, row.names = F, sep = ',')
  
  #choose the varieties we want to analyze
  marker = marker[sample(nrow(marker), max(nvar_vec)), ]
  
  #imputation
  marker = A.mat(marker - 1, return.imputed = T)$imputed
  
  #marker matrix coded in {-1, 0, 1}
  marker2 = marker
  
  #marker matrix coded in {0, 1, 2}
  marker = marker + 1
  
  #do the simulation multiple times
  for(niter in 1:nsim){
    cat('Now at iteration ', niter, '\n',sep = '')
    
    #sample causal/non_causal loci
    ncausal = length(marker_mu)
    causal_label = sample(ncol(marker), ncausal)
    non_causal_label = as.vector(1:ncol(marker))[-causal_label]
    
    #causal marker matrix and non_causal marker matrix
    causal = marker[ ,causal_label] #coded in {0, 1, 2}
    non_causal = marker2[, non_causal_label] #coded in {-1, 0, 1}
    
    #sample causal loci marker effects
    marker_eff = mvrnorm(1, marker_mu, marker_Sig)
    
    #score vector
    score = as.vector(causal %*% marker_eff)
    var_score = var(score)
    
    #compute Gaussian error from observers, given heritability
    var_obs = (1 / h2 - 1) * var_score
    
    
    #additive relationship matrix
    A = A.mat(non_causal, shrink = T) #shrink in order to let the condition number to be large
    
    #generate synthetic ranking data from ranking models
    R_full = matrix(0, max(nobs_vec), max(nvar_vec)) #full ranking data
    for(sampler in c('Thurs', 'PL')){
      
      
      if(sampler == 'Thurs'){
        
        samplerFunc = rThurstone
        #parameter related to Gaussian errors passed to the sampler
        var_param = rep(var_obs, nvar)
      }
      
      if(sampler == 'PL'){
        
        samplerFunc = rPL
        #parameter related to Gaussian errors passed to the sampler
        var_param = sqrt(6 * var_obs / pi^2)
        
        #tune the score so that the location parameter in the gumbel distribution
        #is the original score
        score = score / var_param
      }
      
      #generate full ranking data
      for(i in 1:max(nobs_vec)){
        
        R_full[i, ] = samplerFunc(score, var_param)$ranks
        
      }
      
      
      #loop over all possible number of varieties
      
      for(i in 1:length(nvar_vec)){
        
        nvar = nvar_vec[i]
        
        
        #loop over all possible number of observers
        
        for(j in 1:length(nobs_vec)){
          
          nobs = nobs_vec[j]
          
          varieties = sample(max(nvar_vec), nvar)
          observers = sample(max(nobs_vec), nobs)
          
          #score in this (nvar, nobs) combination
          score_each = score[varieties]
          names(score_each) = 1:nvar #assign labels to varieties, from 1 to nvar
          
          #covariance matrix in the normal prior on scores
          prior_cov = A[varieties, varieties]
          
          #full ranking with nobs observers and nvar varieties
          R_each = R_full[observers, varieties]
          colnames(R_each) = 1:nvar #assign labels to varieties, from 1 to nvar 
          
          #trim the full ranking into partial ranking according to an alpha design
          ratio = ceiling(nobs / nvar)
          #generate alpha designs that can accomodate nobs observers
          design = matrix(0, nvar * ratio, 3) #each row is the allocation of varieties
          for(k in 1:ratio){
            
            design_temp = design.alpha(1:nvar, 3, 3)$sketch
            design_temp = rbind(design_temp[[1]], design_temp[[2]], design_temp[[3]])
            
            design[((k - 1) * nvar + 1):(k * nvar) , ] = design_temp

          }
          design = design[sample(nvar * ratio, nobs), ]
          design = matrix(as.numeric(design), nobs, 3)
          
          #select varieties according to the alpha design in R_each
          for(k in 1:nobs){
            
            R_each[k, ][-design[k, ]] = 0 #varieties not ranked are set to 0
            ranks_temp = R_each[k, ][R_each[k, ] != 0] #ranks of varieties included
            sorted_ranks = sort(ranks_temp)
            label = as.numeric(names(sorted_ranks)) #ranking in the form label[1] \succ label[2] \succ ...
            R_each[k, ][label] = c(1,2,3) #assign ranks 1,2,3
            
          }
          
          
          #fit the data using different ranking models
          for(model in c('Thurs', 'BT', 'PL', 'MPM')){
            
            #set up the parameters in sgd functions
            prior_mu = rep(1, nvar)
            identity_cov = diag(1, nvar)
            decay = 1.1
            start = rep(1, nvar + nobs) #starting point
            
            
            if(model == 'Thurs'){
              
              sgdFunc = sgdThurs
              rate = 0.25
              
            }
            
            if(model == 'BT'){
              
              sgdFunc = sgdBT
              rate = 0.48
              
            }
            
            if(model == 'PL'){
              
              sgdFunc = sgdPL
              rate = 0.48
              
            }
            
            if(model == 'MPM'){
              
              sgdFunc = sgdMPM
              start = c(start, rep(1, nvar))
              rate = 0.1
            }
            
            #result incorporating relationship matrix
            res1 = sgdFunc(R_each, prior_mu, prior_cov, rate, maxiter = 1000, tol = 1e-8, start, decay)
            names(res1$score) = 1:nvar
            
            #result with identity cov.
            res2 = sgdFunc(R_each, prior_mu, identity_cov, rate, maxiter = 1000, tol = 1e-8, start, decay)
            names(res2$score) = 1:nvar
            
            #compute real ranking according to score_each
            real_ranking = as.numeric(names(sort(score_each, decreasing = T)))
            #estimated ranking
            est_ranking1 = as.numeric(names(sort(res1$score, decreasing = T)))
            est_ranking2 = as.numeric(names(sort(res2$score, decreasing = T)))
            
            #compute real ranks
            real_rank = match(1:nvar, real_ranking)
            
            #estimated ranks
            est_rank1 = match(1:nvar, est_ranking1)
            est_rank2 = match(1:nvar, est_ranking2)
            
            #calculate kendall's correlation coefficient
            cor1 = cor(cbind(est_rank1, real_rank), method = 'kendall')[1, 2]
            cor2 = cor(cbind(est_rank2, real_rank), method = 'kendall')[1, 2]
            
            #save the Kendall's correlation coefficient
            write.table(t(c(sampler, model, nvar, nobs, cor1, cor2)), 
                        file = 'simRes.csv', col.names = F, row.names = F,
                        append = TRUE, sep = ',')
            
            #check the performance of sgd every 100 iterations
            if(niter %% 100 == 1){
              
              #plot
              #file name is ordered as sampler_model_nvar_nobs
              pdf(file = paste(paste(sampler, model, nvar, nobs, sep = '_'), '.pdf', sep = ''),
                  width = 8.5, height = 11)
              par(mfrow = c(2, 1))
              plot(x = 1:length(res1[[1]]), y = res1[[1]], type = 'l', 
                   xlab = 'niter', ylab = 'target value', 
                   main = 'With marker info')
              plot(x = 1:length(res2[[1]]), y = res2[[1]], type = 'l', 
                   xlab = 'niter', ylab = 'target value', 
                   main = 'Without marker info')
              dev.off()
              
            }
            
            
            
          }
     
          
        }
        
        
      }
      
      
      
    }
    
    
    
    
    
    
    
    
  }
  
  
  
  
  
}



simRank(0.8, 9, 10, data, runif(50, 1, 30), diag(1, 50), 2, 266)
