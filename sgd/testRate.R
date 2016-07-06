#test the best learning rate
set.seed(10)
score = runif(99, 1, 100)
names(score) = 1:99 #assign labels
#99 varieties, 99 farmers
data = matrix(0, 99, 99)
source('../sampler/rThurstone.R')
for(i in 1:99){
  data[i, ] = rThurstone(S = score, Svar = 1)$ranks
}

#use alpha design to assign varieties to farmers
set.seed(100)
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

#change the rank so that the order is reserved but the ranks only contain 1,2,3
colnames(data) = 1:99
for(i in 1:99){
  ranks_temp = data[i, ][data[i, ] != 0]
  sorted_ranks = sort(ranks_temp)
  label = as.numeric(names(sorted_ranks))
  data[i, ][label] = c(1,2,3)
}