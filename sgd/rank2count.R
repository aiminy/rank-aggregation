#transform the ranking data into a count matrix
rank2count = function(data){
  
  #inpute data is a vector,
  #where i-th element in the vector is the rank assigned to the i-th item
  #the entry where item not ranked is replaced by 0
  
  #assume there are no ties
  
  nitem = length(data)
  names(data) = 1:length(data) #assign labels
  
  rank = data[data != 0]
  nitem_ranked = length(rank)
  C = matrix(0, nitem, nitem)
  
  #loop over all pairwise comparisons
  for(i in 1:(nitem_ranked - 1)){
    
    for(j in (i + 1):nitem_ranked){
      
      entry = rank[i] - rank[j]
      item_i = as.numeric(names(rank)[i])
      item_j = as.numeric(names(rank)[j])
      
      if(entry > 0){ #item_j wins over item_i
        C[item_j, item_i] = entry
      } else{  #item_i wins over item_j
        C[item_i, item_j] = -entry
      }
      
    }
  }
  
  return(C)
  
}

