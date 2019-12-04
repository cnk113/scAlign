
PCPCA <- function(mat.list, pcs){
  ns <- nrow(mat.list[[1]])
  #  correlation analysis 
  mat_cor <- map(mat.list, ~cor(.))
  avg_mat_cor <- avg_mat(mat_cor)
  pc_mat_cor <- svd(avg_mat_cor)$u

  # find ordering of CPC 
  reconstruct_cor <- map(mat_cor, ~t(pc_mat_cor) %*% . %*% pc_mat_cor)
  cpc_error <- map_dbl(1:pcs, function(i){
    map_dbl(reconstruct_cor, ~pc_cor(., index = i)) %>% mean
  })
  eig_avg__cor <- eig(avg_mat_cor)
  cpc_error <- map_dbl(1:pcs, function(i){
    mean(map_dbl(reconstruct_cor, ~sum(.[i,-i]^2/(eig_avg_cor[i]*eig_avg_cor[-i])))) * ns / (pcs-1)
  })
  reconstruct_eigenvalues <- map(reconstruct_cor, ~diag(.)[order(cpc_error)]) %>% 
    unlist %>% matrix(pcs)
  pc_mat_cor <- pc_mat_cor[,order(cpc_error)]
  return(pc_mat_cor)
}


avg_mat <- function(list_of_mat){
  sum_of_mat <- 0
  for(i in 1:length(list_of_mat)){
    sum_of_mat <- sum_of_mat + list_of_mat[[i]]
  }
  sum_of_mat/length(list_of_mat)
}


pc_cor <- function(matrix, index){
  matrix <- cov2cor(matrix)
  if(length(index) == 1){
    sum((matrix[index,-index])^2)
  }else{
    diag(matrix) <- 0
    sum(matrix[index,]^2)
  }
}
