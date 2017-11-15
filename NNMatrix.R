
#### distance matrix for location i and its neighbors ####
i_dist <- function(i, neighbor_index, s)	dist(s[c(i, neighbor_index[[i - 1]]), ])

get_NN_distM <- function (ind, ind_distM_d) {
  if (ind < M ){l = ind } else {l = M}; 
  M_i <- rep(0, M * (M - 1) / 2);
  if (l == 1) {}
  else{
    M_i[1: (l * (l - 1) / 2)] <- 
      c(ind_distM_d[[ind]])[(l + 1): (l * (l + 1) / 2)]
  }
  return(M_i)
}

get_NN_dist <- function (ind, ind_distM_d) {
  if (ind < M ){l = ind } else {l = M}; 
  D_i <- rep(0, M);
  D_i[1:l] <- c(ind_distM_d[[ind]])[1:l]
  return(D_i)
}

get_NN_ind <- function (ind, ind_distM_i) {
  if (ind < M ){l = ind } else {l = M}; 
  D_i <- rep(0, M);
  D_i[1:l] <- c(ind_distM_i[[ind]])[1:l]
  return(D_i)
}


#### wrap up in function NNMatrix ####

NNMatrix <- function(N, coords.ord, n.indx){
    
    NN_ind <- t(sapply(1: (N - 1), get_NN_ind, n.indx))
    neighbor_dist <- sapply(2:N, i_dist, n.indx, coords.ord)
    NN_distM <- t(sapply(1: (N - 1), get_NN_distM, neighbor_dist))
    NN_dist <- t(sapply(1: (N - 1), get_NN_dist, neighbor_dist))
    
    return(
      list(NN_ind = NN_ind, NN_distM = NN_distM,
           NN_dist = NN_dist))
}

#### Function for checking neighbors ####

Check_Neighbors <- function(coords, NN, NN.matrix, ind){
  
  plot(coords) 
  points(coords[1:ind, , drop = FALSE], col = "grey", pch = 19) 
  
  # neighbors
  if (ind < NN) {dim = ind} else {dim = NN}
  for (j in 1:dim){
    points(coords[NN.matrix$NN_ind[ind - 1, j], , drop = FALSE], 
           col = "orange",  pch = 19)
  }
  points(coords[ind, , drop = FALSE], col = "blue", pch = 19) 
  
}

