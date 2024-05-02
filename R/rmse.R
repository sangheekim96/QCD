# compute rmse with fitted beta

rmse = function(beta, truebeta, lambda) {
  rmse = matrix(NA, nrow=length(lambda), ncol=1)
  for(i in 1:length(lambda)){
    rmse[i,] = sqrt(norm(beta[i,] - matrix(truebeta, nrow=1), "2")/norm(matrix(truebeta, nrow=1), "2"))
  }
  return(rmse)
}
