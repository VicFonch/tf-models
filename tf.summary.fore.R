tf.summary.fore = function(...){
  funs = list(...)
  nfun = length(funs)
  model.names = paste0("model.", 1:nfun)
  models = list()
  for (i in 1:nfun) {

    res = funs[[i]]

    MSE = mean(res^2)
    RMSE = sqrt(MSE)
    MAE = mean(abs(res))

    models[[i]] = list(MSE, RMSE, MAE)

    salida = c(format(MSE, scientific = F, digits = 4),
               format(RMSE, scientific = F, digits = 4),
               format(MAE, scientific = F, digits = 4))
    names(salida) = c("MSE", "RMSE", "MAE")

    cat(model.names[i],":", "\n")
    print(salida)
    cat("\n")

  }
  return(invisible(models))
}

