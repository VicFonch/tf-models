tf.summary = function(...)
{

  funs = list(...)
  nfun = length(funs)
  model.names = paste0("model.", 1:nfun)
  models = list()
  for (i in 1:nfun) {
    if(is.tf.model(funs[[i]]) || is.Arima(funs[[i]])){
      if(is.tf.model(funs[[i]])){
        y = funs[[i]]$model$y
        y.ad = funs[[i]]$fitted
        n = length(funs[[i]]$model$y)
        k = length(funs[[i]]$coefficient$coefficient)
        gl = n - k
        lik = funs[[i]]$log.likelihood
        AIC = funs[[i]]$AIC
        AICc = funs[[i]]$AICc
        BIC = funs[[i]]$BIC
      }else{
        y = funs[[i]]$x
        y.ad = funs[[i]]$fitted
        n = length(funs[[i]]$x)
        k = length(funs[[i]]$coef)
        gl = n - k
        lik = funs[[i]]$loglik
        AIC = funs[[i]]$aic
        AICc = funs[[i]]$aicc
        BIC = funs[[i]]$bic
      }
      res = funs[[i]]$residuals

      DW = stats::acf(res, plot = FALSE, lag.max = 1)$acf[2]
      DW = 2*(1 - DW)
      R2 = sum((y.ad - mean(y))^2)/sum((y - mean(y))^2)
      R2.a = 1 - (1 - R2)*(n - 1)/(n - k - 1)
      MSE = mean(res^2)
      RMSE = sqrt(MSE)
      MAE = mean(abs(res))
      PMAE = mean(abs(res/as.vector(y)))*100

      models[[i]] = list(n = n, gl = gl, DW = DW, R2 = R2, R2.a = R2.a, Loglik = lik, AIC = AIC,
                         AICc = AICc, BIC = BIC, MSE = MSE, RMSE = RMSE, MAE = MAE, PMAE = PMAE)
      salida1 = c(n, gl)
      salida2 = c(format(DW, scientific = F, digits = 4), format(R2, scientific = F, digits = 4),
                  format(R2.a, scientific = F, digits = 4), format(lik, scientific = F, digits = 4),
                  format(AIC, scientific = F, digits = 4), format(AICc, scientific = F, digits = 4),
                  format(BIC, scientific = F, digits = 4),format(MSE, scientific = F, digits = 4),
                  format(RMSE, scientific = F, digits = 4),format(MAE, scientific = F, digits = 4),
                  format(PMAE, scientific = F, digits = 4))
      names(salida1) = c("n", "gl")
      names(salida2) = c("DW", "R2", "R2.a", "LogLik", "AIC", "AICc", "BIC",
                         "MSE", "RMSE", "MAE", "PMAE")

      cat(model.names[i],":", "\n")
      print(salida1)
      print(salida2)
      cat("\n")

    }else stop("Los argumentos solo pueden ser de tipo tf.model o  Arima(forecast)")



  }
  return(invisible(models))
}



