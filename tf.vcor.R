#' @export tf.vcor
tf.vcor = function(object, print = TRUE)
{
  if(is.tf.model(object)){
    cov = object$corr.coef
    if (print){
      print(cov)
      return(invisible(cov))
    }else{
      return(invisible(cov))
    }
  }else{
    if(is.Arima(object)){
      cov = stats::vcov(object)
      n = ncol(cov)

      var = c(rep(0,n))
      for (i in 1:n) {
        var[i] = cov[i, i]
      }

      for (i in 1:n) {
        for (j in 1:n) {
          cov[i,j] = cov[i, j]/(sqrt(var[i])*sqrt(var[j]))
        }
      }
      if (print){
        print(cov)
        return(invisible(cov))
      }else{
        return(invisible(cov))
      }
    }
  }

}

