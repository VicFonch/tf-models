tf.tsplot = function(...)
{
  y = list(...)
  l = length(y)
  l1 = length(y[[1]])
  if(l > 1){
    for (i in 2:l) {
      if(length(y[[i]]) != l1) stop("Argumantos difieren en tsmaño")
    }
    x = ts(data.frame(...))
  }else{
    x = ts(...)
  }
  plot(ggplot2::autoplot(x))
}


