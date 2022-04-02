tf.ltf = function (y, x, order = c(k = 20, b = 0), res.order = c(p = 0, d = 0, q = 0),
                          res.Sorder = c(P = 0, D = 0, Q = 0, s = NA), summary = TRUE,
                          include.mean = TRUE, fri.plot = TRUE, fre.plot = FALSE, x.plot)
{
  k = order[1]
  b = order[2]
  p = res.order[1]
  d = res.order[2]
  q = res.order[3]
  P = res.Sorder[1]
  D = res.Sorder[2]
  Q = res.Sorder[3]
  s = res.Sorder[4]
  if (is.na(s)) s = frequency(y)
  ncol = ncol(as.matrix(x))

  if(!is.ts(y)) stop("y debe ser una serie de tiempo")
  if (k <= 0) stop("El orden k  debe ser mayor que 0")
  if (b < 0) stop("El orden b debe ser mayor que 0 ")
  if (d < 0) stop("Los ordenes de differenciacion deben ser mayores o iguales 0.")
  if (D < 0) stop("Los ordenes de differenciacion deben ser mayores o iguales 0.")
  if (p < 0) stop("Orden mal especificado.")
  if (q < 0) stop("Orden mal especificado.")
  if (P < 0) stop("Orden estacional mal especificado.")
  if (Q < 0) stop("Orden estacional mal especificado.")

  if(!is.ts(x)) x = ts(x)
  if(d > 0 || D > 0){
    if(d > 0 && D > 0){
      x = diff(x, differences = d )
      x = diff(x, differences = D, lag = s)
      y = diff(y, differences = d)
      y = diff(y, differences = D, lag = s)
    }else{
      if(d > 0){
        x = diff(x, differences = d)
        y = diff(y, differences = d)
      }else{
        x = diff(x, differences = D, lag = s)
        y = diff(y, differences = D, lag = s)
      }
    }
  }

  if(ncol > 1){
    if (length(y) != nrow(x)){
      dif = abs(length(y) - nrow(x))
      if (length(y) < nrow(x)){
        x = x[1:(nrow(x) - dif), ]
      }else{
        y = y[1:(length(y) - dif)]
      }
    }
  }else{
    if (length(y) != length(x)){
      dif = abs(length(y) - length(x))
      if (length(y) < length(x)){
        x = x[1:(length(x) - dif)]
      }else{
        y = y[1:(length(y) - dif)]
      }
    }
  }

  n = length(y)

  y = y[(b + k):n]
  if (ncol > 1){
    X = x[k:(n - b), 1]
    for (i in 1:ncol) {
      if(k > 1){
        for (j in 1:(k - 1)) {
          X = cbind(X, x[(k - j):(n - b - j), i])
        }
      }
      if (i == ncol) break
      X = cbind(X, x[k:(n - b), (i + 1)])
    }
  }else{
    X = x[k:(n - b)]
    if(k > 1){
      for (i in 1:(k - 1)) {
        X = cbind(X, x[(k - i):(n - b - i)])
      }
    }
  }

  if (include.mean){
    model = forecast::Arima(y, order = c(p, 0, q), seasonal = list(order = c(P, 0, Q), period = s),
                            xreg = X, include.mean = TRUE, method = "CSS")
  }else{
    model = forecast::Arima(y, order = c(p, 0, q), seasonal = list(order = c(P, 0, Q), period = s),
                            xreg = X, include.mean = FALSE, method = "CSS")
  }



  coef = model$coef
  cov.coef = model$var.coef
  n.coef = length(coef)
  se = sqrt(diag(cov.coef))
  if (include.mean){
    n.arima = p + q + P + Q + 1
  }else{
    n.arima = p + q + P + Q
  }
  n.v = k

  if (n.arima > 0) {
    coef.arima = coef[1:n.arima]
    cov.coef = cov.coef[n.arima:n.v, n.arima:n.v]
    se.arima = se[1:n.arima]
    t.a.stadistic = coef.arima/se.arima
    p.a.value = c(rep(0, n.arima))
    for (i in 1:n.arima) {
      p.a.value[i] = pt(q = t.a.stadistic[i], df = n - n.coef)/2
    }
    c.arima = cbind(coef.arima, se.arima, t.a.stadistic, p.a.value)
  }

  if (ncol > 1){
    V = coef[(n.arima + 1):n.coef]
    SE.V = se[(n.arima + 1):n.coef]
    v = V[1:k]
    se.v = SE.V[1:k]
    for (i in 1:(ncol - 1)) {
      v = cbind(v, V[(i*k + 1):((i + 1)*k)])
      se.v = cbind(se.v, SE.V[(i*k + 1):((i + 1)*k)])
    }
    t.v.stadistic = v/se.v
    p.v.value = matrix(0, n.v, ncol)
    for (i in 1:(n.v)) {
      for (j in 1:ncol) {
        p.v.value[i,j] = pt(q = abs(t.v.stadistic[i,j]), df = n - n.coef, lower.tail = FALSE)/2
      }
    }

    c.v = list()
    for (i in 1:ncol) {
      c.v[[i]] = cbind(v[, i], se.v[, i], t.v.stadistic[, i], p.v.value[, i])
      colnames(c.v[[i]]) = c("v", "se.v", "t.v.stadistic", "p.v.value")
      rownames(c.v[[i]]) = paste0("v.", 1:k)
    }

  }else{
    v = coef[(n.arima + 1):n.coef]
    names(v) = paste0("v.", 1:n.v)
    if(k > 1){
      colnames(X) = paste0("v.", 1:n.v)
    }else names(X) = "v.1"

    se.v = se[(n.arima + 1):n.coef]
    t.v.stadistic = v/se.v
    p.v.value = c(rep(0, n.v))
    for (i in 1:(n.v)) {
      p.v.value[i] = pt(q = abs(t.v.stadistic[i]), df = n - n.coef, lower.tail = FALSE)/2
    }
    c.v = cbind(v, se.v, t.v.stadistic, p.v.value)
  }


  if (fri.plot){
    if(ncol > 1){
      if (missing(x.plot)){
        plot(ggplot(data = as.data.frame(v[ ,1]), aes(y = v[ ,1], x = 1:n.v)) + geom_bar(stat = "identity")+
               geom_hline(yintercept = 0)+
               labs(x = "coefficients",y = "FRI"))
      }else{
        plot(ggplot(data = as.data.frame(v[ ,x.plot]), aes(y = v[ ,x.plot], x = 1:n.v)) + geom_bar(stat = "identity")+
               geom_hline(yintercept = 0)+
               labs(x = "coefficients",y = "FRI"))
      }
    }else{
      plot(ggplot(data = as.data.frame(v), aes(y = v, x = 1:n.v)) + geom_bar(stat = "identity")+
             geom_hline(yintercept = 0)+
             labs(x = "coefficients",y = "FRI"))
    }

  }

  if(fre.plot){
    if (ncol > 1){
      if (missing(x.plot)){
        fre = c(rep(0, n.v))
        for (i in 1:n.v) {
          for (j in 1:i) {
            fre[i] = fre[i] + v[j, 1]
          }
        }
      }else{
        fre = c(rep(0, n.v))
        for (i in 1:n.v) {
          for (j in 1:i) {
            fre[i] = fre[i] + v[j, x.plot]
          }
        }
      }
    }else{
      fre = c(rep(0, n.v))
      for (i in 1:n.v) {
        for (j in 1:i) {
          fre[i] = fre[i] + v[j]
        }
      }
    }

    plot(ggplot(data = as.data.frame(fre), aes(y = fre, x = 1:n.v)) + geom_bar(stat = "identity")+
           geom_hline(yintercept = 0)+
           labs(x = "coefficients",y = "FRE"))
  }

  if (ncol > 1){
    if(summary){
      res = paste0("Funcion Respuesta a Impulsos & s.e & t.test de la variable X." , 1:ncol)
      for (i in 1:ncol) {
        cat(res[i], "\n")
        print(t(c.v[[i]]), digits = 3)
        cat("\n")
      }
      cat("\n")
      if (n.arima > 0){
        cat("Coeficientes ARIMA & s.e & t.test.:", "\n")
        print(t(c.arima), digits = 3)
      }
    }
  }else{
    if (summary){
      cat("Funcion Respuesta a Impulsos & s.e & t.test.:", "\n")
      print(t(c.v), digits = 3)
      cat("\n")
      if (n.arima > 0){
        cat("Coeficientes ARIMA & s.e & t.test.:", "\n")
        print(t(c.arima), digits = 3)
      }
    }
  }

  if(k > 1){
    cor = cor(X)
  }else cor = NULL

  if(k > 1){
    fitted = ts(c(y[1:(k - 1)], model$fitted))
    residuals = ts(c(model$residuals[1:(k - 1)], model$residuals))
  }else{
    fitted = model$fitted
    residuals = model$residuals
  }

  if (n.arima > 0){
    t.test = list(t.v.stadistic = t.v.stadistic, p.v.value = p.v.value,
                  t.a.stadistic = t.a.stadistic, p.a.value = p.a.value)
    return(invisible(list(v = v, se.v = se.v, cov.v = cov.coef, coef.arima = coef.arima,
                          se.arima = se.arima, t.test = t.test,fitted = fitted,
                          residuals = residuals, cor = cor)))
  }else{
    t.test = list(t.v.stadistic = t.v.stadistic, p.v.value = p.v.value)
    return(invisible(list(v = v, se.v = se.v, cov.v = cov.coef, t.test = t.test,
                          fitted = fitted, residuals = residuals, cor = cor)))
  }

}

