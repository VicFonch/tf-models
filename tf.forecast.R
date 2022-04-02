#' @export tf.forecast
tf.forecast = function(object, h = 12, arima.x = c(p,d,q), Sarima.x = c(P,D,Q,s),
                       x.predict ,plot = TRUE)
{

  fore.1 = function(y, x, res, c, d.coef, d.ret, o.coef, o.ret, q.coef, q.ret,
                    tf.order, arima.order){

    a = tf.order[1]
    b = tf.order[2]
    m = tf.order[3]
    p = arima.order[1]
    q = arima.order[3]
    P = arima.order[4]
    Q = arima.order[6]
    s = arima.order[7]
    n = length(x)

    y = y - c

    y.t = 0
    if(a > 0 || p > 0 || P > 0){
      for (i in d.ret) {
        y.t = y.t + d.coef[i + 1]*y[n - i]
      }
    }

    x.t = 0
    for (i in o.ret) {
      x.t = x.t + o.coef[i + 1]*x[n - b - i]
    }

    a.t = 0
    n.res = length(res)
    if(a > 0 || q > 0 || Q > 0){
      for (i in q.ret) {
        if(n - i < n.res) a.t = a.t + q.coef[i + 1]*res[n - i]
      }
    }

    fore = y.t + x.t + a.t + c

    return(fore)
  }

  fore.2 = function(y, x, res, c, d.coef, d.ret, o.coef, o.ret, q.coef,
                    q.ret, tf.order, arima.order, ncol){

    a = tf.order[[1]]
    b = tf.order[[2]]
    m = tf.order[[3]]
    p = arima.order[1]
    q = arima.order[2]
    P = arima.order[3]
    Q = arima.order[4]
    s = arima.order[5]
    n = nrow(x)

    y = y - c

    n.a = 0
    for (i in 1:length(a)) {
      if(a[i] > 0) n.a = n.a + 1
    }

    if(n.a > 0 || P > 0 || p > 0){
      y.t = 0
      for (i in d.ret) {
        y.t = y.t + d.coef[i + 1]*y[n - i]
      }
    }

    x.t = 0
    for (i in 1:ncol) {
      for (j in o.ret[[i]]) {
        x.t = x.t + o.coef[[i]][j + 1]*x[(n - b[i] - j), i]
      }
    }

    a.t = 0
    n.res = length(res)
    if(n.a > 0 || Q > 0 || q > 0){
      for (i in q.ret) {
        if(n - i < n.res) a.t = a.t + q.coef[i + 1]*res[n - i]
      }
    }

    fore = y.t + x.t + a.t + c

    return(fore)

  }

  #FUNCIONEEEEEEEE

  if(is.tf.model(object)){

    coefficients = object$coefficients$coefficients
    y = object$model$y
    y.real = y
    x = object$model$x
    xreg = object$model$xreg
    det.seasonal = object$model$y.seasonal
    det.trend = object$model$y.trend
    type.deter = object$model$y.type.deter

    p = object$model$arima[1]
    d = object$model$arima[2]
    q = object$model$arima[3]
    P = object$model$arima[4]
    D = object$model$arima[5]
    Q = object$model$arima[6]
    s = object$model$arima[7]
    arima.order = object$model$arima

    n = length(y)
    ncol = ncol(as.matrix(x))

    if(det.seasonal || det.trend > 0){
      if(det.seasonal && det.trend > 0){

        det.s = tf.decompose(y, type = type.deter)$seasonal
        for (i in 1:s) {
          if(det.s[n] == det.s[i]) {
            seson = i
            break
          }
        }
        X = NULL
        for (i in (0:det.trend)) X = cbind(X, (1:n)^i)
        ind = solve.qr(qr(X), y)
        det.t = as.vector(X%*%ind)

        if(type.deter == "additive"){
          y = y - det.s - det.t
        }else y = y/(det.s*det.t)

      }else{
        if(det.seasonal){
          det.s = tf.decompose(y, type = type.deter)$seasonal
          for (i in 1:s) {
            if(det.s[n] == det.s[i]) {
              seson = i
              break
            }
          }
          if(type.deter == "additive"){
            y = y - det.s
          }else y = y/det.s

        }else{
          X = NULL
          for (i in (0:det.trend)) X = cbind(X, (1:n)^i)
          ind = solve.qr(qr(X), y)
          det.t = as.vector(X%*%ind)
          if(type.deter == "additive"){
            y = y - det.t
          }else y = y/det.t
        }
      }
    }

    if(!is.ts(x)) x = ts(x)

    if(d > 0){
      if(d > 1){
        diff.pol = polynom::polynomial(c(1, -1))
        cont = diff.pol
        for (i in 2:d) {
          diff.pol = diff.pol*cont
        }
      }else{
        diff.pol = polynom::polynomial(c(1, -1))
      }
    }else diff.pol = 1

    if(D > 0){
      if(D > 1){
        Sdiff.pol = polynom::polynomial(c(1, rep(0, s - 1), -1))
        cont = Sdiff.pol
        for (i in 2:D) {
          Sdiff.pol = Sdiff.pol*cont
        }
      }else{
        Sdiff.pol = polynom::polynomial(c(1, rep(0, s - 1), -1))
      }
    }else Sdiff.pol = 1

    diff.pol = diff.pol*Sdiff.pol


    if(ncol > 1){

      px = NULL
      dx = NULL
      qx = NULL
      Px = NULL
      Dx = NULL
      Qx = NULL
      sx = NULL
      if(!missing(arima.x)){
        for (i in 1:ncol) {
          if(!is.list(arima.x)) stop("Elemento arima.x mal especificado")
          if(!is.list(Sarima.x)) stop("Elemento Sarima.x mal especificado")
          if(length(arima.x[[i]]) != 3) stop("Elemento arima.x mal especificado")
          if(length(Sarima.x[[i]]) != 4) stop("Elemento Sarima.x mal especificado")
          px = c(px, arima.x[[i]][1])
          dx = c(dx, arima.x[[i]][2])
          qx = c(qx, arima.x[[i]][3])
          Px = c(Px, Sarima.x[[i]][1])
          Dx = c(Dx, Sarima.x[[i]][2])
          Qx = c(Qx, Sarima.x[[i]][3])
          sx = c(sx, Sarima.x[[i]][4])
        }
        for (i in 1:ncol) {
          if(px[i] < 0) stop("arima.x mal especificado")
          if(qx[i] < 0) stop("arima.x mal especificado")
          if(Px[i] < 0) stop("arima.x mal especificado")
          if(Qx[i] < 0) stop("arima.x mal especificado")
          if (dx[i] < 0) stop("Los ordenes de differenciacion deben ser > 0.")
          if (sx[i] <= 0) stop("Componente (s) mal especificado")
          if (Dx[i] < 0) stop("Los ordenes de differenciacion deben ser > 0.")
        }

      }else{
        for (i in 1:ncol) {
          arim = forecast::auto.arima(x[, i], method = "CSS", d = 2, max.D = 1,
                                      max.p = 3, max.P = 2, max.q = 3, max.Q = 2)
          arim = arim$arma
          px = c(px, arim[1])
          dx = c(dx, arim[6])
          qx = c(qx, arim[2])
          Px = c(Px, arim[3])
          Dx = c(Dx, arim[7])
          Qx = c(Qx, arim[4])
          sx = c(sx, arim[5])
        }
      }

      a = object$model$tf[[1]]
      b = object$model$tf[[2]]
      m = object$model$tf[[3]]
      tf.order = list(a, b ,m)

      if(is.null(object$coefficients$Intercept)){
        coefficients = c(0, coefficients)
        c = 0
        intercep = FALSE
      }else{
        c = coefficients[1]
        intercep = TRUE
      }

      omega = list()
      delta = list()
      o.pol = list()
      d.pol = list()

      omega[[1]] = coefficients[2:(m[1] + 2)]
      if(m[i] > 0){
        o.pol[[1]] = polynomial(c(omega[[1]][1], -coefficients[3:(m[1] + 2)]))
      }else{
        o.pol[[1]] = omega[[1]]
      }
      its = m[1] + 2
      for (i in 2:ncol) {
        omega[[i]] = coefficients[(its + 1):(its + m[i] + 1)]
        if(m[i] > 0){
          o.pol[[i]] = polynomial(omega[[i]])
        }else{
          o.pol[[i]] = omega[[i]]
        }
        its = its + m[i] + 1
      }

      for (i in 1:ncol) {
        if(a[i] > 0){
          delta[[i]] = coefficients[(its + 1):(its + a[i])]
          d.pol[[i]] = polynomial(c(1, -delta[[i]]))
          if (i == 2) d.acum = d.pol[[1]]*d.pol[[2]]
          if (i > 2) d.acum = d.acum*d.pol[[i]]
          its = its + a[i]
        }else{
          d.pol[[i]] = 1
          if (i == 2) d.acum = d.pol[[1]]*d.pol[[2]]
          if (i > 2) d.acum = d.acum*d.pol[[i]]
        }
      }

      if (p > 0){
        phi = coefficients[(its + 1):(its + p)]
        p.pol = polynomial(c(1, -phi))
        its = its + p
      }else p.pol = 1

      if (q > 0){
        theta = coefficients[(its + 1):(its + q)]
        q.pol = polynomial(c(1, -theta))
        its = its + q
      }else q.pol = 1

      if (P > 0){
        Phi = coefficients[(its + 1):(its + P)]
        P.pol = polynomial(c(1, rep(0, s - 1), -Phi[1]))
        if(P > 1){
          for (i in 2:P) {
            P.pol = polynomial(c(P.pol, rep(0, s - 1), -Phi[i]))
          }
        }
        its = its + P
      }else P.pol = 1

      if (Q > 0){
        Theta = coefficients[(its + 1):(its + Q)]
        Q.pol = polynomial(c(1, rep(0, s - 1), -Theta[1]))
        if(Q > 1){
          for (i in 2:Q) {
            Q.pol = polynomial(c(Q.pol, rep(0, s - 1), -Theta[i]))
          }
        }
      }else Q.pol = 1

      if(!is.null(xreg)) {
        its = its + Q
        nreg = ncol(as.matrix(xreg))
        xreg.coef = coefficients[(its + 1):(its + nreg)]
        if(nreg > 1){
          for (i in 1:nreg) {
            y = y - xreg.coef[i]*xreg[, i]
          }
        }else{
          y = y - xreg.coef*xreg
        }
      }

      if(intercep){
        for (i in 1:ncol) {
          if(a[i] > 0){
            for (j in 1:a[i]) {
              c = c*(1 - delta[[i]][j])
            }
          }
        }
        if(p > 0){
          for (i in 1:p) {
            c = c*(1 - phi[i])
          }
        }
        if(P > 0){
          for (i in 1:P) {
            c = c*(1 - Phi[i])
          }
        }
      }

      if(P > 0 || p > 0){
        if(P > 0 && p > 0){
          d.coef = d.acum*P.pol*p.pol*diff.pol
          d.ret = NULL
          for (i in 2:length(d.coef)) {
            if(d.coef[i] != 0){
              d.ret = c(d.ret, i - 1)
            }
          }
          d.coef = -d.coef

          o.coef = list()
          o.ret = list()
          for (i in 1:ncol) {
            o.coef[[i]] = o.pol[[i]]
            for (j in 1:ncol) {
              if(i != j){
                o.coef[[i]] = o.coef[[i]]*d.pol[[j]]
              }
            }
            o.coef[[i]] = o.coef[[i]]*P.pol*p.pol*diff.pol
            tam = 0
            for (j in 1:length(o.coef[[i]])) {
              if(o.coef[[i]][j] != 0) tam = tam + 1
            }
            o.ret[[i]] = rep(0, tam)
            cont = 1
            for (j in 1:length(o.coef[[i]])) {
              if(o.coef[[i]][j] != 0){
                o.ret[[i]][cont] = j - 1
                cont = cont + 1
              }
            }
          }

        }else{
          if(P > 0){
            d.coef = d.acum*P.pol*diff.pol
            d.ret = NULL
            for (i in 2:length(d.coef)) {
              if(d.coef[i] != 0){
                d.ret = c(d.ret, i - 1)
              }
            }
            d.coef = -d.coef

            o.coef = list()
            o.ret = list()
            for (i in 1:ncol) {
              o.coef[[i]] = o.pol[[i]]
              for (j in 1:ncol) {
                if(i != j){
                  o.coef[[i]] = o.coef[[i]]*d.pol[[j]]
                }
              }
              o.coef[[i]] = o.coef[[i]]*P.pol*diff.pol
              tam = 0
              for (j in 1:length(o.coef[[i]])) {
                if(o.coef[[i]][j] != 0) tam = tam + 1
              }
              o.ret[[i]] = rep(0, tam)
              cont = 1
              for (j in 1:length(o.coef[[i]])) {
                if(o.coef[[i]][j] != 0){
                  o.ret[[i]][cont] = j - 1
                  cont = cont + 1
                }
              }
            }

          }else{
            d.coef = d.acum*p.pol*diff.pol
            d.ret = NULL
            for (i in 2:length(d.coef)) {
              if(d.coef[i] != 0){
                d.ret = c(d.ret, i - 1)
              }
            }
            d.coef = -d.coef

            o.coef = list()
            o.ret = list()
            for (i in 1:ncol) {
              o.coef[[i]] = o.pol[[i]]
              for (j in 1:ncol) {
                if(i != j){
                  o.coef[[i]] = o.coef[[i]]*d.pol[[j]]
                }
              }
              o.coef[[i]] = o.coef[[i]]*p.pol*diff.pol

              for (j in 1:length(o.coef[[i]])) {
                if(o.coef[[i]][j] != 0) tam = tam + 1
              }
              o.ret[[i]] = rep(0, tam)
              cont = 1
              for (j in 1:length(o.coef[[i]])) {
                if(o.coef[[i]][j] != 0){
                  o.ret[[i]][cont] = j - 1
                  cont = cont + 1
                }
              }
            }
          }
        }
      }else{
        d.coef = d.acum*diff.pol
        if(length(d.acum) > 1){
          d.ret = NULL
          for (i in 2:length(d.coef)) {
            if(d.coef[i] != 0){
              d.ret = c(d.ret, i - 1)
            }
          }
        }else d.ret = 0
        d.coef = -d.coef

        o.coef = list()
        o.ret = list()
        for (i in 1:ncol) {
          o.coef[[i]] = o.pol[[i]]
          for (j in 1:ncol) {
            if(i != j){
              o.coef[[i]] = o.coef[[i]]*d.pol[[j]]
            }
          }
          o.coef[[i]] = o.coef[[i]]*diff.pol

          for (j in 1:length(o.coef[[i]])) {
            if(o.coef[[i]][j] != 0) tam = tam + 1
          }
          o.ret[[i]] = rep(0, tam)
          cont = 1
          for (j in 1:length(o.coef[[i]])) {
            if(o.coef[[i]][j] != 0){
              o.ret[[i]][cont] = j - 1
              cont = cont + 1
            }
          }
        }
      }

      if(Q > 0 || q > 0){
        if(Q > 0 && q > 0){
          q.coef = d.acum*Q.pol*q.pol
          q.ret = NULL
          for (i in 2:length(q.coef)) {
            if(q.coef[i] != 0){
              q.ret = c(q.ret, i - 1)
            }
          }
        }else{
          if(Q > 0){
            q.coef = d.acum*Q.pol
            q.ret = NULL
            for (i in 2:length(q.coef)) {
              if(q.coef[i] != 0){
                q.ret = c(q.ret, i - 1)
              }
            }
          }else{
            q.coef = d.acum*q.pol
            q.ret = 1:(length(q.coef) - 1)
          }
        }
      }else{
        q.coef = d.acum
        if(length(d.acum) > 1){
          q.ret = 1:(length(q.coef) - 1)
        }else{
          q.ret = 0
        }
      }


      forecast = NULL
      var.fore = NULL
      lower.95 = NULL
      upper.95 = NULL
      lower.80 = NULL
      upper.80 = NULL
      y.prov = y
      y.fore = y.real
      y.lower.95 = y.real
      y.upper.95 = y.real
      y.lower.80 = y.real
      y.upper.80 = y.real
      x.fore = x
      res = object$residuals
      for (i in 1:h) {
        modelx = list()
        for (j in 1:ncol) {
          modelx[[j]] = forecast::Arima(ts(x.fore[, j], frequency = sx[j]), order = c(px[j], dx[j], qx[j]),
                                   seasonal = list(order = c(Px[j], Dx[j], Qx[j]), period = sx[j]),
                                   include.mean = TRUE, method = "CSS")

          f = forecast::forecast(modelx[[j]], h = 1)$mean
          if(j != 1){
            x.prov = cbind(x.prov, c(x.fore[, j], f))
          }else{
            x.prov = c(x.fore[, j], f)
          }
        }

        x.fore = x.prov
        fore = fore.2(y.prov, x.fore, res, c, d.coef, d.ret, o.coef, o.ret, q.coef,
                      q.ret, tf.order, arima.order, ncol)

        y.prov = c(y.prov, fore)

        if(det.seasonal || det.trend > 0){
          if(det.seasonal && det.trend > 0){

            fore.t = predict(det.t, h = 1)$mean
            det.t = c(det.t, fore.t)
            fore.s = det.s[seson + i]

            if(type.deter == "additive"){
              fore = fore + fore.s + fore.t
            }else fore = fore*fore.s*fore.t

          }else{
            if(det.seasonal){

              fore.s = det.s[seson + i]

              if(type.deter == "additive"){
                fore = fore + fore.s
              }else fore = fore*fore.s
            }else{

              fore.t = predict(det.t, h = 1)$mean
              det.t = c(det.t, fore.t)

              if(type.deter == "additive"){
                fore = fore + fore.t
              }else fore = fore*fore.t
            }
          }
        }

        y.fore = c(y.fore, fore)
        forecast = c(forecast, fore)

        n.var = ar.ols(x = res, order.max = i, aic = FALSE, intercept = FALSE, demean = FALSE)$ar
        w.var = list()
        for (j in 1:ncol) {
          w.var[[j]] = ar.ols(x = modelx[[j]]$residuals, order.max = i, aic = FALSE, intercept = FALSE,
                       demean = FALSE)$ar
        }

        if(i > 1){
          n.var = var(res)*sum(n.var^2)
          for (j in 1:ncol) {
            w.var[[j]] = var(modelx[[j]]$residuals)*sum(w.var[[j]]^2)
          }
        }else{
          n.var = var(res)*n.var^2
          for (j in 1:ncol) {
            w.var[[j]] = var(modelx[[j]]$residuals)*w.var[[j]]^2
          }
        }

        sw.var = w.var[[1]]
        for (j in 2:ncol) {
          sw.var = sw.var + w.var[[j]]
        }
        var = sw.var + n.var
        var.fore = c(var.fore, var)

        lower.95 = c(lower.95, fore - 1.959964*sqrt(var.fore[i]))
        upper.95 = c(upper.95, fore + 1.959964*sqrt(var.fore[i]))
        lower.80 = c(lower.80, fore - 1.281552*sqrt(var.fore[i]))
        upper.80 = c(upper.80, fore + 1.281552*sqrt(var.fore[i]))

        y.lower.95 = c(y.lower.95, lower.95[i])
        y.upper.95 = c(y.upper.95, upper.95[i])
        y.lower.80 = c(y.lower.80, lower.80[i])
        y.upper.80 = c(y.upper.80, upper.80[i])
      }


    }else{

      if(!missing(arima.x)){
        if(length(arima.x) != 3) stop("Elemento arima.x mal especificado")
        if(length(Sarima.x) != 4) stop("Elemento Sarima.x mal especificado")
        px = arima.x[1]
        dx = arima.x[2]
        qx = arima.x[3]
        Px = Sarima.x[1]
        Dx = Sarima.x[2]
        Qx = Sarima.x[3]
        sx = Sarima.x[4]
        if(px < 0) stop("arima.x mal especificado")
        if(qx < 0) stop("arima.x mal especificado")
        if(Px < 0) stop("arima.x mal especificado")
        if(Qx < 0) stop("arima.x mal especificado")
        if (dx < 0) stop("Los ordenes de differenciacion deben ser > 0.")
        if (sx <= 0) stop("Componente (s) mal especificado")
        if (Dx < 0) stop("Los ordenes de differenciacion deben ser > 0.")
      }else{
        arim = forecast::auto.arima(x, method = "CSS", d = 2, max.D = 1,
                            max.p = 3, max.P = 2, max.q = 3, max.Q = 2)
        arim = arim$arma
        px = arim[1]
        dx = arim[6]
        qx = arim[2]
        Px = arim[3]
        Dx = arim[7]
        Qx = arim[4]
        sx = arim[5]
      }

      a = object$model$tf[1]
      b = object$model$tf[2]
      m = object$model$tf[3]
      tf.order = c(a, b, m)

      if(is.null(object$coefficients$Intercept)){
        coefficients = c(0, coefficients)
        c = 0
        intercep = FALSE
      }else{
        c = coefficients[1]
        intercep = TRUE
      }

      omega = coefficients[2:(m + 2)]
      if(m > 0) {
        o.pol = polynomial(omega)
      }else{
        o.pol = omega
      }

      if(a > 0){
        delta = coefficients[(3 + m):(2 + m + a)]
        d.pol = polynomial(c(1, -delta))
      }else{
        delta = 1
        d.pol = delta
      }

      if(p > 0){
        phi = coefficients[(3 + m + a):(2 + m + a + p)]
        p.pol = polynomial(c(1, -phi))
      }else p.pol = 1

      if(q > 0) {
        theta = coefficients[(3 + m + a + p):(2 + a + m + p + q)]
        q.pol = polynomial(c(1, -theta))
      }else q.pol = 1

      if(P > 0) {
        Phi = coefficients[(3 + a + m + p + q):(2 + a + m + p + q + P)]
        P.pol = polynomial(c(1, rep(0, s - 1), -Phi[1]))
        if(P > 1){
          for (i in 2:P) {
            P.pol = polynomial(c(P.pol, rep(0, s - 1), -Phi[i]))
          }
        }
      }else P.pol = 1

      if(Q > 0) {
        Theta = coefficients[(3 + a + m + p + q + P):(2 + a + m + p + q + P + Q)]
        Q.pol = polynomial(c(1, rep(0, s - 1), -Theta[1]))
        if(Q > 1){
          for (i in 2:Q) {
            Q.pol = polynomial(c(Q.pol, rep(0, s - 1), -Theta[i]))
          }
        }
      }else Q.pol = 1

      if(!is.null(xreg)) {
        nreg = ncol(as.matrix(xreg))
        xreg.coef = coefficients[(3 + a + m + p + q + P + Q):(2 + a + m + p + q + P + Q + nreg)]
        if(nreg > 1){
          for (i in 1:nreg) {
            y = y - xreg.coef[i]*xreg[, i]
          }
        }else{
          y = y - xreg.coef*xreg
        }
      }

      if(intercep){
        if(a > 0){
          for (i in 1:a) {
            c = c*(1 - delta[i])
          }
        }
        if(p > 0){
          for (i in 1:p) {
            c = c*(1 - phi[i])
          }
        }
        if(P > 0){
          for (i in 1:P) {
            c = c*(1 - Phi[i])
          }
        }
      }

      if(P > 0 || p > 0){
        if(P > 0 && p > 0){
          d.coef = d.pol*P.pol*p.pol*diff.pol
          d.ret = NULL
          for (i in 2:length(d.coef)) {
            if(d.coef[i] != 0){
              d.ret = c(d.ret, i - 1)
            }
          }
          d.coef = -d.coef

          o.coef = o.pol*P.pol*p.pol*diff.pol
          o.ret = NULL
          for (i in 1:length(o.coef)) {
            if(o.coef[i] != 0){
              o.ret = c(o.ret, i - 1)
            }
          }
        }else{
          if(P > 0){
            d.coef = P.pol*d.pol*diff.pol
            d.ret = NULL
            for (i in 2:length(d.coef)) {
              if(d.coef[i] != 0){
                d.ret = c(d.ret, i - 1)
              }
            }
            d.coef = -d.coef

            o.coef = P.pol*o.pol*diff.pol
            o.ret = NULL
            for (i in 1:length(o.coef)) {
              if(o.coef[i] != 0){
                o.ret = c(o.ret, i - 1)
              }
            }
          }else{
            d.coef = p.pol*d.pol*diff.pol
            d.ret = NULL
            for (i in 2:length(d.coef)) {
              if(d.coef[i] != 0){
                d.ret = c(d.ret, i - 1)
              }
            }
            d.coef = -d.coef

            o.coef = p.pol*o.pol*diff.pol
            o.ret = NULL
            for (i in 1:length(o.coef)) {
              if(o.coef[i] != 0){
                o.ret = c(o.ret, i - 1)
              }
            }
          }
        }
      }else{
        d.coef = d.pol*diff.pol
        if(length(d.coef) > 1){
          d.ret = NULL
          for (i in 2:length(d.coef)) {
            if(d.coef[i] != 0){
              d.ret = c(d.ret, i - 1)
            }
          }
        }else d.ret = 0
        d.coef = -d.coef

        o.coef = o.pol*diff.pol
        o.ret = NULL
        for (i in 1:length(o.coef)) {
          if(o.coef[i] != 0){
            o.ret = c(o.ret, i - 1)
          }
        }
      }


      if(Q > 0 || q > 0){
        if(Q > 0 && q > 0){
          q.coef = Q.pol*q.pol*d.pol
          q.ret = NULL
          for (i in 2:length(q.coef)) {
            if(q.coef[i] != 0){
              q.ret = c(q.ret, i - 1)
            }
          }
        }else{
          if(Q > 0){
            q.coef = Q.pol*d.pol
            q.ret = NULL
            for (i in 2:length(q.coef)) {
              if(q.coef[i] != 0){
                q.ret = c(q.ret, i - 1)
              }
            }
          }else{
            q.coef = q.pol*d.pol
            q.ret = 1:(length(q.coef) - 1)
          }
        }
      }else{
        q.coef = d.pol
        if(a > 0){
          q.ret = 1:(length(q.coef) - 1)
        }else q.ret = 0

      }

      forecast = NULL
      var.fore = NULL
      lower.95 = NULL
      upper.95 = NULL
      lower.80 = NULL
      upper.80 = NULL
      y.prov = y
      y.fore = y.real
      y.lower.95 = y.real
      y.upper.95 = y.real
      y.lower.80 = y.real
      y.upper.80 = y.real
      x.fore = x
      res = object$residuals
      num = 1
      for (i in 1:h) {
        if(missing(x.predict)){
          modelx = forecast::Arima(x.fore, order = c(px, dx, qx),
                                   seasonal = list(order = c(Px, Dx, Qx), period = sx),
                                   include.mean = TRUE, method = "CSS")

          f = forecast::forecast(modelx, h = 1)$mean
        }else{
          f = x.predict[num]
          num = num +  1
        }

        x.fore = ts(c(x.fore, f))

        fore = fore.1(y.prov, x.fore, res, c, d.coef, d.ret, o.coef, o.ret, q.coef,
                      q.ret, tf.order, arima.order)

        y.prov = c(y.prov, fore)

        if(det.seasonal || det.trend > 0){
          if(det.seasonal && det.trend > 0){

            fore.t = predict(det.t, h = 1)$mean
            det.t = c(det.t, fore.t)
            fore.s = det.s[seson + i]

            if(type.deter == "additive"){
              fore = fore + fore.s + fore.t
            }else fore = fore*fore.s*fore.t

          }else{
            if(det.seasonal){

              fore.s = det.s[seson + i]

              if(type.deter == "additive"){
                fore = fore + fore.s
              }else fore = fore*fore.s
            }else{

              fore.t = predict(det.t, h = 1)$mean
              det.t = c(det.t, fore.t)

              if(type.deter == "additive"){
                fore = fore + fore.t
              }else fore = fore*fore.t
            }
          }
        }

        y.fore = c(y.fore, fore)
        forecast = c(forecast, fore)

        if(missing(x.predict)){
          n.var = ar.ols(x = res, order.max = i, aic = FALSE, intercept = FALSE, demean = FALSE)$ar
          w.var = ar.ols(x = modelx$residuals, order.max = i, aic = FALSE, intercept = FALSE,
                         demean = FALSE)$ar

          if(i > 1){
            n.var = var(res)*sum(n.var^2)
            w.var = var(modelx$residuals)*sum(w.var^2)
          }else{
            n.var = var(res)*n.var^2
            w.var = var(modelx$residuals)*w.var^2
          }

          var = w.var + n.var
          var.fore = c(var.fore, var)
          lower.95 = c(lower.95, fore - 1.959964*sqrt(var.fore[i]))
          upper.95 = c(upper.95, fore + 1.959964*sqrt(var.fore[i]))
          lower.80 = c(lower.80, fore - 1.281552*sqrt(var.fore[i]))
          upper.80 = c(upper.80, fore + 1.281552*sqrt(var.fore[i]))

          y.lower.95 = c(y.lower.95, lower.95[i])
          y.upper.95 = c(y.upper.95, upper.95[i])
          y.lower.80 = c(y.lower.80, lower.80[i])
          y.upper.80 = c(y.upper.80, upper.80[i])
        }
      }

    }

  }else stop('object debe ser un elemento de tipo tf.model')

  if(missing(x.predict)){
    if(plot) plot(tf.tsplot(y.upper.95, y.lower.95, y.upper.80, y.lower.80, y.fore))
    return(invisible(list(y.fore = y.fore, y.upper.95 = y.upper.95, y.upper.80 = y.upper.80,
                          y.lower.95 = y.lower.95, y.lower.80 = y.lower.80, forecast = forecast,
                          var.fore = var.fore, upper.95 = upper.95, upper.80 = upper.80,
                          lower.95 = lower.95, lower.80 = lower.80, y = y.real, x = x)))
  }else{
    if(plot) plot(tf.tsplot(y.fore))
    return(invisible(list(y.fore = y.fore, forecast = forecast, y = y.real, x = x)))
  }

}


