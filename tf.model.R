#' @export tf.estimation
tf.model = function(y, x, order = c(a, b, m),
                    res.order = c(p = 0, d = 0, q = 0),
                    res.Sorder = c(P = 0, D = 0, Q = 0, s = NA),
                    include.mean = FALSE, xreg = NULL,
                    y.seasonal = FALSE, y.trend = 0, y.type.deter,
                    summary = TRUE)
{

  p = res.order[1]
  d = res.order[2]
  q = res.order[3]
  P = res.Sorder[1]
  D = res.Sorder[2]
  Q = res.Sorder[3]
  s = res.Sorder[4]
  if (is.na(s)) s = frequency(y)
  if(p < 0) stop("res.order mal especificado")
  if(q < 0) stop("res.order mal especificado")
  if(P < 0) stop("res.Sorder mal especificado")
  if(Q < 0) stop("res.Sorder mal especificado")
  if (d < 0) stop("Los ordenes de differenciacion deben ser > 0.")
  if (s <= 0) stop("Componente (s) mal especificado")
  if (D < 0) stop("Los ordenes de differenciacion deben ser > 0.")
  arima.order = c(p,q,P,Q,s)


  if (missing(order)) stop("Ausencia del campo order")

  if(!is.null(xreg)) nreg = ncol(as.matrix(xreg))
  ncol = ncol(as.matrix(x))
  if((ncol > 1 && !is.list(order)) || (ncol == 1 && is.list(order))){
    stop("Variable mal especificada: La variable order debe se de tipo list si ncol(x) > 1 y tener la misma cantidad de vectores que x")
  }
  if(is.list(order) && ncol != length(order)) stop("Tamaño de variable order incompatible con ncol(x)")

  if (!is.ts(y)) stop("y debe ser una serie de tiempo")

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

  xreg.real = xreg
  x.real = x
  y.real = y
  y.detr = y.real
  n.real = length(y.real)

  if(!is.null(xreg)){
    nreg = ncol(as.matrix(xreg))
    if(nreg > 1){
      if(nrow(xreg) != n.real) stop("Tamaño de xreg no compatible")
    }else{
      if(length(xreg) != n.real) stop("Tamaño de xreg no compatible")
    }
  }

  if(y.seasonal || y.trend > 0){
    if(missing(y.type.deter)) y.type.deter = "additive"
    if(y.trend < 0 || !is.numeric(y.trend)) stop("Variable det.trend mal especificada")

    if(y.seasonal && y.trend > 0){

      if (frequency(y) <= 1 || length(na.omit(y)) < 2*frequency(y)) stop("La serie de tiempo no tiene frequencia mayor o igual 2")
      det.s = tf.decompose(y, type = y.type.deter)$seasonal
      X = NULL
      for (i in (0:y.trend)) X = cbind(X, (1:n.real)^i)
      ind = solve.qr(qr(X), y)
      det.t = as.vector(X%*%ind)

      if(y.type.deter == "additive"){
        y = y - det.s - det.t
      }else y = y/(det.s*det.t)

    }else{
      if(y.seasonal){

        if (frequency(y) <= 1 || length(na.omit(y)) < 2*frequency(y)) stop("La serie de tiempo no tiene frequencia mayor o igual 2")
        det.s = tf.decompose(y, type = y.type.deter)$seasonal

        if(y.type.deter == "additive"){
          y = y - det.s
        }else y = y/det.s

      }else{

        X = NULL
        for (i in (0:y.trend)) X = cbind(X, (1:n.real)^i)
        ind = solve.qr(qr(X), y)
        det.t = as.vector(X%*%ind)

        if(y.type.deter == "additive"){
          y = y - det.t
        }else y = y/det.t
      }
    }
    y.detr = y
  }

  if (!is.ts(x)) x = ts(x)
  if (d > 0 || D > 0){
    if(d > 0 && D > 0){

      if(!is.null(xreg)) xreg = diff(xreg, differences = d)
      if(!is.null(xreg)) xreg = diff(xreg, differences = D, lag = s)
      x = diff(x, differences = d)
      x = diff(x, differences = D, lag = s)
      y = diff(y, differences = d)
      y = diff(y, differences = D, lag = s)

      if(d > 1){
        coef.diff = polynom::polynomial(c(1, -1))
        pol.prov = coef.diff
        for (i in 2:d) {
          coef.diff = coef.diff*pol.prov
        }
      }else{
        coef.diff = polynom::polynomial(c(1, -1))
      }
      if(D > 1){
        coef.Sdiff = polynom::polynomial(c(1, rep(0, s - 1), -1))
        pol.prov = coef.Sdiff
        for (i in 2:D) {
          coef.Sdiff = coef.Sdiff*pol.prov
        }
      }else{
        coef.Sdiff = polynom::polynomial(c(1, rep(0, s - 1), -1))
      }

      coef = coef.diff*coef.Sdiff

      pot = NULL
      for (i in 1:length(coef)) {
        if(coef[i] != 0){
          pot = c(pot, i - 1)
        }
      }

      pot = pot[2:length(pot)]
      ret.max = max(pot) + 1

      filter = rep(0, length(y))
      for (i in pot) {
        filter = filter + coef[i + 1]*y.real[(ret.max - i):(n.real - i)]
      }

    }else{
      if(d > 0){

        if(!is.null(xreg)) xreg = diff(xreg, differences = d)
        x = diff(x, differences = d)
        y = diff(y, differences = d)

        if(d > 1){
          coef = polynom::polynomial(c(1, -1))
          pol.prov = coef
          for (i in 2:d) {
            coef = coef*pol.prov
          }
        }else{
          coef = polynom::polynomial(c(1, -1))
        }

        filter = rep(0, length(y))
        for (i in 1:(length(coef) - 1)) {
          filter = filter + coef[i + 1]*y.real[(length(coef) - i):(n.real - i)]
        }

      }else{

        if(!is.null(xreg)) xreg = diff(xreg, differences = D, lag = s)
        x = diff(x, differences = D, lag = s)
        y = diff(y, differences = D, lag = s)

        if(D > 1){
          coef = polynom::polynomial(c(1, rep(0, s - 1), -1))
          pol.prov = coef
          for (i in 2:D) {
            coef = coef*pol.prov
          }
        }else{
          coef = polynom::polynomial(c(1, rep(0, s - 1), -1))
        }

        tam = 0
        for (i in 1:length(coef)) {
          if(coef[i] != 0){
            tam = tam + 1
            ret.max = i
          }
        }
        filter = rep(0, length(y))
        for (i in 1:(tam - 1)) {
          filter = filter + coef[i*s + 1]*y.real[(ret.max - i*s):(length(y.real) - i*s)]
        }
      }
    }
  }

  n = length(y)

  #FUNCIONEEEEEEEEES

  funcionCSS = function(init.par, y, x, tf.order, arima.order, ncol, include.mean,xreg , n){
    if(ncol > 1){
      res = residual.2(init.par, y, x, tf.order, arima.order, ncol, include.mean, xreg, n)
    }else{
      res = residual.1(init.par, y, x, tf.order, arima.order, ncol, include.mean, xreg, n)
    }
    se = sqrt(var(res))
    return(-sum(dnorm(res, mean = rep(0, length(res)), sd = se, log = TRUE)))
  }

  residual.1 = function(init.par, y, x, tf.order, arima.order, ncol, include.mean, xreg, n){

    a = tf.order[1]
    b = tf.order[2]
    m = tf.order[3]
    p = arima.order[1]
    q = arima.order[2]
    P = arima.order[3]
    Q = arima.order[4]
    s = arima.order[5]

    if(include.mean){
      y = y - init.par[1]
    }else{
      init.par = c(0 ,init.par)
    }

    omega = init.par[2:(m + 2)]
    if(a > 0){
      delta = init.par[(3 + m):(2 + m + a)]
    }

    if(p > 0) phi = init.par[(3 + m + a):(2 + m + a + p)]
    if(q > 0) theta = init.par[(3 + m + a + p):(2 + a + m + p + q)]
    if(P > 0) Phi = init.par[(3 + a + m + p + q):(2 + a + m + p + q + P)]
    if(Q > 0) Theta = init.par[(3 + a + m + p + q + P):(2 + a + m + p + q + P + Q)]
    if(!is.null(xreg)) {
      nreg = ncol(as.matrix(xreg))
      if(nreg > 1){
        for (i in 1:nreg) {
          reg = init.par[2 + a + m + p + q + P + Q + i]
          y = y - reg*xreg[, i]
        }
      }else{
        reg = init.par[3 + a + m + p + q + P + Q]
        y = y - reg*xreg
      }
    }

    if(a > 0){
      w.res = filter(x, delta, method = "r", init = rep(mean(x), a))
    }else{
      w.res = x
    }

    ist = b + m + 1
    n.res = length(w.res)
    d.res = y[ist:n.res] - omega[1]*w.res[(ist - b):(n.res - b)]
    if (m > 0) {
      for (i in 1:m) {
        d.res = d.res - omega[i + 1]*w.res[(ist - i - b):(n.res - i - b)]
      }
    }

    n.res = length(d.res)
    p.res = d.res[(p + 1):n.res]
    if (p > 0) {
      for (i in 1:p) {
        p.res = p.res - phi[i]*d.res[(p + 1 - i):(n.res - i)]
      }
    }
    if (q > 0){
      p.res = filter(p.res, theta, method = "r", init = rep(0, q))
    }

    n.res = length(p.res)
    P.res = p.res[(P*s + 1):n.res]
    if (P > 0) {
      for (i in 1:P) {
        P.res = P.res - Phi[i]*p.res[(P*s + 1 - i*s):(n.res - i*s)]
      }
    }
    if (Q > 0) {
      f1 = rep(0, s*Q)
      for (i in 1:Q) {
        f1[i*s] = Theta[i]
      }
      P.res = filter(P.res, f1, method = "r", init = rep(0, s*Q))
    }

    res = P.res
    return(res)
  }

  fitted.1 = function(y, x, coefficients, tf.order, arima.order, include.mean, xreg, res){

    a = tf.order[1]
    b = tf.order[2]
    m = tf.order[3]
    p = arima.order[1]
    q = arima.order[2]
    P = arima.order[3]
    Q = arima.order[4]
    s = arima.order[5]
    n = length(y)

    if(include.mean){
      y = y - coefficients[1]
    }else{
      coefficients = c(0 ,coefficients)
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
    }
    if(q > 0) {
      theta = coefficients[(3 + m + a + p):(2 + a + m + p + q)]
      q.pol = polynomial(c(1, -theta))
    }
    if(P > 0) {
      Phi = coefficients[(3 + a + m + p + q):(2 + a + m + p + q + P)]
      P.pol = polynomial(c(1, rep(0, s - 1), -Phi[1]))
      if(P > 1){
        for (i in 2:P) {
          P.pol = polynomial(c(P.pol, rep(0, s - 1), -Phi[i]))
        }
      }
    }
    if(Q > 0) {
      Theta = coefficients[(3 + a + m + p + q + P):(2 + a + m + p + q + P + Q)]
      Q.pol = polynomial(c(1, rep(0, s - 1), -Theta[1]))
      if(Q > 1){
        for (i in 2:Q) {
          Q.pol = polynomial(c(Q.pol, rep(0, s - 1), -Theta[i]))
        }
      }
    }
    if(!is.null(xreg)){
      nreg = ncol(as.matrix(xreg))
      reg = coefficients[(3 + a + m + p + q + P + Q):(2 + a + m + p + q + P + Q + nreg)]
    }

    if(P > 0 || p > 0){
      if(P > 0 && p > 0){
        d.coef = d.pol*P.pol*p.pol
        d.ret = NULL
        for (i in 2:length(d.coef)) {
          if(d.coef[i] != 0){
            d.ret = c(d.ret, i - 1)
          }
        }
        d.coef = -d.coef

        o.coef = o.pol*P.pol*p.pol
        o.ret = NULL
        for (i in 1:length(o.coef)) {
          if(o.coef[i] != 0){
            o.ret = c(o.ret, i - 1)
          }
        }
      }else{
        if(P > 0){
          d.coef = P.pol*d.pol
          d.ret = NULL
          for (i in 2:length(d.coef)) {
            if(d.coef[i] != 0){
              d.ret = c(d.ret, i - 1)
            }
          }
          d.coef = -d.coef

          o.coef = P.pol*o.pol
          o.ret = NULL
          for (i in 1:length(o.coef)) {
            if(o.coef[i] != 0){
              o.ret = c(o.ret, i - 1)
            }
          }
        }else{
          d.coef = p.pol*d.pol
          d.coef = -d.coef
          d.ret = 1:(length(d.coef) - 1)

          o.coef = p.pol*o.pol
          o.ret = 0:(length(o.coef) - 1)
        }
      }
    }else{
      d.coef = -d.pol
      if(a > 0){
        d.ret = 1:(length(d.coef) - 1)
      }else d.ret = 0

      o.coef = o.pol
      o.ret = 0:(length(o.coef) - 1)
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

    ret.max = max(max(d.ret) + 1, max(o.ret) + b + 1, max(q.ret) + 1)
    n.max = length(y[ret.max:n])
    n.res = length(res)
    res = c(rep(mean(res), n - n.res), res)

    y.t = rep(0, n.max)
    if(a > 0 || p > 0 || P > 0){
      for (i in d.ret) {
        y.t = y.t + d.coef[i + 1]*y[(ret.max - i):(n - i)]
      }
    }

    x.t = rep(0, n.max)
    for (i in o.ret) {
      x.t = x.t + o.coef[i + 1]*x[(ret.max - b - i):(n - b - i)]
    }

    a.t = rep(0, n.max)
    if(a > 0 || q > 0 || Q > 0){
      for (i in q.ret) {
        a.t = a.t + q.coef[i + 1]*res[(ret.max - i):(n - i)]
      }
    }

    fit = y.t + x.t + a.t + coefficients[1]

    if(!is.null(xreg)){
      if(nreg > 1){
        for (i in 1:nreg) fit = fit + reg[i]*xreg[ret.max:n, i]
      }else fit = fit + reg*xreg[ret.max:n]
    }

    return(fit)
  }

  residual.2 = function(init.par, y, x, tf.order, arima.order, ncol, include.mean, xreg, n){

    a = tf.order[[1]]
    b = tf.order[[2]]
    m = tf.order[[3]]
    p = arima.order[1]
    q = arima.order[2]
    P = arima.order[3]
    Q = arima.order[4]
    s = arima.order[5]

    if(include.mean){
      y = y - init.par[1]
    }else{
      init.par = c(0 ,init.par)
    }

    omega = list()
    delta = list()

    omega[[1]] = init.par[2:(m[1] + 2)]
    its = m[1] + 2
    for (i in 2:ncol) {
      omega[[i]] = init.par[(its + 1):(its + m[i] + 1)]
      its = its + m[i] + 1
    }

    for (i in 1:ncol) {
      if(a[i] > 0){
        delta[[i]] = init.par[(its + 1):(its + a[i])]
        its = its + a[i]
      }
    }

    if (p > 0){
      phi = init.par[(its + 1):(its + p)]
      its = its + p
    }
    if (q > 0){
      theta = init.par[(its + 1):(its + q)]
      its = its + q
    }
    if (P > 0){
      Phi = init.par[(its + 1):(its + P)]
      its = its + P
    }
    if (Q > 0){
      Theta = init.par[(its + 1):(its + Q)]
    }

    if(!is.null(xreg)) {
      its = its + Q
      nreg = ncol(as.matrix(xreg))
      if(nreg > 1){
        for (i in 1:nreg) {
          reg = init.par[its + i]
          y = y - reg*xreg[, i]
        }
      }else{
        reg = init.par[its + 1]
        y = y - reg*xreg
      }
    }

    its = rep(0, ncol)
    for (i in 1:ncol) {
      its[i] = b[i] + m[i] + 1
    }
    its = max(its)

    d.res = y[its:n]
    for (i in 1:ncol) {

      w.res = x[, i]
      if (a[i] > 0) {
        w.res = filter(x[, i], delta[[i]], method = "r", init = rep(mean(x[, i]), a[i]))
      }

      n.res = length(w.res)
      d.res = d.res - omega[[i]][1]*w.res[(its - b[i]):(n - b[i])]
      if (m[i] > 0) {
        for (ii in 1:m[i]) {
          d.res = d.res - omega[[i]][ii + 1]*w.res[(its - ii - b[i]):(n.res - ii - b[i])]
        }
      }
    }

    n.res = length(d.res)
    p.res = d.res[(p + 1):n.res]
    if (p > 0) {
      for (i in 1:p) {
        p.res = p.res - phi[i]*d.res[(p + 1 - i):(n.res - i)]
      }
    }
    if (q > 0){
      p.res = filter(p.res, theta, method = "r", init = rep(0, q))
    }

    n.res = length(p.res)
    P.res = p.res[(P*s + 1):n.res]
    if (P > 0) {
      for (i in 1:P) {
        P.res = P.res - Phi[i]*p.res[(P*s + 1 - i*s):(n.res - i*s)]
      }
    }
    if (Q > 0) {
      f1 = rep(0, s*Q)
      for (i in 1:Q) {
        f1[i*s] = Theta[i]
      }
      P.res = filter(P.res, f1, method = "r", init = rep(0, s*Q))
    }

    res = P.res
    return(res)
  }

  fitted.2 = function(y, x, coefficients, tf.order, arima.order, include.mean, xreg, res){

    a = tf.order[[1]]
    b = tf.order[[2]]
    m = tf.order[[3]]
    p = arima.order[1]
    q = arima.order[2]
    P = arima.order[3]
    Q = arima.order[4]
    s = arima.order[5]
    n = length(y)

    if(include.mean){
      y = y - coefficients[1]
    }else{
      coefficients = c(0 ,coefficients)
    }

    omega = list()
    delta = list()
    o.pol = list()
    d.pol = list()

    omega[[1]] = coefficients[2:(m[1] + 2)]
    if(m[1] > 0){
      o.pol[[1]] = polynomial(omega[[1]])
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
    }
    if (q > 0){
      theta = coefficients[(its + 1):(its + q)]
      q.pol = polynomial(c(1, -theta))
      its = its + q
    }
    if (P > 0){
      Phi = coefficients[(its + 1):(its + P)]
      P.pol = polynomial(c(1, rep(0, s - 1), -Phi[1]))
      if(P > 1){
        for (i in 2:P) {
          P.pol = polynomial(c(P.pol, rep(0, s - 1), -Phi[i]))
        }
      }
      its = its + P
    }
    if (Q > 0){
      Theta = coefficients[(its + 1):(its + Q)]
      Q.pol = polynomial(c(1, rep(0, s - 1), -Theta[1]))
      if(Q > 1){
        for (i in 2:Q) {
          Q.pol = polynomial(c(Q.pol, rep(0, s - 1), -Theta[i]))
        }
      }
    }
    if(!is.null(xreg)){
      its = its + Q
      nreg = ncol(as.matrix(xreg))
      reg = coefficients[(its + 1):length(coefficients)]
    }

    if(P > 0 || p > 0){
      if(P > 0 && p > 0){
        d.coef = d.acum*P.pol*p.pol
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
          o.coef[[i]] = o.coef[[i]]*P.pol*p.pol
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
          d.coef = d.acum*P.pol
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
            o.coef[[i]] = o.coef[[i]]*P.pol
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
          d.coef = d.acum*p.pol
          d.coef = -d.coef
          d.ret = 1:(length(d.coef) - 1)

          o.coef = list()
          o.ret = list()
          for (i in 1:ncol) {
            o.coef[[i]] = o.pol[[i]]
            for (j in 1:ncol) {
              if(i != j){
                o.coef[[i]] = o.coef[[i]]*d.pol[[j]]
              }
            }
            o.coef[[i]] = o.coef[[i]]*p.pol
            o.ret[[i]] = 0:(length(o.coef[[i]]) - 1)
          }
        }
      }
    }else{
      d.coef = -d.acum
      if(d.acum == 1){
        d.ret = 0
      }else{
        d.ret = 1:(length(d.coef) - 1)
      }

      o.coef = list()
      o.ret = list()
      for (i in 1:ncol) {
        o.coef[[i]] = o.pol[[i]]
        for (j in 1:ncol) {
          if(i != j){
            o.coef[[i]] = o.coef[[i]]*d.pol[[j]]
          }
        }
        o.ret[[i]] = 0:(length(o.coef[[i]]) - 1)
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
      if(d.acum == 1){
        q.ret = 0
      }else{
        q.ret = 1:(length(q.coef) - 1)
      }
    }

    cont = NULL
    for (i in 1:ncol) {
      cont = c(cont, max(o.ret[[i]]) + b[i] + 1)
    }

    ret.max = max(max(d.ret) + 1, max(cont), max(q.ret) + 1)
    n.max = length(y[ret.max:n])
    n.res = length(res)
    res = c(rep(mean(res), n - n.res), res)

    n.a = 0
    for (i in 1:length(a)) {
      if(a[i] > 0) n.a = n.a + 1
    }

    y.t = rep(0, n.max)
    if(n.a > 0 || P > 0 || p > 0){
      for (i in d.ret) {
        y.t = y.t + d.coef[i + 1]*y[(ret.max - i):(n - i)]
      }
    }

    x.t = rep(0, n.max)
    for (i in 1:ncol) {
      for (j in o.ret[[i]]) {
        x.t = x.t + o.coef[[i]][j + 1]*x[(ret.max - b[i] - j):(n - b[i] - j), i]
      }
    }

    a.t = rep(0, n.max)
    if(n.a > 0 || Q > 0 || q > 0){
      for (i in q.ret) {
        a.t = a.t + q.coef[i + 1]*res[(ret.max - i):(n - i)]
      }
    }

    fit = y.t + x.t + a.t + coefficients[1]

    if(!is.null(xreg)){
      if(nreg > 1){
        for (i in 1:nreg) fit = fit + reg[i]*xreg[ret.max:n, i]
      }else fit = fit + reg*xreg[ret.max:n]
    }

    return(fit)
  }

  #FUNCIONEEEEEEEEES

  if(ncol > 1){

    a = order[[1]][1]
    b = order[[1]][2]
    m = order[[1]][3]
    for (i in 2:ncol) {
      a = c(a, order[[i]][1])
      b = c(b, order[[i]][2])
      m = c(m, order[[i]][3])
    }
    for (i in 1:ncol) {
      if(a[i] < 0) stop("order mal especificado")
      if(b[i] < 0) stop("order mal especificado")
      if(m[i] < 0) stop("order mal especificado")
    }
    tf.order = list(a, b, m)
    names(tf.order) = c("a", "b", "m")


    if (include.mean){
      init.par = mean(y)
    }else{
      init.par = NULL
    }

    for (i in 1:ncol) {
      init.par = c(init.par, rep(0.1, m[i] + 1))
    }
    for (i in 1:ncol) {
      if(a[i] > 0){
        init.par = c(init.par, rep(0.1, a[i]))
      }
    }

    if (p > 0) init.par = c(init.par, rep(0.1, p))
    if (q > 0) init.par = c(init.par, rep(0.01, q))
    if (P > 0) init.par = c(init.par, rep(0.01, P))
    if (Q > 0) init.par = c(init.par, rep(0.01, Q))
    if (!is.null(xreg)) init.par = c(init.par, rep(0.1, nreg))

    model = suppressWarnings(nlm(funcionCSS, init.par, hessian = TRUE,
                                 y, x, tf.order, arima.order, ncol, include.mean, xreg, n))
    coef = model$estimate

    omega = NULL
    delta = NULL
    for(i in 1:ncol){
      pas1 = paste0("omega.", i)
      omega = c(omega, paste(pas1, paste0(1:(m[i] + 1)), sep = "."))
      if(a[i] > 0){
        pas2 = paste0("delta.", i)
        delta = c(delta, paste(pas2, paste0(1:a[i]), sep = "."))
      }
    }

    coef.names = c(if(include.mean) "intercept", omega, delta,
                   if(p > 0) paste0("AR.", 1:p),
                   if(q > 0) paste0("MA.", 1:q), if(P > 0) paste0("sAR.", 1:P),
                   if(Q > 0) paste0("sMA.", 1:Q),
                   if(!is.null(xreg)) paste0("xreg.", 1:nreg))
    names(coef) = coef.names

    tf.numerador = list()
    tf.denominador = list()

    tf.numerador[[1]] = coef[2:(m[1] + 2)]
    its = m[1] + 2
    if(ncol > 1){
      for (i in 2:ncol) {
        tf.numerador[[i]] = coef[(its + 1):(its + m[i] + 1)]
        its = its + m[i] + 1
      }
    }

    for (i in 1:ncol) {
      if(a[i] > 0){
        tf.denominador[[i]] = coef[(its + 1):(its + a[i])]
        its = its + a[i]
      }else tf.denominador[[i]] = NA
    }
    if (p > 0){
      AR = coef[(its + 1):(its + p)]
      its = its + p
    }else AR = NA
    if (q > 0){
      MA = coef[(its + 1):(its + q)]
      its = its + q
    }else MA = NA
    if (P > 0){
      sAR = coef[(its + 1):(its + P)]
      its = its + P
    }else sAR = NA
    if (Q > 0){
      sMA = coef[(its + 1):(its + Q)]
    }else sMA = NA
    if(!is.null(xreg)) {
      its = its + Q
      reg = coef[(its + 1):(its + nreg)]
    }else reg = NA

    names(tf.numerador) = paste0("omega.", 1:ncol)
    names(tf.denominador) = paste0("delta.", 1:ncol)


    if(!include.mean){

      coefficients = list(coefficients = coef,
                          tf.numerador = tf.numerador,
                          tf.denominador = tf.denominador,
                          AR = AR, MA = MA, sAR = sAR, sMA = sMA, xreg = reg)
    }else{
      coefficients = list(coefficients = coef, Intercept = coef[1],
                          tf.numeradores = tf.numerador,
                          tf.denominadores = tf.denominador,
                          AR = AR, MA = MA, sAR = sAR, sMA = sMA, xreg = reg)
    }


    res = residual.2(coef, y, x, tf.order, arima.order, ncol, include.mean, xreg, n)
    fitted = fitted.2(y, x, coef, tf.order, arima.order, include.mean, xreg, res)

    n.res = length(res)
    n.fit = length(fitted)
    if(n.res > n.fit){
      dif = n.res - n.fit
      res[1:dif] = 0
      bool = TRUE
    }
    if(n.res < n.fit){
      dif = n.fit - n.res
      fitted[1:dif] = 0
      bool = FALSE
    }
    if(n.res == n.fit) bool = FALSE

    fitted = ts(c(y[1:(n - length(fitted))], fitted), frequency = s, start = start(y.real))
    if(d > 0 || D > 0){
      fitted = fitted - filter
      fitted = ts(c(y.detr[1:(n.real - length(fitted))], fitted), frequency = s, start = start(y.real))
    }

    if(y.seasonal || y.trend > 0){
      if(y.seasonal && y.trend > 0){
        if(y.type.deter == "additive"){
          fitted = fitted + det.s + det.t
        }else fitted = fitted*det.s*det.t
      }else{
        if(y.seasonal){
          if(y.type.deter == "additive"){
            fitted = fitted + det.s
          }else fitted = fitted*det.s
        }else{
          if(y.type.deter == "additive"){
            fitted = fitted + det.t
          }else fitted = fitted*det.t
        }
      }
    }

    var.res = var(res)
    se.res = sqrt(var(res))

    lik = sum(dnorm(res, mean = 0, sd = se.res, log = TRUE))
    k = length(coef)
    AIC = 2*k - 2*lik
    AICc = AIC + 2*k*(k + 1)/(n.fit - k - 1)
    BIC = -2*lik + k*log(n.fit)

    if(bool){
      res[1:dif] = 0
      res = as.ts(c(rep(0, n.real - n.res), res))
    }else{
      res = as.ts(c(rep(0, n.real - n.res), res))
    }

    cov.coef = solve(model$hessian)
    var.coef = diag(cov.coef)

    if(prod(var.coef) < 0){
      warning("Varianza negativa")
      corr.coef = NULL
      var.coef = NULL
      se.coef = NULL
      cov.coef = NULL
    }else{
      rownames(cov.coef) = coef.names
      colnames(cov.coef) = coef.names
      names(var.coef) = coef.names

      se.coef = sqrt(var.coef)
      names(se.coef) = coef.names
      corr.coef = cov.coef
      n.cov = ncol(cov.coef)
      for (i in 1:n.cov) {
        for (j in 1:n.cov) {
          corr.coef[i,j] = cov.coef[i, j]/(sqrt(var.coef[i])*sqrt(var.coef[j]))
        }
      }
      rownames(corr.coef) = coef.names
      colnames(corr.coef) = coef.names
    }

  }else{

    a = order[1]
    b = order[2]
    m = order[3]
    if(a < 0) stop("order mal especificado")
    if(b < 0) stop("order mal especificado")
    if(m < 0) stop("order mal especificado")
    tf.order = c(a, b, m)

    if (include.mean){
      init.par = mean(y)
    }else{
      init.par = NULL
    }

    init.par = c(init.par, rep(0.1, m + 1))
    if (a > 0) init.par = c(init.par, rep(0.1, a))
    if (p > 0) init.par = c(init.par, rep(0.1, p))
    if (q > 0) init.par = c(init.par, rep(0.01, q))
    if (P > 0) init.par = c(init.par, rep(0.01, P))
    if (Q > 0) init.par = c(init.par, rep(0.01, Q))
    if (!is.null(xreg)) init.par = c(init.par, rep(0.1, nreg))

    model = suppressWarnings(nlm(funcionCSS, init.par, hessian = TRUE,
                             y, x, tf.order, arima.order, ncol, include.mean, xreg, n))

    coef = model$estimate

    coef.names = c(if(include.mean) "intercept", paste0("omega.", 0:m),
                   if(a > 0) paste0("delta.", 1:a), if(p > 0) paste0("AR.", 1:p),
                   if(q > 0) paste0("MA.", 1:q), if(P > 0) paste0("sAR.", 1:P),
                   if(Q > 0) paste0("sMA.", 1:Q),
                   if(!is.null(xreg)) paste0("xreg.", 1:nreg))
    coefficients =
    names(coef) = coef.names

    tf.numerador = coef[2:(m + 2)]
    if(a > 0){
      tf.denominador = coef[(3 + m):(2 + m + a)]
    }else tf.denominador = NA
    if(p > 0){
      AR = coef[(3 + m + a):(2 + m + a + p)]
    }else AR = NA
    if(q > 0){
      MA = coef[(3 + m + a + p):(2 + a + m + p + q)]
    }else MA = NA
    if(P > 0){
      sAR = coef[(3 + a + m + p + q):(2 + a + m + p + q + P)]
    }else sAR = NA
    if(Q > 0){
      sMA = coef[(3 + a + m + p + q + P):(2 + a + m + p + q + P + Q)]
    }else sMA = NA
    if(!is.null(xreg)) {
      reg = coef[(3 + a + m + p + q + P + Q):(2 + a + m + p + q + P + Q + nreg)]
    }else reg = NA

    if(!include.mean){
      coefficients = list(coefficients = coef,
                          tf.numerador = tf.numerador,
                          tf.denominador = tf.denominador,
                          AR = AR, MA = MA, sAR = sAR, sMA = sMA, xreg = reg)
    }else{
      coefficients = list(coefficients = coef, Intercept = coef[1],
                          tf.numerador = tf.numerador,
                          tf.denominador = tf.denominador,
                          AR = AR, MA = MA, sAR = sAR, sMA = sMA, xreg = reg)
    }

    res = residual.1(coef, y, x, tf.order, arima.order, ncol, include.mean, xreg, n)
    fitted = fitted.1(y, x, coef, tf.order, arima.order, include.mean, xreg, res)

    n.res = length(res)
    n.fit = length(fitted)
    if(n.res > n.fit){
      dif = n.res - n.fit
      res[1:dif] = 0
      bool = TRUE
    }
    if(n.res < n.fit){
      dif = n.fit - n.res
      fitted[1:dif] = 0
      bool = FALSE
    }
    if(n.res == n.fit) bool = FALSE

    fitted = ts(c(y[1:(n - n.fit)], fitted), frequency = s, start = start(y.real))
    if(d > 0 || D > 0){
      fitted = fitted - filter
      fitted = ts(c(y.detr[1:(n.real - length(fitted))], fitted), frequency = s, start = start(y.real))
    }
    if(y.seasonal || y.trend > 0){
      if(y.seasonal && y.trend > 0){
        if(y.type.deter == "additive"){
          fitted = fitted + det.s + det.t
        }else fitted = fitted*det.s*det.t
      }else{
        if(y.seasonal){
          if(y.type.deter == "additive"){
            fitted = fitted + det.s
          }else fitted = fitted*det.s
        }else{
          if(y.type.deter == "additive"){
            fitted = fitted + det.t
          }else fitted = fitted*det.t
        }
      }
    }

    var.res = var(res)
    se.res = sqrt(var.res)

    lik = sum(dnorm(res, mean = 0, sd = se.res, log = TRUE))
    k = length(coef)
    AIC = 2*k - 2*lik
    AICc = AIC + 2*k*(k + 1)/(n.fit - k - 1)
    BIC = -2*lik + k*log(n.fit)

    if(bool){
      res[1:dif] = 0
      res = as.ts(c(rep(0, n.real - n.res), res))
    }else{
      res = as.ts(c(rep(0, n.real - n.res), res))
    }

    cov.coef = solve(model$hessian)
    var.coef = diag(cov.coef)

    if(prod(var.coef) < 0){
      warning("Varianza negativa")
      corr.coef = NULL
      var.coef = NULL
      se.coef = NULL
      cov.coef = NULL
    }else{
      rownames(cov.coef) = coef.names
      colnames(cov.coef) = coef.names
      names(var.coef) = coef.names
      se.coef = sqrt(var.coef)
      names(se.coef) = coef.names
      corr.coef = cov.coef
      n.cov = ncol(cov.coef)
      for (i in 1:n.cov) {
        for (j in 1:n.cov) {
          corr.coef[i,j] = cov.coef[i, j]/(sqrt(var.coef[i])*sqrt(var.coef[j]))
        }
      }
      rownames(corr.coef) = coef.names
      colnames(corr.coef) = coef.names
    }
  }

  t.statistic = coef/se.coef
  p.value = NULL
  for (i in 1:length(coef)) {
    p.value = c(p.value, pt(abs(t.statistic[i]), df = length(y.real) - length(coef),
                    lower.tail = FALSE)/2)
  }

  t.test = cbind(coef, se.coef, t.statistic, p.value)
  colnames(t.test) = c("coef", "s.e", "statistic", "p.value")
  if (summary){
    cat("Coefficients & s.e & t.test.:", "\n")
    print(t(t.test), digits = 3)
    cat("\n")
  }


  if(missing(y.type.deter)) y.type.deter = NA
  model = list(y = y.real, x = x.real, xreg = xreg.real, tf = tf.order, arima = c(p,d,q,P,D,Q,s),
               y.seasonal = y.seasonal, y.trend = y.trend, y.type.deter = y.type.deter)

  structure(list(coefficients = coefficients, fitted = fitted, var.coef = var.coef, se.coef = se.coef,
            cov.coef = cov.coef, corr.coef = corr.coef, residuals = res, var.res = var.res,
            se.res = se.res, t.test = t(t.test), log.likelihood = lik, AIC = AIC, AICc = AICc, BIC = BIC,
            model = model),
            class = "tf.model")

}
