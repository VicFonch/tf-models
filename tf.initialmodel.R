#' @export tf.initialmodel
tf.initialmodel = function(y, x, order = c(a, b, m), diff, Sdiff, fricoef)
{
  if (missing(order)) stop("El vector order es necesario")
  if(!is.ts(y)) stop("y debe ser una serie de tiempo")
  ncol = ncol(as.matrix(x))

  if(ncol > 1 && !is.list(order) || ncol == 1 && is.list(order)){
    stop("Variable mal especificada:
         La variable order debe se de tipo list si ncol(x) > 1 y tener la misma cantidad de vectores que x")
  }

  if(missing(diff)) diff = 0
  if(missing(Sdiff)) Sdiff = c(0, frequency(y))

  if (diff != 0 || Sdiff[1] != 0){
    if(!is.ts(x)) x = ts(x)
    if(diff != 0 && Sdiff[1] != 0){

      d = diff
      D = Sdiff[1]
      s = Sdiff[2]

      if (d <= 0) stop("Los ordenes de differenciacion deben ser > 0. Si no quire aplicar
                         differencias a las series no llene los parametros diff y Sdiff")
      if (length(Sdiff) != 2) stop("El parametro Sdiff debe ser un vector de tamaño 2")
      if (s <= 0) stop("El componente s de Sdiff debe ser > 0")
      if (D <= 0) stop("Los ordenes de differenciacion deben ser > 0. Si no quire aplicar
                         differencias a las series no llene los parametros diff y Sdiff")

      x = diff(x, differences = d)
      x = diff(x, differences = D, lag = s)
      y.real = y
      y = diff(y, differences = d)
      y = diff(y, differences = D, lag = s)

    }else{
      if(diff != 0){

        d = diff

        if (d <= 0) stop("Los ordenes de differenciacion deben ser > 0. Si no quire aplicar
                         differencias a las series no llene los parametros diff y Sdiff")

        x = diff(x, differences = d)
        y.real = y
        y = diff(y, differences = d)

      }else{
        if (length(Sdiff) != 2) stop("El parametro Sdiff debe ser un vector de tamaño 2")

        D = Sdiff[1]
        s = Sdiff[2]

        if (s <= 0) stop("El componente s de Sdiff debe ser > 0")
        if (D <= 0) stop("Los ordenes de differenciacion deben ser > 0. Si no quire aplicar
                         differencias a las series no llene los parametros diff y Sdiff")

        x = base::diff(x, differences = D, lag = s)
        y.real = y
        y = base::diff(y, differences = D, lag = s)
      }
    }
  }

  if (ncol == 1){
    for (i in 1:3){
      if (order[i] < 0){
        stop("Valor negativo en el vector orden")
      }
    }
    a = order[1]
    b = order[2]
    m = order[3]

    if (length(order) != 3) stop("order debe ser un vector de tamaño 3")
    if (!is.numeric(order)) stop("order debe ser un vector numerico")
    if(!missing(fricoef)){
      if (!is.numeric(fricoef)) stop("fricoef debe ser un vector numérico")
      if (length(fricoef) < a + b + m + 1){
        stop("La cantidad de elementos en fricoef para estimar el modelo es insuficiente.
         De otra forma no llene el campo fricoef y el algoritmo se encargarÃ¡ del problema")
      }
    }
    if (!is.ts(x)) stop("x debe ser de tipo ts")

    if (length(x) != length(y)) {
      dif = abs(length(x) - length(y))
      if (length(x) > length(y)) {
        x = x[1:(length(x) - dif)]
      }else{
        y = y[1:(length(y) - dif)]
      }
    }
    n = length(y)

    if (missing(fricoef)){
      fri = suppressWarnings(tf.ltf(y, x, summary = FALSE,
                             res.order = c(1,0,0),res.Sorder = c(1,0,0,NA), fri.plot = FALSE,
                             fre.plot = FALSE)$v)
      fricoef = fri
    }

    if (b > 0) fricoef[1:b] = 0
    if(a > 0) fricoef = c(rep(0, a), fricoef)

    cant.param = a + m + 1

    matriz1 = matrix(0,cant.param, a)
    matriz2 = matrix(0,cant.param, m + 1)

    if(a > 0){
      cont = 2
      for (i in 2:cant.param) {
        for (j in 1:a) {
          matriz1[i, j] = fricoef[a + b + cont - j]
        }
        cont = cont + 1
      }
    }

    for (i in 1:(m + 1)) {
      for (j in 1:(cant.param)) {
        if (i == j) {
          matriz2[i,j] = 1
        }
      }
    }

    if(a > 0){
      matriz = cbind(matriz1, matriz2)
    }else matriz = matriz2
    v = fricoef[(b + a + 1):(b + a + ncol(matriz))]

    if (a != 0){
      ar.coef = solve(matriz,v)[1:a]
      ma.coef = solve(matriz,v)[(a + 1):(cant.param)]
    }else{
      ar.coef = NA
      ma.coef = solve(matriz, v)[1:(cant.param)]
    }

#AKkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
  }else{
    if (missing(order)) stop("El vector order es necesario")
    for (i in 1:ncol) {
      if (length(order[[i]]) != 3) stop("order debe contener solo vectores de tamÃ±o 3")
      if (!is.numeric(order[[i]])) stop("order debe contener solo vectores numericos")
    }

    max.order = order[[1]][1] + order[[1]][2] + order[[1]][3] + 1
    for (i in 2:ncol) {
      if (max.order < (order[[i]][1] + order[[i]][2] + order[[i]][3] + 1)){
          max.order = order[[i]][1] + order[[i]][2] + order[[i]][3] + 1
      }
    }

    if (!missing(fricoef)){
      if (!is.data.frame(fricoef)) stop("El elemento fricoef debe ser de tipo data.frame")
      if (ncol(fricoef) != ncol) stop("El elemento fricoef debe tener la misma cantidad de columnas
                                      que x")
      if (length(fricoef) < max.order){
        stop("La cantidad de elementos en fricoef para estimar el modelo es insuficiente.
              De otra forma no llene el campo fricoef y el algoritmo se encargarÃ¡ del problema")
      }
    }

    if (length(y) != nrow(x)){
      dif = abs(length(y) - nrow(x))
      if (length(y) < nrow(x)){
        x = x[1:(nrow(x) - dif), ]
      }else{
        y = y[1:(length(y) - dif)]
      }
    }
    n = length(y)

    if (missing(fricoef)){
      fri = suppressWarnings(tf.ltf(y, x, summary = FALSE,
                       res.order = c(1,0,0),res.Sorder = c(1,0,0,NA), fri.plot = FALSE,
                       fre.plot = FALSE))
      fri = fri$v
    }

    fricoef = list()
    for (i in 1:ncol) {
      if(order[[i]][1] > 0){
        fricoef[[i]] = c(rep(0, order[[i]][1]), fri[, i])
      }else{
        fricoef[[i]] =  fri[, i]
      }
    }

    cant.param = c(rep(0, ncol))
    matriz1 = list()
    matriz2 = list()
    matriz = list()
    ar.coef = list()
    ma.coef = list()
    for (i in 1:ncol) {
      if (order[[i]][2] > 0) fricoef[[i]][1:(order[[i]][1] + order[[i]][2])] = 0

      cant.param[i] = order[[i]][1] + order[[i]][3] + 1
      matriz1 = matrix(0,cant.param[i], order[[i]][1])
      matriz2 = matrix(0,cant.param[i], order[[i]][3] + 1)

      if(order[[i]][1] > 0){
        cont = 2
        for (ii in 2:cant.param[i]) {
          for (jj in 1:order[[i]][1]) {
            matriz1[ii, jj] = fricoef[[i]][order[[i]][1] + order[[i]][2] + cont - jj]
          }
          cont = cont + 1
        }
      }

      for (ii in 1:(order[[i]][3] + 1)) {
        for (jj in 1:cant.param[i]) {
          if (ii == jj) {
            matriz2[ii,jj] = 1
          }
        }
      }

      if(order[[i]][1] > 0){
        matriz[[i]] = cbind(matriz1, matriz2)
      }else matriz[[i]] = matriz2

      v = fricoef[[i]][(order[[i]][2] + order[[i]][1] + 1):(order[[i]][2] + order[[i]][1] + ncol(matriz[[i]]))]

      if (order[[i]][1] != 0){
        ar.coef[[i]] = solve(matriz[[i]],v)[1:order[[i]][1]]
        ma.coef[[i]] = solve(matriz[[i]],v)[(order[[i]][1] + 1):(cant.param[i])]
      }else{
        ar.coef[[i]] = NA
        ma.coef[[i]] = solve(matriz[[i]],v)[1:(cant.param[i])]
      }
    }

    names(ar.coef) = paste0("delta.", 1:ncol)
    names(ma.coef) = paste0("omega.", 1:ncol)
    names(matriz) = paste0("x.", 1:ncol)
  }

  return(invisible(list(Coef.ar = ar.coef, Coef.ma = ma.coef,
                        matriz = matriz, fricoef = fri)))

}

