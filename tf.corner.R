#' @export tf.corner
tf.corner = function (y, x, row = 11, col = 7, diff, Sdiff = c(D, s), fri)
{
  if ((length(row) != 1) || (length(col) != 1)) stop("El número de filas debe ser un valor entero")
  if (((floor(row) - row) != 0) || ((floor(col) - col) != 0)) stop("El número de columnas debe ser un valor entero")
  if (!missing(fri) && !is.numeric(fri)) stop("El vector fri debe ser un vector numerico")
  if(!is.ts(y)) stop("y debe ser una serie de tiempo")

  if(!is.ts(x)) x = ts(x)
  if (!missing(diff) || !missing(Sdiff)){
    if(!missing(diff) && !missing(Sdiff)){

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
      y = diff(y, differences = d)
      y = diff(y, differences = D, lag = s)

    }else{
      if(!missing(diff)){
        d = diff

        if (d <= 0) stop("Los ordenes de differenciacion deben ser > 0. Si no quire aplicar
                         differencias a las series no llene los parametros diff y Sdiff")

        x = diff(x, differences = d)
        y = diff(y, differences = d)
      }else{
        if (length(Sdiff) != 2) stop("El parametro Sdiff debe ser un vector de tamaño 2")

        D = Sdiff[1]
        s = Sdiff[2]

        if (s <= 0) stop("El componente s de Sdiff debe ser > 0")
        if (D <= 0) stop("Los ordenes de differenciacion deben ser > 0. Si no quire aplicar
                         differencias a las series no llene los parametros diff y Sdiff")

        x = diff(x, differences = D, lag = s)
        y = diff(y, differences = D, lag = s)
      }
    }
  }

  lag = row + col + 1

  if (length(x) != length(y)){
    dif = abs(length(x) - length(y))
    if (length(x) > length(y)){
      x = as.ts(x[1:(length(x) - dif)])
    }else{
      y = as.ts(y[1:(length(y) - dif)])
    }
  }
  leng.min = length(y)

  if (missing(fri)){
    v.i = tf.prewhiten(y, x, nfri = lag)$fri
    v.i = v.i/max(abs(v.i))
  }else{
    v.i = fri/max(abs(fri))
  }


  tbl = matrix(0, row, col)
  tbl[, 1] = v.i[1:row]

  for (j in 2:col) {
    for (i in 1:row) {
      cmx = diag(rep(v.i[i], j))
      for (ii in 2:j) {
        for (jj in 1:(ii - 1)) {
          idx = i + jj
          cmx[ii, jj] = v.i[idx]
        }
      }
      for (jj in 2:j) {
        for (ii in 1:(jj - 1)) {
          idx = i - jj + 1
          if (idx > 0)
            cmx[ii, jj] = v.i[idx]
        }
      }
      tbl[i, j] = det(cmx)
    }
  }

  for (i in 1:row) {
    for (j in 1:col) {
      if(is.na(tbl[i, j])){
        tbl[i,j] = 0
      }
    }
  }

  stbl = tbl
  crit = 2/sqrt(leng.min)
  tbl = cbind(c(1:row) - 1, tbl)
  colnames(tbl) = c("r->", paste(c(1:col)))
  cat("Tabla de Corner: ", "\n")
  print(round(tbl, 3))

  for (i in 1:row) {
    for (j in 1:col) {
      if (abs(tbl[i, j + 1]) <= crit) {
        stbl[i, j] = "O"
      }
      else {
        stbl[i, j] = "X"
      }
    }
  }

  cat("\n")
  cat("Tabla simplificada: ", "\n")
  stbl = cbind(c(1:row) - 1, stbl)
  colnames(stbl) = c("r->", paste(c(1:col)))
  print(stbl)
  return(invisible(list(Corner = tbl, Simply.Corner = stbl)))

}

