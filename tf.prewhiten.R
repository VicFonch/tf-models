#' @title Estimación de la función respuestas a impulsos por preblanqueo
#' @description Esta función aproxima la variable dependiente por un proceso autorregresivo de
#' orden deseado, en caso de que no se identifique el órden, este se identificará a través del criterio
#' AIC. Calcula las series alpha y beta, de las cuales alpha es la serie de residuos del proceso  dicho
#' anteriormente,y beta es la serie temporal referente a la variable dependiente y, a la cual se le aplica
#' el filtro autorregresivo estimado con x y con estas series se estima la función de respuesta a
#' impulsos y a escalones. Además de esto se ofrece los valores de la función de correlaciones cruzadas y
#' datos acerca del modelo autorregresivo estimado en x
#' @param x Serie de tiempo que toma el lugar de variable expicativa o exógena
#' @param y Serie de tiempo que toma el lugar de variable dependiente o  endógena
#' @param fri.plot valor lógico: TRUE = hacer un grafico de barras para la funcion de respuesta
#' a impulsos
#' @param fre.plot valor lógico: TRUE = hacer un grafico de barras para la funcion de respuesta
#' a escalones
#' @param ccf.plot valor logico: TRUE = hacer un grafico de las correlaciones cruzadas entre las
#' variables preblanqueadas
#' @param maxorder orden máximo del polinomio autorregresivo
#' @param nfri cantidad de coeficientes de la funcion respuesa a impulsos y a escalones
#' @return Dependiendo de opciones :
#' Grafico de barras de la funcion respuesta a impulsos.
#' Grafico de barras de funcion respuesta a escalon.
#' Lista de resultados.
#' @export tf.prewhiten
#' @examples
#' x = ts(rnorm(100, 100, 12), 12, 1950)
#' y = ts(rnorm(100, 85, 15), 12, 1950)
#' p = tf.prewhiten(x, y)


tf.prewhiten = function(y, x, d = 0, D = 0, max.order, nfri = 20,
                        fri.plot = FALSE, fre.plot = FALSE, ccf.plot = FALSE)
{

  if(!is.ts(y)) stop("y debe ser una serie de tiempo")
  if(!is.ts(x)) stop("x debe ser una serie de tiempo")
  if (!missing(max.order) && max.order < 0) stop("maxorder debe ser > 0")
  if (nfri < 0) stop("nfri debe ser > 0")
  if (!is.logical(fri.plot)) stop("fri.plot debe ser un valor lógico")
  if (!is.logical(fre.plot)) stop("fre.plot debe ser un valor lógico")

  if(d > 0 || D > 0){
    s1 = frequency(y)
    s2 = frequency(x)
    if(s1 != s2){
      stop("Las series tienen frecuencias distintas")
    }else s = s1
    if (d < 0) stop("Los ordenes de differenciacion deben ser > 0.")
    if (D < 0) stop("Los ordenes de differenciacion deben ser > 0.")

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

  if (length(x) != length(y)){
    dif = abs(length(x) - length(y))
    if (length(x) > length(y)){
      x = x[1:(length(x) - dif)]
    }else{
      y = y[1:(length(y) - dif)]
    }
  }

  cont = 0
  desicion = c(fri.plot, fre.plot, ccf.plot)
  for(i in 1:length(desicion)) {
    if (desicion[i]){
      cont = cont + 1
      if (cont == 2) {
        stop("Esta funci�n solo soporta un gráfico")
      }
    }
  }

  if (!missing(max.order)) {
    modelo = stats::ar.ols(x = x, aic = FALSE ,order.max = max.order, demean = FALSE, intercept = FALSE)
  }else{
    modelo = stats::ar.ols(x = x, demean = FALSE, intercept = FALSE)
  }

  alpha = modelo$resid[(modelo$order + 1):length(modelo$resid)]
  beta = ts(rep(0,length(alpha)))
  coefi = c(1, -modelo$ar)

  for (i in 1:length(alpha)) {
    for (j in 1:length(coefi)) {
      beta[i] = beta[i] + coefi[j]*y[i + length(coefi) - j]
    }
  }

  alpha = ts(alpha)
  beta = ts(beta)

  if (ccf.plot){
    co.ccf = tf.ccf(alpha, beta, lag.max = nfri, plot = TRUE)$ccf[(nfri + 1):(2*nfri)]
    fri = as.vector(co.ccf*(sqrt(var(beta))/sqrt(var(alpha))))
  }else{
    co.ccf = tf.ccf(alpha, beta, lag.max = nfri, plot = FALSE)$ccf[(nfri + 1):(2*nfri)]
    fri = as.vector(co.ccf*(sqrt(var(beta))/sqrt(var(alpha))))
  }


  fre = c(rep(0, length(fri)))
  for (i in 1:length(fri)) {
    for (j in 1:i) {
      fre[i] = fre[i] + fri[j]
    }
  }

  data = data.frame(fri, fre)


  if (fri.plot){

    plot(ggplot(data = data, aes(y = fri, x = 1:length(fri))) + geom_bar(stat = "identity")+
           geom_hline(yintercept = 0)+
           labs(x = "coefficients",y = "FRI"))

  }
  if(fre.plot){

    plot(ggplot(data = data, aes(y = fre, x = 1:length(fre))) + geom_bar(stat = "identity")+
           geom_hline(yintercept = 0)+
           labs(x = "coefficients",y = "FRE"))


  }
  modelo = list(ar.order = modelo$order, n.serie = modelo$n.used, frequency = modelo$frequency,
                coef = modelo$ar, Coef.var = modelo$asy.se.coef$ar, residuals = modelo$resid)


  return(invisible(list(alpha = alpha, beta = beta ,fri = fri ,fre = fre ,ccf = co.ccf,
                        model = modelo)))


}



