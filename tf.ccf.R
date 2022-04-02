#' @export tf.ccf
tf.ccf = function(y, x, lag.max, seasonal, seasonal.auto = FALSE, prewhiten = FALSE,
                  arima.prewhiten = list(order = c(p,d,q), Sorder = c(P,D,Q,s)),
                  summary = FALSE, plot = TRUE)
{
  x = na.omit(x)
  y = na.omit(y)
  if(!is.ts(x)) stop("x debe ser una serie de tiempo")
  if(!is.ts(y)) stop("y debe ser una serie de tiempo")
  if (is.mts(x) || is.mts(y)) stop("Esta funcion es solo para series de tiempo univariantes")
  if (!missing(lag.max) && lag.max < 0) stop("lag.max debe ser > 0")
  if (!is.logical(seasonal.auto)) stop("seasonal debe ser un valor logico")
  if (!is.logical(plot)) stop("plot debe ser ")
  if (!missing(seasonal) && !is.numeric(seasonal)) stop("seasonal debe ser un valor numerico")
  if (!missing(seasonal) && length(seasonal) > 1) stop("seasonal de ver un numero positivo")
  if (!missing(seasonal) && seasonal < 0) stop("seasonal debe ser un numero positivo")

  if(prewhiten){
    start = start(x)
    frequency = frequency(x)
    if(missing(arima.prewhiten)){
      modelo = stats::ar.ols(x = x, demean = FALSE, intercept = TRUE, aic = TRUE)
      x = modelo$resid[(modelo$order + 1):length(modelo$resid)]
      y = y[(modelo$order + 1):length(y)]
      x = ts(x, start = start, frequency = frequency)
    }else{
      order = arima.prewhiten[[1]]
      Sorder = list(order = arima.prewhiten[[2]][1:3], period = arima.prewhiten[[2]][4])
      modelo = forecast::Arima(x, order = order, seasonal = Sorder, include.mean = TRUE,
                               method = "CSS")
      x = modelo$residuals
      x = ts(x, start = start, frequency = frequency)
    }
  }
  tam.min = min(length(x), length(y))

  if (!missing(seasonal) || (seasonal.auto)) {
    if (!missing(seasonal) && (seasonal.auto)){

      desc = decompose(y)$seasonal
      for (i in 1:tam.min) {
        if (desc[1] == desc[i + 1]){
          s.y = i
          break
        }
      }

      desc = decompose(x)$seasonal
      for (i in 1:tam.min) {
        if (desc[1] == desc[i + 1]){
          s.x = i
          break
        }
      }

      if (s.x != s.y){
        stop("Las series difieren en estacionalidad")
      }else{
        s = s.x
      }

      if (s != seasonal){
        stop("seasonal y seasonal.auto difieren(es recomendable que elija entre
             seasonal y seasonal.auto)")
      }
    }else{
      if (seasonal.auto) {

        desc = decompose(y)$seasonal
        for (i in 1:tam.min) {
          if (desc[1] == desc[i + 1]){
            s.y = i
            break
          }
        }

        desc = decompose(x)$seasonal
        for (i in 1:tam.min) {
          if (desc[1] == desc[i + 1]){
            s.x = i
            break
          }
        }

        if (s.x != s.y){
          stop("Las series difieren en estacionalidad")
        }else{
          s = s.x
        }

      }else{
        s = seasonal
      }
    }
  }else{
    s = 0
  }


  if (s == 0){

    if (!missing(lag.max)){

      if (tam.min <= lag.max){
        lag.max = tam.min - 1
      }

      ccf = stats::ccf(as.vector(x), as.vector(y), lag.max = lag.max, plot = FALSE)
      lag.value = ccf$lag
      ccf.value = ccf$acf
      data = data.frame(ccf.value = ccf$acf, lag.value = ccf$lag)


    }else{
      ccf = stats::ccf(as.vector(x), as.vector(y), plot = FALSE)
      lag.value = ccf$lag
      ccf.value = ccf$acf
      data = data.frame(ccf.value, lag.value)
    }

  }else{
    if (!missing(lag.max)) {
      if (tam.min <= lag.max){
        lag.max = tam.min - 1
      }

      Slag.max = 2*floor(lag.max/s) + 1
      if (Slag.max < 1) stop("Retardos insuficientes")


      ccf = stats::ccf(as.vector(x), as.vector(y), lag.max = lag.max, plot = FALSE)
      ccf.value = c(rep(0, Slag.max))
      lag.value = c(rep(0, Slag.max))
      ccf.value[floor(lag.max/s) + 1] = ccf$acf[lag.max + 1]

      cont = s - 1
      for (i in 1:(t = floor(Slag.max/2))) {
        ccf.value[t + 1 + i] = ccf$acf[(lag.max + 2):(2*lag.max + 1)][i*s]
        lag.value[t + 1 + i] = ccf$lag[(lag.max + 2):(2*lag.max + 1)][i*s]
        ccf.value[t + 1 - i] = ccf$acf[1:lag.max][lag.max - cont]
        lag.value[t + 1 - i] = ccf$lag[1:lag.max][lag.max - cont]
        cont = cont + s
      }
      data = data.frame(ccf.value, lag.value)

    }else{
      lag.max = tam.min - 1

      Slag.max = (2*floor(lag.max/s) + 1)
      ccf = stats::ccf(as.vector(x), as.vector(y), lag.max = lag.max, plot = FALSE)
      ccf.value = c(rep(0, Slag.max))
      lag.value = c(rep(0, Slag.max))
      ccf.value[floor(lag.max/s) + 1] = ccf$acf[lag.max + 1]

      cont = s - 1
      for (i in 1:(t = floor(Slag.max/2))) {
        ccf.value[t + 1 + i] = ccf$acf[(lag.max + 2):(2*lag.max + 1)][i*s]
        lag.value[t + 1 + i] = ccf$lag[(lag.max + 2):(2*lag.max + 1)][i*s]
        ccf.value[t + 1 - i] = ccf$acf[1:lag.max][lag.max - cont]
        lag.value[t + 1 - i] = ccf$lag[1:lag.max][lag.max - cont]
        cont = cont + s
      }

      data = data.frame(ccf.value, lag.value)

    }
  }

  ccf.out = list(ccf = as.vector(ccf.value), lag = as.vector(lag.value))

  if (summary){
    ccf.data = cbind(ccf.out$ccf, ccf.out$lag)
    colnames(ccf.data) = c("ccf", "lag")

    print(ccf.data, digits = 3)

    cat("Se basa en el supuesto que las series no están correlacionadas de forma cruzada y", "\n",
        "que una de las series es de ruido blanco.")
  }

  if(plot){

    IC = sqrt(tam.min)

    plot(ggplot(data = data, aes(y = ccf.value, x = lag.value)) + geom_bar(stat = "identity")+
           geom_hline(yintercept = 2/IC, linetype="dashed")+
           geom_hline(yintercept = -2/IC, linetype="dashed")+
           geom_hline(yintercept = 0)+
           labs(x = "LAG",y = "CCF"))

  }

  return(invisible(ccf.out))

}
