#' @export tf.acf
tf.acf = function(x, lag.max, seasonal, seasonal.auto = FALSE, plot = TRUE,
                summary = TRUE, type.test = "Ljung-Box")
{
  x = na.omit(x)
  if(!is.ts(x)) stop("x debe ser una serie de tiempo")
  if (is.mts(x)) stop("Esta funcion es solo para series de tiempo univariantes")
  if (!missing(lag.max) && lag.max < 0) stop("lag.max debe ser > 0")
  if (!is.logical(seasonal.auto)) stop("seasonal debe ser un valor logico")
  if (!is.logical(plot)) stop("plot debe ser ")
  if (!missing(seasonal) && !is.numeric(seasonal)) stop("seasonal debe ser un valor numerico")
  if (!missing(seasonal) && length(seasonal) > 1) stop("seasonal de ver un numero positivo")
  if (!missing(seasonal) && seasonal < 0) stop("seasonal debe ser un numero positivo")

  if (!missing(seasonal) || (seasonal.auto == TRUE)) {
    if (!missing(seasonal) && (seasonal.auto == TRUE)){

      desc = decompose(x)$seasonal
      for (i in 1:length(desc)) {
        if (desc[1] == desc[i + 1]){
          s = i
          break
        }
      }

      if (s != seasonal){
        stop("seasonal y seasonal.auto difieren(es recomendable que elija entre
             seasonal y seasonal.auto)")
      }
    }else{
      if (seasonal.auto) {

        desc = decompose(x)$seasonal
        for (i in 1:length(desc)) {
          if (desc[1] == desc[i + 1]){
            s = i
            break
          }
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
      acf.value = stats::acf(x, lag.max = (lag.max + 1), plot = FALSE)$acf[2:(lag.max + 1)]
      lag.value = 1:length(acf.value)

      statistic = data.frame(rep(0,lag.max))
      df = data.frame(rep(0,lag.max))
      p.value = data.frame(rep(0,lag.max))

      if (type.test == "Ljung-Box"){
        for (i in 1:lag.max) {

          box = stats::Box.test(x = x, lag = i, type = "Ljung-Box")
          statistic[i,1] = box$statistic
          df[i,1] = box$parameter
          p.value[i,1] = box$p.value

        }
      }else{

        for (i in 1:lag.max) {
          box = stats::Box.test(x = x, lag = i, type = "Box-Pierce")
          statistic[i,1] = box$statistic
          df[i,1] = box$parameter
          p.value[i,1] = box$p.value
        }
      }

    }else{
      acf.value = stats::acf(x,plot = FALSE)$acf
      acf.value = acf.value[2:length(acf.value)]
      lag.value = 1:length(acf.value)
      lag.max = length(acf.value)

      statistic = data.frame(rep(0,lag.max))
      df = data.frame(rep(0,lag.max))
      p.value = data.frame(rep(0,lag.max))

      if (type.test == "Ljung-Box"){
        for (i in 1:lag.max) {
          box = stats::Box.test(x = x, lag = i, type = "Ljung-Box")
          statistic[i,1] = box$statistic
          df[i,1] = box$parameter
          p.value[i,1] = box$p.value
        }

      }else{
        for (i in 1:lag.max) {
          box = stats::Box.test(x = x, lag = i, type = "Box-Pierce")
          statistic[i,1] = box$statistic
          df[i,1] = box$parameter
          p.value[i,1] = box$p.value
        }
      }
    }

  }else{
    if (!missing(lag.max)) {

      Slag.max = floor(lag.max/s)
      if (Slag.max < 1){
        stop("Retardos insuficientes")
      }

      acf.value = c(rep(0, Slag.max))
      lag.value = c(rep(0, Slag.max))
      y = stats::acf(x, lag.max = (lag.max + 1), plot = FALSE)$acf[2:(lag.max + 1)]

      for (i in 1:Slag.max) {
        acf.value[i] = y[i*s]
        lag.value[i] = i*s
      }

      statistic = data.frame(rep(0,Slag.max))
      df = data.frame(rep(0,Slag.max))
      p.value = data.frame(rep(0,Slag.max))

      if (type.test == "Ljung-Box"){
        for (i in 1:Slag.max) {
          box = stats::Box.test(x = x, lag = i, type = "Ljung-Box")
          statistic[i,1] = box$statistic
          df[i,1] = box$parameter
          p.value[i,1] = box$p.value
        }
      }else{
        for (i in 1:Slag.max) {
          box = stats::Box.test(x = x, lag = i, type = "Box-Pierce")
          statistic[i,1] = box$statistic
          df[i,1] = box$parameter
          p.value[i,1] = box$p.value
        }
      }

    }else{
      Slag.max = floor((length(x)/s) - 1)
      y = stats::acf(x, lag.max = length(x), plot = FALSE)$acf[2:length(x)]
      acf.value = c(rep(0, Slag.max))
      lag.value = c(rep(0, Slag.max))

      cont = 1
      while (length(x) > cont*s) {
        lag.value[cont] = cont*s
        acf.value[cont] = y[cont*s]
        cont = cont + 1
      }

      statistic = data.frame(rep(0,Slag.max))
      df = data.frame(rep(0,Slag.max))
      p.value = data.frame(rep(0,Slag.max))

      if (type.test == "Ljung-Box"){
        for (i in 1:Slag.max) {
          box = stats::Box.test(x = x, lag = i, type = "Ljung-Box")
          statistic[i,1] = box$statistic
          df[i,1] = box$parameter
          p.value[i,1] = box$p.value
        }
      }else{
        for (i in 1:Slag.max) {
          box = stats::Box.test(x = x, lag = i, type = "Box-Pierce")
          statistic[i,1] = box$statistic
          df[i,1] = box$parameter
          p.value[i,1] = box$p.value
        }
      }
    }
  }

  if (length(lag.value) != length(as.ts(df))){
    if (length(lag.value) > length(as.ts(df))){
      lag.value = lag.value[1:(length(lag.value) - 1)]
      acf.value = acf.value[1:(length(acf.value) - 1)]
    }else{
      df = df[1:(length(as.ts(df)) - 1), 1]
      p.value = p.value[1:(length(as.ts(p.value)) - 1), 1]
      statistic = statistic[1:(length(as.ts(statistic)) - 1), 1]
    }
  }

  data = data.frame(lag.value, acf.value, df, statistic, p.value)

  if (type.test == "Ljung-Box"){

    colnames(data) = c("lag", "acf", "(L-B) df","(L-B) statistic", "(L-B) p.value")
  }else{

    colnames(data) = c("lag","acf", "(B-P) df", "(B-P) statistic", "(B-P) p.value")
  }

  if (summary){
    print(data)
  }

  if(plot){

    data.plot = data.frame(lag.value, acf.value)

    plot(ggplot(data = data.plot, aes(y = acf.value, x = lag.value)) + geom_bar(stat = "identity")+
           geom_hline(yintercept = 2/sqrt(length(x)), linetype="dashed")+
           geom_hline(yintercept = -2/sqrt(length(x)), linetype="dashed")+
           geom_hline(yintercept = 0)+
           labs(x = "LAG",y = "ACF"))

  }


  if (type.test == "Box-Pierce"){

    Box.Pierce = data.frame(statistic, df, p.value)
    colnames(Box.Pierce) = c("statistic", "df", "p.value")

    if (s == 0){
      return(invisible(list(acf = acf.value, lag = lag.value, Box.Pierce = Box.Pierce)))
    }else{
      return(invisible(list(acf = acf.value, lag = lag.value, seasonal = s, Box.Pierce = Box.Pierce)))
    }

  }else{

    Ljung.Box = data.frame(statistic, df, p.value)
    colnames(Ljung.Box) = c("statistic", "df", "p.value")

    if (s == 0){
      return(invisible(list(acf = acf.value, lag = lag.value, Ljung.Box = Ljung.Box)))
    }else{
      return(invisible(list(aacf = acf.value, lag = lag.value, seasonal = s, Ljung.Box = Ljung.Box)))
    }

  }

}

