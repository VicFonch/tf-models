#' @export tf.decompose
tf.decompose = function (x, type = c("additive", "multiplicative"))
{
  type = match.arg(type)
  n = length(x)
  f = frequency(x)
  s = start(x)
  if (f <= 1 || length(na.omit(x)) < 2*f) stop("La serie de tiempo no tiene frequencia mayor o igual 2")

  filter = if (f%%2 != 0) {
    c(0.5, rep(1, f - 1), 0.5)/f
  }else rep(1, f)/f

  trend = filter(x, filter)
  season = if (type == "additive"){
    x - trend
  }else x/trend

  periods = n%/%f
  index = seq.int(1, n, by = f) - 1
  figure = f

  for (i in 1:f) figure[i] = mean(season[index + i], na.rm = TRUE)

  figure = if (type == "additive"){
    figure - mean(figure)
  }else figure/mean(figure)

  seasonal = ts(rep(figure, periods + 1)[seq_len(n)], start = s,
                 frequency = f)

  t = if (type == "additive"){
    x - seasonal
  }else x/seasonal

  t = c(rep(t[1], 3), t, rep(t[n], 2))
  filter = c(1:3, 2:1)/9
  trend = ts(filter(t, filter)[3:(n + 2)], start = s, frequency = f)

  random = if (type == "additive"){
    x - seasonal - trend
  } else x/(seasonal*trend)

  structure(list(x = x, seasonal = seasonal, trend = trend,
            random = random, figure = figure, type = type),
            class = "decomposed.ts")
}

