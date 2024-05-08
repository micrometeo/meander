library(modelbased)

synthesizer <- function(length=36000, step=0.1, amplitude=1, frequency=1/60, amplitude.2=0, frequency.2=1/30, noise.fraction=0.1) {
  noise <- rnorm(length, 0, amplitude*noise.fraction)
  t     <- (0:(length-1)) * step
  s     <- amplitude*sin(2*pi*frequency*t) + amplitude.2*sin(2*pi*frequency.2*t) + noise
  return(s)
}

zero.cross <- function(s,lag.max=3600, plot.show=FALSE) {
  c <- acf(s,lag.max=lag.max, plot=plot.show)$acf[,1,1]
  z <- zero_crossings(c)
  return(z)
}


# Read one hourly file
get.raw <- function(file.name) {

  # Get data, including the 1s resolution time stamp
  data <- readBin(file.name, "raw", 10e6)
  df <- as.data.frame(t(matrix(as.numeric(data), nrow=10)))
  two.compl <- function(x, y) { (x + 2^8*y) -> z; ifelse(z < 2^15, z, z - 2^16) }
  out <- data.frame(
    second = two.compl(df$V1, df$V2),
    u      = two.compl(df$V3, df$V4) / 100.0,
    v      = two.compl(df$V5, df$V6) / 100.0,
    w      = two.compl(df$V7, df$V8) / 100.0,
    t      = two.compl(df$V9, df$V10) / 100.0
  )
  out <- out[out$second < 5000,]
  out$vel <- sqrt(out$u^2 + out$v^2)
  out$dir <- 180/pi * atan2(-out$u, -out$v)
  out$dir[out$dir < 0] <- out$dir[out$dir < 0] + 360.0

  # Retime sequence
  n <- length(out$second)
  out$second <- 3600*(0:(n-1))/n

  # Yield result
  return(out)

}


get.1 <- function() {
  d <- get.raw("20160308.12R")
  return(list(
    time = as.POSIXct("2016-03-08 12:00:00", tz="UTC"),
    data = d
  ))
}


get.2 <- function() {
  d <- get.raw("20160308.18R")
  return(list(
    time = as.POSIXct("2016-03-08 18:00:00", tz="UTC"),
    data = d
  ))
}


get.3 <- function() {
  d <- get.raw("20160801.10R")
  return(list(
    time = as.POSIXct("2016-08-01 10:00:00", tz="UTC"),
    data = d
  ))
}


get.4 <- function() {
  d <- get.raw("20160801.23R")
  return(list(
    time = as.POSIXct("2016-08-01 23:00:00", tz="UTC"),
    data = d
  ))
}


# Filter data (assuming retiming was successful)
dFZ03 <- function(d, tau, delta.t=0.1) {

  # Calculate the filter constant
  rCos <- cos(-2*pi*delta.t/tau)
  rAlpha <- 2 - rCos - sqrt(rCos**2 - 4*rCos + 3)

  # Reserve workspace
  n <- length(d$second)
  forward.u <- numeric(n)
  forward.v <- numeric(n)
  forward.w <- numeric(n)
  filtered.u <- numeric(n)
  filtered.v <- numeric(n)
  filtered.w <- numeric(n)

  # Apply the forward-only McMillen filter using revised (dFZ03) coefficient
  forward.u[1] <- d$u[1]
  forward.v[1] <- d$v[1]
  forward.w[1] <- d$w[1]
  for(i in 2:n) {
    forward.u[i] <- rAlpha * forward.u[i-1] + (1 - rAlpha) * d$u[i]
    forward.v[i] <- rAlpha * forward.v[i-1] + (1 - rAlpha) * d$v[i]
    forward.w[i] <- rAlpha * forward.w[i-1] + (1 - rAlpha) * d$w[i]
  }

  # Filter backward 'a la McMillen, using dFZ03 revised coefficient again
  filtered.u[n] = forward.u[n]
  filtered.v[n] = forward.v[n]
  filtered.w[n] = forward.w[n]
  for(i in (n-1):1) {
    filtered.u[i] = rAlpha * filtered.u[i+1] + (1 - rAlpha) * forward.u[i]
    filtered.v[i] = rAlpha * filtered.v[i+1] + (1 - rAlpha) * forward.v[i]
    filtered.w[i] = rAlpha * filtered.w[i+1] + (1 - rAlpha) * forward.w[i]
  }

  e <- data.frame(
    second = d$second,
    u      = filtered.u,
    v      = filtered.v,
    w      = filtered.w
  )
  e$vel <- sqrt(e$u^2 + e$v^2)
  e$dir <- 180/pi * atan2(-e$u, -e$v)
  e$dir[e$dir < 0] <- e$dir[e$dir < 0] + 360.0

  return(e)

}


plot.data <- function(data.num) {

  if(data.num==1) {
    d <- get.1()
  }
  else if(data.num==2) {
    d <- get.2()
  }
  else if(data.num==3) {
    d <- get.3()
  }
  else if(data.num==4) {
    d <- get.4()
  }

  plot.postfix <- format(d$time, "%Y-%m-%d_%H.png", tz="UTC")
  d <- d$data
  e <- dFZ03(d, 60)

  png(file=paste("plots/Components_",plot.postfix,sep=""), width=16, height=12, units="in", res=96*4)
  par(mfrow=c(2,2))
  plot(d$second/60, d$u, type="l", xaxt="n", xlab="Time (min)", ylab="u (m/s)", col="lightgrey")
  axis(1, xaxp=c(0,60,6))
  lines(e$second/60, e$u, col="black", lwd=2)
  plot(d$second/60, d$v, type="l", xaxt="n", xlab="Time (min)", ylab="v (m/s)", col="lightgrey")
  axis(1, xaxp=c(0,60,6))
  lines(e$second/60, e$v, col="black", lwd=2)
  plot(d$second/60, d$w, type="l", xaxt="n", xlab="Time (min)", ylab="w (m/s)", col="lightgrey")
  axis(1, xaxp=c(0,60,6))
  lines(e$second/60, e$w, col="black", lwd=2)
  plot(d$second/60, d$dir, type="p", cex=0.2, xaxt="n", xlab="Time (min)", ylab="Dir (°)")
  axis(1, xaxp=c(0,60,6))
  par(mfrow=c(1,1))
  dev.off()

}


plot.trend <- function(data.num) {

  if(data.num==1) {
    d <- get.1()
  }
  else if(data.num==2) {
    d <- get.2()
  }
  else if(data.num==3) {
    d <- get.3()
  }
  else if(data.num==4) {
    d <- get.4()
  }

  plot.postfix <- format(d$time, "%Y-%m-%d_%H.png", tz="UTC")
  d <- d$data
  e.3600 <- dFZ03(d, 3600)
  e.0900 <- dFZ03(d, 0900)
  e.0450 <- dFZ03(d, 0450)

  png(file=paste("plots/Trend_",plot.postfix,sep=""), width=16, height=12, units="in", res=96*4)
  par(mfrow=c(2,2))
  plot(d$second/60, d$u, type="l", xaxt="n", xlab="Time (min)", ylab="u (m/s)", col="lightgrey")
  axis(1, xaxp=c(0,60,6))
  lines(e.3600$second/60, e.3600$u, col="black", lwd=2)
  lines(e.0900$second/60, e.0900$u, col="red", lwd=2)
  lines(e.0450$second/60, e.0450$u, col="blue", lwd=2)
  plot(d$second/60, d$v, type="l", xaxt="n", xlab="Time (min)", ylab="v (m/s)", col="lightgrey")
  axis(1, xaxp=c(0,60,6))
  lines(e.3600$second/60, e.3600$v, col="black", lwd=2)
  lines(e.0900$second/60, e.0900$v, col="red", lwd=2)
  lines(e.0450$second/60, e.0450$v, col="blue", lwd=2)
  plot(d$second/60, d$w, type="l", xaxt="n", xlab="Time (min)", ylab="w (m/s)", col="lightgrey")
  axis(1, xaxp=c(0,60,6))
  lines(e.3600$second/60, e.3600$w, col="black", lwd=2)
  lines(e.0900$second/60, e.0900$w, col="red", lwd=2)
  lines(e.0450$second/60, e.0450$w, col="blue", lwd=2)
  plot(d$second/60, d$dir, type="p", cex=0.2, xaxt="n", xlab="Time (min)", ylab="Dir (°)")
  axis(1, xaxp=c(0,60,6))
  par(mfrow=c(1,1))
  dev.off()

}


analyze.data <- function(data.num) {

  if(data.num==1) {
    d <- get.1()
  }
  else if(data.num==2) {
    d <- get.2()
  }
  else if(data.num==3) {
    d <- get.3()
  }
  else if(data.num==4) {
    d <- get.4()
  }

  # Remove trend
  f <- dFZ03(d$data, 900)
  d$data$u <- d$data$u - f$u
  d$data$v <- d$data$v - f$v
  d$data$w <- d$data$w - f$w

  # Filter data
  d <- d$data
  n <- 60
  u.25 <- numeric(n)
  u.50 <- numeric(n)
  u.75 <- numeric(n)
  v.25 <- numeric(n)
  v.50 <- numeric(n)
  v.75 <- numeric(n)
  w.25 <- numeric(n)
  w.50 <- numeric(n)
  w.75 <- numeric(n)
  for(i in 1:n) {
    e <- dFZ03(d, i)
    val.u <- fivenum(diff(zero.cross(e$u)))
    u.25[i] <- val.u[2]
    u.50[i] <- val.u[3]
    u.75[i] <- val.u[4]
    val.v <- fivenum(diff(zero.cross(e$v)))
    v.25[i] <- val.v[2]
    v.50[i] <- val.v[3]
    v.75[i] <- val.v[4]
    val.w <- fivenum(diff(zero.cross(e$w)))
    w.25[i] <- val.w[2]
    w.50[i] <- val.w[3]
    w.75[i] <- val.w[4]
    print(i)
  }

  f <- data.frame(
    tau = 1:n,
    u.25 = u.25,
    u.50 = u.50,
    u.75 = u.75,
    v.25 = v.25,
    v.50 = v.50,
    v.75 = v.75,
    w.25 = w.25,
    w.50 = w.50,
    w.75 = w.75
  )

}


plot.stats <- function(d.1, d.2, d.3, d.4, plot.name) {
  u.1 <- d.1$u.50
  v.1 <- d.1$v.50
  u.2 <- d.2$u.50
  v.2 <- d.2$v.50
  u.3 <- d.3$u.50
  v.3 <- d.3$v.50
  u.4 <- d.4$u.50
  v.4 <- d.4$v.50
  m.1 <-pmin(u.1,v.1,na.rm=TRUE)
  m.2 <-pmin(u.2,v.2,na.rm=TRUE)
  m.3 <-pmin(u.3,v.3,na.rm=TRUE)
  m.4 <-pmin(u.4,v.4,na.rm=TRUE)
  m <- max(c(m.1,m.2,m.3,m.4),na.rm=TRUE)/10
  png(file=plot.name, width=8, height=6, units="in", res=96*4)
  plot(1:60, m.1/10, type="l", ylim=c(0,m), xlab="Filter time width (s)", ylab="Median difference between consecutive zero crossings (s)", col="black")
  lines(1:60,m.2/10, col="red")
  lines(1:60,m.3/10, col="blue")
  lines(1:60,m.4/10, col="green")
  dev.off()
}
