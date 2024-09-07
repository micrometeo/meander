library(phonTools)

# SSB is a binary file format used by ARPA Lombardia to store ultrasonic
# anemometer raw data in a form allowing fast access in Fortran.
read.ssb <- function(file.name, dt = 0.1) {
  
  # Get file descriptor
  to.read <- file(file.name, "rb")
  
  # Get magic block
  magic <- readBin(to.read, "character", size=6)
  if(magic != "ssb_v0") {
    close(to.read)
    return(NULL)
  }
  # OK, this is likely to be an SSB V0 file, as expected
  
  # Skip reserved byte (by reading and not using it; this is actually the 2nd reserved,
  # the first was already consumed by the character read as terminating 0)
  res <- readBin(to.read, integer(), size=1, n=1, endian="little")

  # Get date
  iYear  <- readBin(to.read, integer(), size=2, n=1, endian="little", signed=FALSE)
  iMonth <- readBin(to.read, integer(), size=1, n=1)
  iDay   <- readBin(to.read, integer(), size=1, n=1)

  # Get number of data
  iTotalData <- readBin(to.read, integer(), size=4, n=1)
  iNominalData <- round(24*3600/dt, 0)
  if(iTotalData < iNominalData) {
    # Should never happen, data having been already censored against this condition
    close(to.read)
    return(NULL)
  }
  ivNumData  <- readBin(to.read, integer(), size=4, n=24)

  # Gather data
  ivTimeStamp <- readBin(to.read, integer(), size=2, n=iTotalData)
  ivU         <- readBin(to.read, integer(), size=2, n=iTotalData)
  ivV         <- readBin(to.read, integer(), size=2, n=iTotalData)
  ivW         <- readBin(to.read, integer(), size=2, n=iTotalData)
  ivT         <- readBin(to.read, integer(), size=2, n=iTotalData)

  daily.set <- data.frame(
    sec = ivTimeStamp[1:iNominalData],
    u   = ivU[1:iNominalData] / 100.0,
    v   = ivV[1:iNominalData] / 100.0,
    w   = ivW[1:iNominalData] / 100.0,
    t   = ivT[1:iNominalData] / 100.0
  )

  hourly.data <- iNominalData %/% 24
  hour.begin  <- (0:23)*hourly.data + 1
  hour.end    <- hour.begin + hourly.data - 1

  close(to.read)

  out <- list(
    dt         = dt,
    num.data   = ivNumData,
    hour.begin = hour.begin,
    hour.end   = hour.end,
    daily.set  = daily.set
  )
  
  return(out)
  
}


get.hour <- function(data.set, hour) {
  
  if(hour < 1 | hour > 24) return(NULL)
  
  # Get data for this hour
  u <- data.set$daily.set$u[data.set$hour.begin[hour]:data.set$hour.end[hour]]
  v <- data.set$daily.set$v[data.set$hour.begin[hour]:data.set$hour.end[hour]]
  w <- data.set$daily.set$w[data.set$hour.begin[hour]:data.set$hour.end[hour]]
  t <- data.set$daily.set$t[data.set$hour.begin[hour]:data.set$hour.end[hour]]
  
  # Yield result
  d <- data.frame(u, v, w, t)
  return(d)
  
}


compress.zeros <- function(idx, delta.idx=50) {
  
  n <- length(idx)
  
  if(n > 1) {
    
    # Label each index by the cluster ID it belongs to
    i.blk <- 1
    lbl <- numeric(n)
    lbl[1] <- i.blk
    for(i in 2:n) {
      if(idx[i]-idx[i-1] >= delta.idx) {
        i.blk <- i.blk + 1
      }
      lbl[i] <- i.blk
    }
    
    # Get the approximate center for each cluster
    cluster.min <- aggregate(1:n, by=list(lbl), FUN=min)$x
    cluster.max <- aggregate(1:n, by=list(lbl), FUN=max)$x
    idx <- idx[cluster.min + round((cluster.max-cluster.min)/2, 0)]
    
  }

  return(idx)
  
}


# Find the zero crossings (if any) of the ACF "ac" at time lags "tau" (s)
find.zero.crossings <- function(tau, ac) {
  
  n <- length(ac)
  
  # Indices at which passages occur
  ac.m <- ac[1:(n-1)]
  ac.p <- ac[2:n]
  z.c.idx <- compress.zeros(which(ac.m*ac.p <= 0), delta.idx=150)
  
  # Translate index passages to time
  if(length(z.c.idx) > 0) {
    z.c <- tau[z.c.idx]
  }
  else {
    z.c <- NULL
  }
  
  return(list(tm=z.c,idx=z.c.idx))
  
}


estimate.period <- function(tau, ac) {
  
  cross <- find.zero.crossings(tau, ac)
  tm  <- cross$tm
  idx <- cross$idx

  if(length(tm) >= 2) {
    
    # Some oscillation occurs, at a time scale significantly smaller
    # than averaging time (1h): it makes sense to proceed
    a <- abs(min(ac[idx[1]:idx[2]]))
    T <- 4*tm[1]
    p <- abs(log(a)/(2*tm[1]))

  }
  else {
    p <- NA
    T <- NA
  }
  
  return(list(p=p, T=T))
  
}


# Despite its name, this function performs not only trend removal,
# but also period estimation on all components
detrend <- function(d, dt=0.1) {
  
  # This is the assumed time stamp:
  time <- seq(from=0, to=3600-dt, by=dt)
  
  # Compute linear regressions
  l.u <- lm(d$u~time)
  l.v <- lm(d$v~time)
  l.w <- lm(d$w~time)
  l.t <- lm(d$t~time)
  
  # Compute detrended components by subtracting each data the linear model
  u.mult <- as.numeric(l.u$coefficients[2])
  u.off  <- as.numeric(l.u$coefficients[1])
  u.d <- d$u - u.mult*time - u.off
  v.mult <- as.numeric(l.v$coefficients[2])
  v.off  <- as.numeric(l.v$coefficients[1])
  v.d <- d$v - v.mult*time - v.off
  w.mult <- as.numeric(l.w$coefficients[2])
  w.off  <- as.numeric(l.w$coefficients[1])
  w.d <- d$w - w.mult*time - w.off
  t.mult <- as.numeric(l.t$coefficients[2])
  t.off  <- as.numeric(l.t$coefficients[1])
  t.d <- d$t - t.mult*time - t.off
  
  # Compute autocorrelations
  m <- length(u.d) %/% 2  # Maximum lag not too large
  tau <- time[1:m]
  u.a <- fastacf(u.d, lag.max=m, show=FALSE, correct=TRUE)$acf
  v.a <- fastacf(v.d, lag.max=m, show=FALSE, correct=TRUE)$acf
  w.a <- fastacf(w.d, lag.max=m, show=FALSE, correct=TRUE)$acf
  t.a <- fastacf(t.d, lag.max=m, show=FALSE, correct=TRUE)$acf
  
  # Approximately fit autocorrelations to theoretical model
  per.u <- estimate.period(tau, u.a)
  per.v <- estimate.period(tau, v.a)
  per.w <- estimate.period(tau, w.a)
  per.t <- estimate.period(tau, t.a)

  # Compose return value
  out <- list(
    u.mult = u.mult,
    v.mult = v.mult,
    w.mult = w.mult,
    t.mult = t.mult,
    u.off  = u.off,
    v.off  = v.off,
    w.off  = w.off,
    t.off  = t.off,
    u      = u.d,
    v      = v.d,
    w      = w.d,
    t      = t.d,
    u.a    = u.a,
    v.a    = v.a,
    w.a    = w.a,
    t.a    = t.a,
    per.u  = per.u,
    per.v  = per.v,
    per.w  = per.w,
    per.t  = per.t
  )
  return(out)
  
}


# Data processing
process.period <- function(station.altitude=163.0, zr=10) {
  
  # Retrieve file names from data directory, and process them in turn
  files <- get.file.names()
  n.files <- length(files)
  n.hours <- n.files * 24
  time.stamp <- rep(as.POSIXct("2000-01-01 00:00:00", tz="UTC"), times=n.hours)
  file.name <- rep("", times=n.hours)
  p.u <- numeric(n.hours)
  p.v <- numeric(n.hours)
  p.w <- numeric(n.hours)
  p.t <- numeric(n.hours)
  T.u <- numeric(n.hours)
  T.v <- numeric(n.hours)
  T.w <- numeric(n.hours)
  T.t <- numeric(n.hours)
  meander.level <- numeric(n.hours)
  vel <- numeric(n.hours)
  vel.scalar <- numeric(n.hours)
  S   <- numeric(n.hours)
  us  <- numeric(n.hours)
  H0  <- numeric(n.hours)
  dir <- numeric(n.hours)
  cur.hour <- 0
  for(f in files) {
    
    # Get raw data
    daily.set <- read.ssb(f, dt=0.1)
    print(f)
    
    # Compute this day's time stamp, in string form
    f.len <- nchar(f)
    date.stamp <- substr(f, f.len-13, f.len-4)

    # Main loop: process hours
    for(hour in 1:24) {
      
      cur.hour <- cur.hour + 1
      
      # Get current hour for processing, and inform user we're here
      data <- get.hour(daily.set, hour)

      # Remove linear trend, compute autocorrelation, and fit model
      data.d <- detrend(data)
      
      # Save data
      p.u[cur.hour] <- data.d$per.u$p
      p.v[cur.hour] <- data.d$per.v$p
      p.w[cur.hour] <- data.d$per.w$p
      p.t[cur.hour] <- data.d$per.t$p
      T.u[cur.hour] <- data.d$per.u$T
      T.v[cur.hour] <- data.d$per.v$T
      T.w[cur.hour] <- data.d$per.w$T
      T.t[cur.hour] <- data.d$per.t$T

      # Eddy covariance
      
      u <- data$u
      v <- data$v
      w <- data$w
      t <- data$t
      avg.u <- mean(u)
      avg.v <- mean(v)
      avg.w <- mean(w)
      avg.t <- mean(t)
      u.d <- data.d$u
      v.d <- data.d$v
      w.d <- data.d$w
      t.d <- data.d$t
      uu <- cov(u.d,u.d)
      vv <- cov(v.d,v.d)
      ww <- cov(w.d,w.d)
      uv <- cov(u.d,v.d)
      uw <- cov(u.d,w.d)
      vw <- cov(v.d,w.d)
      ut <- cov(u.d,t.d)
      vt <- cov(v.d,t.d)
      wt <- cov(w.d,t.d)
      vel[cur.hour] <- sqrt(avg.u^2 + avg.v^2)
      vel.scalar[cur.hour] <- mean(sqrt(u^2 + v^2))
      time.stamp[cur.hour] <- as.POSIXct(paste(date.stamp,sprintf("%2.2d",cur.hour),sep=" "),tz="UTC")
      dir[cur.hour] <- atan2(-avg.u, -avg.v) * 180/pi
      file.name[cur.hour] <- f
      
      # Form relevant matrices
      wind.avg   <- matrix(data=c(avg.u, avg.v, avg.w), nrow=3, ncol=1, byrow=TRUE)
      wind.cov   <- matrix(data=c(uu, uv, uw, uv, vv, vw, uw, vw, ww), nrow=3, ncol=3, byrow=TRUE);
      wind.t.cov <- matrix(data=c(ut, vt, wt), nrow=3, ncol=1, byrow=TRUE);
      
      # Build and apply first rotation matrix
      denom <- sqrt(avg.u^2 + avg.v^2)
      if(!is.null(denom) && !is.na(denom) && denom > 1.e-3) {
        sin.a <- avg.v / denom
        cos.a <- avg.u / denom
      }
      else {
        sin.a <- 0.0
        cos.a <- 1.0
      }
      R01 <- matrix(c(cos.a, sin.a, 0, -sin.a, cos.a, 0, 0, 0, 1), nrow=3, ncol=3, byrow=TRUE)
      wind.avg <- R01 %*% wind.avg
      wind.cov <- R01 %*% wind.cov %*% t(R01)
      wind.t.cov <- R01 %*% wind.t.cov
      
      # Build and apply second rotation matrix
      denom <- sqrt(wind.avg[1]^2 + wind.avg[3]^2)
      if(!is.null(denom) && !is.na(denom) && denom > 1.e-3) {
        sin.b <- wind.avg[3] / denom
        cos.b <- wind.avg[1] / denom
      }
      else {
        sin.b <- 0.0
        cos.b <- 1.0
      }
      R12 <- matrix(c(cos.b, 0, sin.b, 0, 1, 0, -sin.b, 0, cos.b), nrow=3, ncol=3, byrow=TRUE)
      wind.avg <- R12 %*% wind.avg
      wind.cov <- R12 %*% wind.cov %*% t(R12)
      wind.t.cov <- R12 %*% wind.t.cov
      
      # Compute u* and H0
      Ta     <- avg.t + 273.15
      Pa     <- 1013.0 * exp(-0.0342*station.altitude/Ta)
      RhoCp  <- 350.25*Pa/Ta
      u.star <- (wind.cov[1,3]^2 + wind.cov[2,3]^2)^0.25
      H0[cur.hour] <- RhoCp * wind.t.cov[3]
      us[cur.hour] <- u.star
      
      # Compute Obukhov length and stability parameter
      L           <- -Ta / (0.4*9.81) * u.star^3/wind.t.cov[3]
      S[cur.hour] <- zr/L
      
    }
    
  }
  
  # Calculate meandering factors
  m.u <- 2*pi / (p.u*T.u)
  m.v <- 2*pi / (p.v*T.v)
  m.w <- 2*pi / (p.w*T.w)
  m.t <- 2*pi / (p.t*T.t)
  m   <- pmin(m.u, m.v, m.t)
  
  # Classify the hours in meandering or not depending on m.x (x=u,v,t)
  meander <- (m >= 1) & (vel.scalar <= 1.5)
  meander[is.na(meander)] <- FALSE
  meandering <- ifelse(meander, "meandering", "non-meandering")
  
  # Classify the hours in slow or fast, and stable or not
  slow.wind <- ifelse(vel.scalar <= 1, "slow wind (< 1m/s)", "fast wind")
  stable    <- ifelse(S > 0.1, "stable", ifelse(S < -0.1, "unstable", "neutral"))
  
  out <- data.frame(
    time.stamp,
    file.name,
    p.u,
    p.v,
    p.w,
    p.t,
    T.u,
    T.v,
    T.w,
    T.t,
    m.u,
    m.v,
    m.w,
    m.t,
    m,
    vel,
    dir,
    vel.scalar,
    S,
    us,
    H0,
    meandering,
    slow.wind,
    stable
  )
  return(out)
    
}


# 'pp' is the output of "process.period"
filter.ssb <- function(pp, idx.ssb, dt=0.1) {
  
  # Get the SSB file
  ssb.file <- pp$file.name[idx.ssb]
  data.set <- read.ssb(ssb.file, dt)
  
  # Deduce hour, and get it
  base.idx <- (idx.ssb %/% 24) * 24
  hour <- idx.ssb - base.idx
  hour.set <- get.hour(data.set, hour)
  
  # Calculate the filter constant
  tau <- min(c(pp$T.u[idx.ssb], pp$T.v[idx.ssb]))
  delta.t <- data.set$dt
  rCos <- cos(-2*pi*delta.t/tau)
  rAlpha <- 2 - rCos - sqrt(rCos**2 - 4*rCos + 3)
  
  # Reserve workspace
  n <- length(hour.set$u)
  forward.u <- numeric(n)
  forward.v <- numeric(n)
  filtered.u <- numeric(n)
  filtered.v <- numeric(n)
  
  # Apply the forward-only McMillen filter using revised (dFZ03) coefficient
  u <- hour.set$u
  v <- hour.set$v
  forward.u[1] <- u[1]
  forward.v[1] <- v[1]
  for(i in 2:n) {
    forward.u[i] <- rAlpha * forward.u[i-1] + (1 - rAlpha) * u[i]
    forward.v[i] <- rAlpha * forward.v[i-1] + (1 - rAlpha) * v[i]
  }
  
  # Filter backward 'a la McMillen, using dFZ03 revised coefficient again
  filtered.u[n] = forward.u[n]
  filtered.v[n] = forward.v[n]
  for(i in (n-1):1) {
    filtered.u[i] = rAlpha * filtered.u[i+1] + (1 - rAlpha) * forward.u[i]
    filtered.v[i] = rAlpha * filtered.v[i+1] + (1 - rAlpha) * forward.v[i]
  }
  
  second <- seq(from=0, to=3600-dt, by=dt)
  u.f    <- filtered.u
  v.f    <- filtered.v
  
  # Compose virtual trajectories
  x   <- cumsum(u)*dt
  y   <- cumsum(v)*dt
  x.f <- cumsum(u.f)*dt
  y.f <- cumsum(v.f)*dt
  
  # Calculate limiting square (for plotting)
  if(sum(!is.na(x.f)) > 0) {
    mn <- min(c(x,y,x.f,y.f))
    mx <- max(c(x,y,x.f,y.f))
  }
  else {
    mn <- min(c(x,y))
    mx <- max(c(x,y))
  }
  
  # Yield result
  out <- list(
    sw  = mn,
    ne  = mx,
    x   = x,
    y   = y,
    x.f = x.f,
    y.f = y.f
  )
  return(out)
  
}
