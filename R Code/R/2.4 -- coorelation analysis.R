# Header
{
  # Load required packages
  {
    require(rstudioapi)
    require(readxl)
    require(writexl)
  }
  # Clear workspace and set working directory
  {
    rm(list=objects())
    this_fp <- function(){return(rstudioapi::getActiveDocumentContext()$path)}
    this_dir <- function(){return(dirname(this_fp()))}
    setwd(this_dir())
    set.seed(333)
  }
}
# Read in data
{
  # Raw abundance data
  Beluga_RAW_Nt <- readxl::read_xlsx("../Input Data/Beluga.xlsx")  
  Beluga_RAW_Nt <- data.frame(
    Year = Beluga_RAW_Nt$YR,
    Nt = Beluga_RAW_Nt$Nhat,
    SD = Beluga_RAW_Nt$SD_Nhat
  )
  # Model results from unfiltered abundance data
  dir <- "20250131"
  
  Chains <- list()
  for (i in 1:4){
    load(paste0("../Output Data/Vital Rates/", dir, "/chain", i, ".Rdata"))
    Chains[[i]] <- MCMC_out
  }
  # Indexes of availability
  prey_ff <- readxl::read_xlsx("../Output Data/Salmon_Availability.xlsx")
}
# Helper functions
{
  # Get posterior samples for 
    # per adult rate of reproduction
    # net survival rate (across the entire population) 
  get_posterior_samples <- function(param, ylim, Chains, nsamp, burnin){
    samples <- matrix(NA, nrow = nsamp, ncol = 14)
    for (i in 1:nsamp){
      chain <- sample(1:4, 1)
      sample <- sample((burnin+1):nrow(Chains[[chain]]$Pr), 1)
      if (param == "Pr"){
        samples[i,] <- Chains[[chain]]$Pr[sample,]
      }
      if (param == "Ps"){
        samples[i,] <- (
          Chains[[chain]]$Ps[sample,,1]*Chains[[chain]]$Na[sample,,1] +
          Chains[[chain]]$Ps[sample,,2]*Chains[[chain]]$Na[sample,,2] +
          Chains[[chain]]$Ps[sample,,3]*Chains[[chain]]$Na[sample,,3]
        )/(Chains[[chain]]$Na[sample,,1] + Chains[[chain]]$Na[sample,,2] + Chains[[chain]]$Na[sample,,3])
      }
    }
    samples <- samples[,(2005:2018) >= ylim[1] & (2005:2018) <= ylim[2]]
    return(samples)
  }
  # Get index of prey availability
  get_index <- function(prey_ff, ylim, species, river){
    print(
      prey_ff$Year[
        prey_ff$Year >= ylim[1] &
          prey_ff$Year <= ylim[2] &
          prey_ff$Species==species &
          prey_ff$River==river
      ]
    )
    return(
      prey_ff$Biomass_MT[
        prey_ff$Year >= ylim[1] &
          prey_ff$Year <= ylim[2] &
          prey_ff$Species==species &
          prey_ff$River==river
      ]
    )
  }
  # Get inverse cdf
  get_inverse_cdf <- function(param, index){
    pdf <- numeric(nrow(param))
    for (i in 1:nrow(param)) {
      pdf[i] <- cor(param[i,1:length(index)], index)
    }
    x <- seq(-1, 1, length =  100)
    y <- numeric(100)
    for (i in 1:length(x)){
      y[i] <- sum(pdf > x[i])/length(pdf)
    }
    return(data.frame(
      x = x,
      y = y
    ))
  }
  # Plot cdfs
  plot_inverse_cdfs <- function(cdf1, cdf2, cdf3, cdf4, names, par){
    par(mfrow=c(2,2))
    plot(cdf1$x, cdf1$y, ty = "l", main = names[1], xlab = NA, ylab = "probability")
    for (m in c(0, 0.2, 0.4, 0.6, 0.8, 1)){
      lines(
        c(-2, 2),
        c(m, m),
        col = "lightgrey"
      )
    }
    lines(c(0,0), c(0,cdf1$y[abs(cdf1$x)==min(abs(cdf1$x))]))
    if (par == "Pr"){
      lines(c(-20,cdf1$x[cdf1$y==cdf1$y[abs(cdf1$x)==min(abs(cdf1$x))]]), c(cdf1$y[abs(cdf1$x)==min(abs(cdf1$x))],cdf1$y[abs(cdf1$x)==min(abs(cdf1$x))]))
    }
    if (par == "Ps"){
      lines(c(-20,0), c(cdf1$y[abs(cdf1$x)==min(abs(cdf1$x))],cdf1$y[abs(cdf1$x)==min(abs(cdf1$x))]))
    }
    lines(cdf1$x, cdf1$y, col = "blue", lwd = 2)
    text(0.75, 0.9, labels = paste("p = ",  sprintf("%.2f", cdf1$y[abs(cdf1$x) == min(abs(cdf1$x))])), cex = 1.1) 
    plot(cdf2$x, cdf2$y, ty = "l", main = names[2], xlab = NA, ylab = "probability")
    for (m in c(0, 0.2, 0.4, 0.6, 0.8, 1)){
      lines(
        c(-2, 2),
        c(m, m),
        col = "lightgrey"
      )
    }
    lines(c(0,0), c(0,cdf2$y[abs(cdf2$x)==min(abs(cdf2$x))]))
    if (par == "Pr"){
      lines(c(-20,cdf2$x[cdf2$y==cdf2$y[abs(cdf2$x)==min(abs(cdf2$x))]]), c(cdf2$y[abs(cdf2$x)==min(abs(cdf2$x))],cdf2$y[abs(cdf2$x)==min(abs(cdf2$x))]))
    }
    if (par == "Ps"){
      lines(c(-20,0), c(cdf2$y[abs(cdf2$x)==min(abs(cdf2$x))],cdf2$y[abs(cdf2$x)==min(abs(cdf2$x))]))
    }
    lines(cdf2$x, cdf2$y, col = "blue", lwd = 2)
    text(0.75, 0.9, labels = paste("p = ",  sprintf("%.2f", cdf2$y[abs(cdf2$x) == min(abs(cdf2$x))])), cex = 1.1)
    plot(cdf3$x, cdf3$y, ty = "l", main = names[3], xlab = NA, ylab = "probability")
    for (m in c(0, 0.2, 0.4, 0.6, 0.8, 1)){
      lines(
        c(-2, 2),
        c(m, m),
        col = "lightgrey"
      )
    }
    lines(c(0,0), c(0,cdf3$y[abs(cdf3$x)==min(abs(cdf3$x))]))
    lines(c(-20,cdf3$x[cdf3$y==cdf3$y[abs(cdf3$x)==min(abs(cdf3$x))]]), c(cdf3$y[abs(cdf3$x)==min(abs(cdf3$x))],cdf3$y[abs(cdf3$x)==min(abs(cdf3$x))]))
    lines(cdf3$x, cdf3$y, col = "blue", lwd = 2)
    text(0.75, 0.9, labels = paste("p = ",  sprintf("%.2f", cdf3$y[abs(cdf3$x) == min(abs(cdf3$x))])), cex = 1.1)
    
    plot(cdf4$x, cdf4$y, ty = "l", main = names[4], xlab = NA, ylab = "probability")
    for (m in c(0, 0.2, 0.4, 0.6, 0.8, 1)){
      lines(
        c(-2, 2),
        c(m, m),
        col = "lightgrey"
      )
    }
    lines(c(0,0), c(0,cdf4$y[abs(cdf4$x)==min(abs(cdf4$x))]))
    if (par == "Pr"){
      lines(c(-20,cdf4$x[cdf4$y==cdf4$y[abs(cdf4$x)==min(abs(cdf4$x))]]), c(cdf4$y[abs(cdf4$x)==min(abs(cdf4$x))],cdf4$y[abs(cdf4$x)==min(abs(cdf4$x))]))
    }
    if (par =="Ps"){
      lines(c(-20,0), c(cdf4$y[abs(cdf4$x)==min(abs(cdf4$x))],cdf4$y[abs(cdf4$x)==min(abs(cdf4$x))]))
    }
    lines(cdf4$x, cdf4$y, col = "blue", lwd = 2)
    text(0.75, 0.9, labels = paste("p = ",  sprintf("%.2f", cdf4$y[abs(cdf4$x) == min(abs(cdf4$x))])), cex = 1.1)
  }
  # Time series plot
  plot_time_series <- function(years, prey_indexes, post_samples, nsamp = 50, names, param){
    par(mfrow = c(2,2), mar = c(4,4,4,4))
    for (ip in 1:length(prey_indexes)){
      # Get data for given panel
      year <- years[[ip]]
      prey_index <- prey_indexes[[ip]]
      samples <- post_samples[[ip]][1:nsamp,]
      # Rescale data for plotting
      prey_index_scaled <- (prey_index - min(prey_index))/(max(prey_index)-min(prey_index))
      samples_scaled <- (samples - min(samples))/(max(samples) - min(samples))
      # Plotting
      plot(year, prey_index_scaled, yaxt = "n", xaxt = "n", ty = "l", lwd = 2, main = "", xlim = c(2004, 2017), xlab = NA, ylab = NA)
      title(main = names[ip], line = 3)
      for (m in c(0, 0.2, 0.4, 0.6, 0.8, 1)){
        lines(
          c(0, 5000),
          c(m, m),
          col = "lightgrey"
        )
      }
      for (is in 1:nsamp){
        lines(year, samples_scaled[is,], col = rgb(1, 0, 0, alpha = 0.5))
      }
      lines(year, prey_index_scaled, col = "blue", lwd = 2)
      points(year, prey_index_scaled, pch = 19)
      if (param == "Pr"){
        axis(side = 1, labels = as.character(year - 1), at = year, col.axis = "blue")
        axis(side = 3, labels = as.character(year), at = year, col.axis = "red")
      }
      if (param == "Ps"){
        axis(side = 1, labels = as.character(year), at = year)
      }
      y2_lab <- round(c(0, 0.2, 0.4, 0.6, 0.8, 1) * (max(prey_index) - min(prey_index)) + min(prey_index),0)
      axis(side = 2, labels = y2_lab, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), col.axis = "blue")
      mtext("Biomass (metric tons)", side = 2, line = 2.1)
      y2_lab <- round(c(0, 0.2, 0.4, 0.6, 0.8, 1) * (max(samples) - min(samples)) + min(samples),2)
      axis(side = 4, labels = y2_lab, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), col.axis = "red")
      if (param == "Pr"){
        mtext("Per-Adult Reproduction Rate", side = 4, line = 2.1)
      }
      if (param == "Ps"){
        mtext("Rate of Survival to Next Year", side = 4, line = 2.1)
      }
    }
  }
  
  # Time series plot
  plot_prey_index_survival <- function(years, prey_indexes, post_samples, names){
    par(mfrow = c(2,2), mar = c(4,4,4,4))
    for (ip in 1:length(prey_indexes)){
      year <- years[[ip]]
      prey_index <- prey_indexes[[ip]]
      prey_index_scaled <- (prey_index - min(prey_index))/(max(prey_index)-min(prey_index))
      plot(year, prey_index_scaled, yaxt = "n", xaxt = "n", ty = "l", lwd = 2, main = "", xlim = c(2005, 2017), xlab = NA, ylab = NA)
      title(main = names[ip])
      for (m in c(0, 0.2, 0.4, 0.6, 0.8, 1)){
        lines(
          c(0, 5000),
          c(m, m),
          col = "lightgrey"
        )
      }
      lines(year, prey_index_scaled, col = "purple", lwd = 2)
      points(year, prey_index_scaled, pch = 19)
      axis(side = 1, labels = as.character(2005:2017), at = 2005:2017)
      y2_lab <- round(c(0, 0.2, 0.4, 0.6, 0.8, 1) * (max(prey_index) - min(prey_index)) + min(prey_index),0)
      axis(side = 2, labels = y2_lab, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0))
      mtext("Biomass (metric tons)", side = 2, line = 2.1)
      
    }
  }
  
}
## Reproduction
# Lagged prey indexes
Chinook_ND <- get_index(prey_ff, ylim = c(2004, 2016), species = "Chinook", river = "Susitna")
Chinook_CD <- get_index(prey_ff, ylim = c(2004, 2016), species = "Chinook", river = "Kenai")
Sockeye_ND <- get_index(prey_ff, ylim = c(2006, 2015), species = "Sockeye", river = "Susitna")
Sockeye_CD <- get_index(prey_ff, ylim = c(2004, 2015), species = "Sockeye", river = "Kenai") + get_index(prey_ff, ylim = c(2004, 2015), species = "Sockeye", river = "Kasilof")
# Posterior reproductive rates
Pr1 <- get_posterior_samples("Pr", c(2005, 2017), Chains, 10000, 100)
Pr2 <- get_posterior_samples("Pr", c(2005, 2017), Chains, 10000, 100)
Pr3 <- get_posterior_samples("Pr", c(2007, 2016), Chains, 10000, 100)
Pr4 <- get_posterior_samples("Pr", c(2005, 2016), Chains, 10000, 100)
# Inverse cdfs
cdf1 <- get_inverse_cdf(Pr1, Chinook_ND)
cdf2 <- get_inverse_cdf(Pr2, Chinook_CD)
cdf3 <- get_inverse_cdf(Pr3, Sockeye_ND)
cdf4 <- get_inverse_cdf(Pr4, Sockeye_CD)
# Names for plots
names <- c(
  "Northern District Chinook",
  "Central District Chinook",
  "Northern District Sockeye",
  "Central District Sockeye"
)
# Plot cdfs
plot_inverse_cdfs(cdf1, cdf2, cdf3, cdf4, names, par = "Pr") # Graphic for publication
## Survival
# Lagged prey indexes
Chinook_ND <- get_index(prey_ff, ylim = c(2005, 2017), species = "Chinook", river = "Susitna")
Chinook_CD <- get_index(prey_ff, ylim = c(2005, 2017), species = "Chinook", river = "Kenai")
Sockeye_ND <- get_index(prey_ff, ylim = c(2006, 2015), species = "Sockeye", river = "Susitna")
Sockeye_CD <- get_index(prey_ff, ylim = c(2005, 2015), species = "Sockeye", river = "Kenai") + get_index(prey_ff, ylim = c(2005, 2015), species = "Sockeye", river = "Kasilof")
# Posterior reproductive rates
Ps1 <- get_posterior_samples("Ps", c(2005, 2017), Chains, 10000, 100)
Ps2 <- get_posterior_samples("Ps", c(2005, 2017), Chains, 10000, 100)
Ps3 <- get_posterior_samples("Ps", c(2006, 2015), Chains, 10000, 100)
Ps4 <- get_posterior_samples("Ps", c(2005, 2015), Chains, 10000, 100)
# Inverse cdfs
cdf1 <- get_inverse_cdf(Ps1, Chinook_ND)
cdf2 <- get_inverse_cdf(Ps2, Chinook_CD)
cdf3 <- get_inverse_cdf(Ps3, Sockeye_ND)
cdf4 <- get_inverse_cdf(Ps4, Sockeye_CD)
# Names for plots
names <- c(
  "Northern District Chinook",
  "Central District Chinook",
  "Northern District Sockeye",
  "Central District Sockeye"
)
# Plot cdfs
plot_inverse_cdfs(cdf1, cdf2, cdf3, cdf4, names, par = "Ps") # Graphic for publication

plot_reproduction <- function(fp_results){
  rep_rates <- readxl::read_xlsx(fp_results, sheet = "Pr")
  year = 2005:2017
  ROR <- rep_rates$med[1:13]
  ROR_SD <- rep_rates$sd[1:13]
  par(mfrow = c(1,1))
  plot(
    year,
    ROR,
    ylim = c(0,0.25),
    xlab = NA,
    ylab = "Per Adult Reproduction Rate"
    #,main = "CIBW Reproduction Rate"
  )
  for (y in c(0.05, 0.1, 0.15, 0.2, 0.25)){
    lines(
      c(2004, 2019),
      c(y, y),
      col = "lightgrey"
    )
  }
  lines(
    year,
    ROR
  )
  for (i in 1:length(year)){
    yr = year[i]
    lines(
      c(yr, yr),
      c(ROR[i]-ROR_SD[i], ROR[i]+ROR_SD[i]),
      col="red"
    )
  }
  points(
    year,
    ROR,
    pch = 19
  )
}
plot_vital_rates <- function(par, med, sd, ylim){
  par(mfrow = c(1,1))
  if (par == "Ps"){
    ylab = "Rate of Survival to Next Year"
  }
  if (par == "Pr"){
    ylab = "Per-Adult Reproduction Rate"
  }
  plot(2005:2017, med, ty = "l", ylim = ylim, ylab = ylab, xlab = NA)
  if (par == "Ps"){
    for (m in c(.5, 0.6, 0.7, 0.8, 0.9, 0.9, 1)){
      lines(
        c(-4000, 4000),
        c(m, m),
        col = "lightgrey"
      )
    }
  }
  if (par == "Pr"){
    for (m in c(0, 0.1, 0.2, 0.3, 0.4, 0.5)){
      lines(
        c(-4000, 4000),
        c(m, m),
        col = "lightgrey"
      )
    }
  }
  for (i in 1:length(med)){
    year <- (2005:2017)[i]
    lines(
      c(year, year), 
      c(med[i]-sd[i], med[i]+sd[i]),
      col = "red",
      lwd = 2
    )
  }
  lines(2005:2017, med, lwd = 2)
  points(2005:2017, med, lwd = 2, pch = 19)
}
plot_vital_rates(
  par = "Pr",
  med = apply(get_posterior_samples("Pr", c(2005, 2017), Chains, 1000, 100), MARGIN = 2, mean),
  sd = apply(get_posterior_samples("Pr", c(2005, 2017), Chains, 1000, 100), MARGIN = 2, sd),
  ylim = c(0, 0.5) 
)
plot_vital_rates(
  par = "Ps",
  med <- apply(get_posterior_samples("Ps", c(2005, 2017), Chains, 1000, 100), MARGIN = 2, mean),
  sd <- apply(get_posterior_samples("Ps", c(2005, 2017), Chains, 1000, 100), MARGIN = 2, sd),
  ylim <- c(0.5, 1)
)