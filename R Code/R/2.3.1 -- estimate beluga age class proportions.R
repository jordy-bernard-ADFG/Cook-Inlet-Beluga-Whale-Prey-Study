#####
# This script is used to estimate the proportion of belugas by age class by year from Himes-Boor et al 2022 photo-id mark-recapture data
#   Note:
#     Age class 1 is YOY
#     Age class 2 is calves that are not YOY
#     Age class 3 is mature belugas
#####

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

# Read in photo-id mark-recapture data
mr_data <- readxl::read_xlsx("../Input Data/Beluga.xlsx", sheet = "Annual Mark Recapture")
mr_data <- as.data.frame(mr_data)

# Himes-Boor et al 2022 encoding shown below:
#   1 = 0 not seen     
#   2 = P seen alone     
#   3 = J0 seen with YOY     
#   4 = J1- seen with a calf thought to be either a YOY or 1yo
#   5 = J1 seen with 1yo calf (no uncertainty in calf age)     
#   6 = J1+ seen with a calf thought to be 1yo or older)     
#   7 = J2- seen with a calf thought to be a YOY, 1yo, or 2yo    
#   8 = J2 seen with 2yo calf (no uncertainty in calf age)   
#   9 = J2+ seen with a calf thought to be 2yo or older     
#   10 = J3 seen with 3yo calf (no uncertainty in calf age)   
#   11 = J3+ seen with a calf thought to be 3yo or older     
#   12 = J4 seen with a 4yo calf     
#   13 = J4+ seen with a calf thought to be 4yo or older    
#   14 = Unk seen with a calf of indeterminate age       
#   15 = J1+ J0 seen with 2 calves, one a YOY (J0) and one thought to be 1yo or older   
#   16 = J1+ J1- seen with 2 calves, one thought to be 1yo or older and one thought to be 1yo or younger
#   17 = J1+ J1 etc
#   18 = J2- J0	
#   19 = J2- J1-	
#   20 = J2 J0	
#   21 = J2 J1-	
#   22 = J2 J2-	
#   23 = J2+ J0	
#   24 = J2+ J1-	
#   25 = J2+ J1	
#   26 = J2+ J1+	
#   27 = J2+ J2-	
#   28 = J2+ J2	
#   29 = J3 J0	
#   30 = J3 J1-	
#   31 = J3 J1	
#   32 = J3 J1+	
#   33 = J3 J2-	
#   34 = J3+ J0	
#   35 = J3+ J1-	
#   36 = J3+ J1	
#   37 = J3+ J1+	
#   38 = J3+ J2-	
#   39 = J3+ J2	
#   40 = J3+ J2+	
#   41 = J4 J0	
#   42 = J4 J1-	
#   43 = J4 J1	
#   44 = J4 J1+	
#   45 = J4 J2-	
#   46 = J4 J2	
#   47 = J4 J2+	
#   48 = J4+ J0	
#   49 = J4+ J1-     	
#   50 = J4+ J1     	
#   51 = J4+ J1+	     
#   52 = J4+ J2-	     
#   53 = J4+ J2	     
#   54 = J4+ J2+	     
#   55 = J0 Unk	     
#   56 = J1- Unk	     
#   57 = J1 Unk	     
#   58 = J1+ Unk	     
#   59 = J2- Unk	     
#   60 = J2 Unk	     
#   61 = J2+ Unk	     
#   62 = J3 Unk	     
#   63 = J3+ Unk	     
#   64 = J4 Unk	     
#   65 = J4+ Unk	     
#   66 = Unk Unk
#   67 = J1+ J1+
#   68 = J1+ J2-
#   69 = J1+ J2
#   70 = J1+ J2+
#   71 = J2- J2-
#   72 = J2+ J2+ 

# Lookup table to convert Himes-Boor et al 2022 observation codes to the number of belugas observed by age class
LOOKUP <- data.frame(
  code = 1:72, # Himes-Boor et al 2022 code
  symb = NA, # Himes-Boor et al 2022 symbol
  ac1c = NA, # Number belugas age class 1 observed with certainty
  ac1u = NA, # Number belugas age class 1 including unclear observations
  ac2c = NA, # Number belugas age class 2 observed with certainty
  ac2u = NA, # Number belugas age class 2 including unclear observations
  ac3c = NA, # Number belugas age class 3 observed with certainty
  ac3u = NA  # Number belugas age class 3 including unclear observations
)

# Lookup table for codes that are not composites
lookup_1 <- data.frame(
  code = 1:14,
  symb = c("0", "P", "J0", "J1-", "J1", "J1+", "J2-", "J2", "J2+", "J3", "J3+", "J4", "J4+", "Unk"),
  ac1c = c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), # Age class 1 observed with certainty
  ac1u = c(0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1), # Age class 1 including uncertain observations
  ac2c = c(0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0), # Age class 2 observed with certainty
  ac2u = c(0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) # Age class 2 including uncertain observations
)

# Lookup table for composite codes
lookup_2 <- data.frame(
  code = 15:72,
  symb1 = c(
    "J1+","J1+","J1+","J2-","J2-",	
    "J2","J2","J2","J2+","J2+",
    "J2+","J2+","J2+","J2+","J3",
    "J3","J3","J3","J3","J3+",	
    "J3+","J3+","J3+","J3+","J3+",
    "J3+","J4","J4","J4","J4",
    "J4","J4","J4","J4+","J4+", 	
    "J4+","J4+","J4+","J4+","J4+",    
    "J0","J1-","J1","J1+","J2-",    
    "J2","J2+","J3","J3+","J4",	
    "J4+","Unk","J1+","J1+","J1+",
    "J1+","J2-","J2+"
  ),
  symb2 = c(
    "J0","J1-","J1","J0","J1-",	
    "J0","J1-","J2-","J0","J1-",	
    "J1","J1+","J2-","J2","J0",
    "J1-","J1","J1+","J2-","J0",
    "J1-","J1","J1+","J2-","J2",
    "J2+","J0","J1-","J1","J1+",
    "J2-","J2","J2+","J0","J1-",    	
    "J1","J1+","J2-","J2","J2+",	     
    "Unk","Unk","Unk","Unk","Unk",	     
    "Unk","Unk","Unk","Unk","Unk",	     
    "Unk","Unk","J1+","J2-","J2",
    "J2+","J2-","J2+" 
  )
)

# Populate primary lookup table
LOOKUP[1:14,1:6] <- lookup_1
for (code in 15:72){
  symb1 <- lookup_2$symb1[lookup_2$code == code]
  symb2 <- lookup_2$symb2[lookup_2$code == code]
  LOOKUP$symb[code] <- paste(symb1, symb2)
  LOOKUP$ac1c[code] <- lookup_1$ac1c[lookup_1$symb == symb1] + lookup_1$ac1c[lookup_1$symb == symb2]
  LOOKUP$ac1u[code] <- lookup_1$ac1u[lookup_1$symb == symb1] + lookup_1$ac1u[lookup_1$symb == symb2]
  LOOKUP$ac2c[code] <- lookup_1$ac2c[lookup_1$symb == symb1] + lookup_1$ac2c[lookup_1$symb == symb2]
  LOOKUP$ac2u[code] <- lookup_1$ac2u[lookup_1$symb == symb1] + lookup_1$ac2u[lookup_1$symb == symb2]
}
LOOKUP[,7] <- c(0, rep(1,71))
LOOKUP[,8] <- c(0, rep(1,71))

# Compute the number of belugas observed by age class and year
Counts_ac1c <- Counts_ac2c <- Counts_ac3c <- Counts_ac1u <- Counts_ac2u <- Counts_ac3u <- as.data.frame(matrix(nrow = nrow(mr_data), ncol = ncol(mr_data)))
for (i in 1:nrow(mr_data)){
  for (j in 1:ncol(mr_data)){
    Counts_ac1c[i,j] <- LOOKUP$ac1c[LOOKUP$code == mr_data[i,j]]
    Counts_ac2c[i,j] <- LOOKUP$ac2c[LOOKUP$code == mr_data[i,j]]
    Counts_ac3c[i,j] <- LOOKUP$ac3c[LOOKUP$code == mr_data[i,j]]
    Counts_ac1u[i,j] <- LOOKUP$ac1u[LOOKUP$code == mr_data[i,j]]
    Counts_ac2u[i,j] <- LOOKUP$ac2u[LOOKUP$code == mr_data[i,j]]
    Counts_ac3u[i,j] <- LOOKUP$ac3u[LOOKUP$code == mr_data[i,j]]
  }
}

# data structures used in modelling
M_min <- data.frame( # Minimum number of individuals of a given age class observed
  year = 2005:2017,
  ac1 = apply(Counts_ac1c, MARGIN = 2, sum),
  ac2 = apply(Counts_ac2c, MARGIN = 2, sum),
  ac3 = apply(Counts_ac3c, MARGIN = 2, sum)
)
M_max <- data.frame( # Maximum number of individuals of a given age class observed
  year = 2005:2017,
  ac1 = apply(Counts_ac1u, MARGIN = 2, sum),
  ac2 = apply(Counts_ac2u, MARGIN = 2, sum),
  ac3 = apply(Counts_ac3u, MARGIN = 2, sum)
)
  
# plot stuff out
par(mfrow = c(3,1))
plot(M_max$year, M_max$ac1, ty = "l", ylim = c(0,30))
lines(M_min$year, M_min$ac1)
plot(M_max$year, M_max$ac2, ty = "l", ylim = c(0,55))
lines(M_min$year, M_min$ac2)
plot(M_max$year, M_max$ac3, ty = "l", ylim = c(0,250))  
  
#' Monte Carlo routine to estimate the proportion of belugas by age class each year
#' Function factors in uncertainty in the number of individuals observed by age class and uncertainty in the probability of detection
#' @param M_min Minimum number of individuals of a given age class observed
#' @param M_max Maximum number of individuals of a given age class observed 
#' @param p1_hat Detection probability for YOY; The detection probability estimated by Himes-Boor et al 2022 for individuals <= 1 year old is used here, since one was not estimated for YOY individuals.
#' @param p2_hat Detection probability for calves that are not YOY; The detection probability estimated by Himes-Boor et al 2022 for individuals > 1 year old is used here.
#' @param sd_p1_hat Standard deviation of p1_hat
#' @param sd_p2_hat Standard devivation of p2_hat
#' @param niter Number of Monte Carlo Iterations
estimate_proportion_belugas_by_age_class <- function(M_min = M_min[,2:4], M_max = M_max[,2:4], p1_hat = 0.546, p2_hat = 0.473, sd_p1_hat = 0.029, sd_p2_hat = 0.050, niter = 1000){
  
  # Instantiate array to hold output
  out <- array(NA, dim = c(niter, nrow(M_min), ncol(M_min)))
  
  for (imc in 1:niter){ # Start MC integration
    
    # Number of individuals observed by year and age class
    M <- matrix(NA, nrow = nrow(M_min), ncol = ncol(M_max))
    for (iy in 1:nrow(M_min)){
      for (ia in 1:3){
        M[iy, ia] <- runif(1, M_min[iy, ia], M_max[iy, ia])
      }
    }
    
    # Probability of detecting YOY
    p1 <- rnorm(1, p1_hat, sd_p1_hat)
    
    # Probability of detection calves older that YOY
    p2 <- rnorm(1, p2_hat, sd_p2_hat)
    
    #if (p1 < 0) p1 <- 0.00001
    #if (p1 > 1) p1 <- 1
    #if (p2 < 0) p2 <- 0.00001
    #if (p2 > 1) p2 <- 1
    
    # Detection probability correction
    M[,1] <- M[,1]/p1
    M[,2] <- M[,2]/p2
    
    # Proportion of belugas in each age class by year
    q1 <- M[,1]/(M[,1]+M[,2]+M[,3]) 
    q2 <- M[,2]/(M[,1]+M[,2]+M[,3])
    q3 <- M[,3]/(M[,1]+M[,2]+M[,3])
    
    # Hold onto results
    out[imc,,1] <- q1
    out[imc,,2] <- q2
    out[imc,,3] <- q3
    
  }
  return(out)
}

# Run MC routine
MC_out <- estimate_proportion_belugas_by_age_class(M_min = M_min[,2:4], M_max = M_max[,2:4], p1_hat = 0.546, p2_hat = 0.473, sd_p1_hat = 0.029, sd_p2_hat = 0.050, niter = 10000)

# Summarize results
p_med <- as.data.frame(cbind(2005:2017, apply(MC_out, c(2,3), median)))
p_sd <- as.data.frame(cbind(2005:2017, apply(MC_out, c(2,3), sd)))
names(p_med) <- names(p_sd) <- c("year", "ac1", "ac2", "ac3")

# Write results to disk
output <- list(p_med = p_med, p_sd = p_sd)
writexl::write_xlsx(output, path = "../Output Data/Beluga_Age_Class_Proportions.xlsx")

# Plot stuff out
for (i in 1:25){
  s <- sample(1:10000, size = 1)
  MC_draw <- MC_out[s,,]
  if (i == 1){
    plot(2005:2017, MC_draw[,1], ty = "l", ylim = c(0, 0.20), col = rgb(0,0,1,0.5), main = "Proportion YOY", xlab = NA, ylab = "proportion")
  }else{
    lines(2005:2017, MC_draw[,1], ty = "l", col = rgb(0,0,1,0.5))
  }
}
lines(p_med$year, p_med$ac1, ty = "l", lwd=3)
points(p_med$year, p_med$ac1, pch=19, cex=2)
for (i in 1:25){
  s <- sample(1:10000, size = 1)
  MC_draw <- MC_out[s,,]
  if (i == 1){
    plot(2005:2017, MC_draw[,2], ty = "l", ylim = c(0, 0.40), col = rgb(1,0,0,0.5), main = "Proportion Calves Age 1 and Older", xlab = NA, ylab = "proportion")
  }else{
    lines(2005:2017, MC_draw[,2], ty = "l", col = rgb(1,0,0,0.5))
  }
}
lines(p_med$year, p_med$ac2, ty = "l", lwd = 3)
points(p_med$year, p_med$ac2, pch = 19, cex=2)
for (i in 1:25){
  s <- sample(1:10000, size = 1)
  MC_draw <- MC_out[s,,]
  if (i == 1){
    plot(2005:2017, MC_draw[,3], ty = "l", ylim = c(0.5, 1), col = rgb(0,1,0,0.5), main = "Proportion Mature Belugas", xlab = NA, ylab = "proportion")
  }else{
    lines(2005:2017, MC_draw[,3], ty = "l", col = rgb(0,1,0,0.5))
  }
}
lines(p_med$year, p_med$ac3, ty = "l", lwd = 3)
points(p_med$year, p_med$ac3, pch = 19, cex=2)
