# Header
{
  # Load required packages
  {
    require(rstudioapi)
    require(readxl)
    require(writexl)
    require(foreach)
    require(doParallel)
    require(MCMCvis)
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

# Functions
{
  
  # MCMC routine (for multiple chains)
  run_MCMC <- function(data, inits, niter, nthin, nchain, parallel, dir_out, nbatch){
  
    # MCMC routine (for a single chain)
    MCMC <- function(data, inits, niter, nthin, fp_out){
      
      # number of posterior draws to keep
      ns <- ceiling(niter/nthin)
      
      # unpack data
      ny <- nrow(data)
      na <- 3
      Nt_hat <- data$Nt_hat
      SD_Nt_hat <- data$SD_Nt_hat
      Pa_hat <- cbind(
        data$P1_hat,
        data$P2_hat,
        data$P3_hat
      )
      SD_Pa_hat <- cbind(
        data$SD_P1_hat,
        data$SD_P2_hat,
        data$SD_P3_hat
      )
      
      # instantiate data structures to return
      Nt_post <- array(NA, dim = c(ns,ny))
      Na_post <- array(NA, dim = c(ns,ny,na))
      Pa_post <- array(NA, dim = c(ns,ny,na))
      Pr_post <- array(NA, dim = c(ns,ny))
      Ps_post <- array(NA, dim = c(ns,ny,na))
      Phi_post <- array(NA, dim = c(ns, ny))
      LL_post <- array(NA, dim = c(ns))
      
      for (imc in 1:niter){ # start MC iteration
        
        # instantiate proposed state
        Nt_prop <- array(NA, dim = c(ny))
        Na_prop <- array(NA, dim = c(ny,na))
        Pa_prop <- array(NA, dim = c(ny,na))
        Pr_prop <- array(NA, dim = c(ny))
        Ps_prop <- array(NA, dim = c(ny,na))
        Phi_prop <- array(NA, dim = c(ny))
        
        if (imc == 1){ # initialize proposed (parameter) state 
          Na_prop[1,1:3] <- inits$Na
          Pr_prop <- inits$Pr
          Ps_prop <- inits$Ps
          Phi_prop <- inits$Phi
          
        }else{ # propose new state
          for (ia in 1:na){
            Na_prop[1,ia] <- Na[1,ia] + runif(1, -5, 5)
          }
          for (iy in 1:ny){
            Pr_prop[iy] <- Pr[iy] + runif(1, -0.01, 0.01)
            for(ia in 1:na){
              Ps_prop[iy,ia] <- Ps[iy,ia] + runif(1, -0.01, 0.01)
            }
          }
          for (iy in 1:ny){
            Phi_prop[iy] <- Phi[iy] + runif(1, -0.01, 0.01)
          }
          
          
        }
        
        if ( # if proposal lies within support
          sum(
            sum(Na_prop[1,1:3] < 0) | sum(Na_prop[1,1:3] > 1000),
            sum(Pr_prop[1:ny] < 0 | Pr_prop[1:ny] > 0.5), 
            sum(Ps_prop[1:ny,1:na] < 0 | Ps_prop[1:ny,1:na] > 1),  
            sum(Phi_prop[1:ny] < 0 | Phi_prop[1:ny] > 1)
          ) == 0
            
        ){
          
          # determine latent state
          for (iy in 1:(ny-1)){
            
            # equation 2
            Na_prop[iy+1,2] <- Na_prop[iy,2]*Ps_prop[iy,1]+Na_prop[iy,2]*Ps_prop[iy,2]*(1-Phi_prop[iy])
            
            # equation 3
            Na_prop[iy+1,3] <- Na_prop[iy,2]*Ps_prop[iy,2]*Phi_prop[iy] + Na_prop[iy,3]*Ps_prop[iy,3]
            
            # equation 1
            Na_prop[iy+1,1] <- Na_prop[iy+1,3]/2*Pr_prop[iy+1]
            
          }
          
          # equation 4
          Nt_prop[1:ny] <- Na_prop[1:ny,2]+Na_prop[1:ny,3]
          
          # equation 5
          for (ia in 1:3){
            Pa_prop[1:ny,ia] <- Na_prop[1:ny,ia]/(Na_prop[1:ny,1]+Na_prop[1:ny,2]+Na_prop[1:ny,3])
          }
          
          
          # log likelihood
          LL_prop <- 0
          for (iy in which(!is.na(Nt_hat))){ # years where abundance was estimated
            LL_prop <- LL_prop + dnorm(Nt_hat[iy], Nt_prop[iy], SD_Nt_hat[iy], log = T)
          }
          for (ia in 1:na){
            for (iy in which(!is.na(Pa_hat[,1]))){ # years with photo-id mark-recapture data
              LL_prop <- LL_prop + dnorm(Pa_hat[iy,ia], Pa_prop[iy,ia], SD_Pa_hat[iy,ia], log = T)
            }
          }
          
          if (imc == 1){ # automatically accept proposal
            Nt <- Nt_prop
            Na <- Na_prop
            Pa <- Pa_prop
            Pr <- Pr_prop
            Ps <- Ps_prop
            Phi <- Phi_prop
            LL <- LL_prop
            
          }else{ # do Metropolis-Hastings
            p <- runif(1, 0, 1)
            mh <- exp(LL_prop - LL)
            if (!is.nan(LL_prop) & mh >= p){ # accept proposal
              Nt <- Nt_prop
              Na <- Na_prop
              Pa <- Pa_prop
              Pr <- Pr_prop
              Ps <- Ps_prop
              Phi <- Phi_prop
              LL <- LL_prop
              
            }
          }
          
        } # end if proposal lies within support
        
        # hold onto thinned data
        if (imc %% nthin == 1){
          if (imc == 1){
            k <- 1
          }else{
            k <- k + 1
          }
          Nt_post[k,] <- Nt
          Na_post[k,,] <- Na
          Pa_post[k,,] <- Pa
          Pr_post[k,] <- Pr
          Ps_post[k,,] <- Ps
          Phi_post[k,] <- Phi
          LL_post[k] <- LL
          
          print(round(imc/niter,2))
          
        }
        
        # save results so far to disk (in case computer crashes)
        if (imc %% nbatch == 0){
          print(k)
          MCMC_out <- list(
            Nt = Nt_post[1:k,],
            Na = Na_post[1:k,,],
            Pa = Pa_post[1:k,,],
            Pr = Pr_post[1:k,],
            Ps = Ps_post[1:k,,],
            Phi = Phi_post[1:k,],
            LL = LL_post[1:k]
          )
          save(
            MCMC_out,
            file = fp_out
          )
        }
        
        
      } # end MC iteration
      
      # save final results to disk
      MCMC_out <- list(
        Nt = Nt_post,
        Na = Na_post,
        Pa = Pa_post,
        Pr = Pr_post,
        Ps = Ps_post,
        Phi = Phi_post,
        LL = LL_post
      )
      save(
        MCMC_out,
        file = fp_out
      )
      
    } # end MCMC routine (for a single chain)
    
    if (isFALSE(parallel)){ # loop through chains and run MCMC routine
      a <- Sys.time()
      for (ic in 1:nchain){
        MCMC(
          data = data, 
          inits = inits, 
          niter = niter, 
          nthin = nthin, 
          fp_out = paste0(dir_out, "/chain", ic, ".Rdata")
        )
      }
      print(Sys.time()-a)
      
    }else{ # run MCMC in parallel
      a <- Sys.time()
      cl <- makePSOCKcluster(7)
      registerDoParallel(cl)
      foreach(ic=1:nchain) %dopar% MCMC(
        data = data, 
        inits = inits, 
        niter = niter, 
        nthin = nthin, 
        fp_out = paste0(dir_out, "/chain", ic, ".Rdata")
      )
      stopCluster(cl)
      print(Sys.time()-a)
      
    }
  }
  
  # format posteriors (so that the MCMCvis package can be used)
  format_post <- function(dir_out, burnin = NA){
    
    # instantiate list to return
    list_out <- list()
    
    # posterior samples for each chain
    files <- list.files(dir_out)
    files <- files[substr(files, nchar(files)-4, nchar(files)) == "Rdata"]
    
    for (ic in 1:length(files)){ # iterate through chains
      
      # load posterior samples
      load(paste0(dir_out, "/", files[ic]))
      
      # number of samples
      nr <- length(MCMC_out$LL)
      
      # number of parameters
      np <- 0
      for (il in 1:length(MCMC_out)){
        np <- np + length(MCMC_out[[il]])/nr
      }
      
      # instantiate posterior sample matrix
      mat_out <- matrix(NA, nrow = nr, ncol = np)
      
      # fill in matrix
      k <- 1
      for (il in 1:length(MCMC_out)){
        dims <- dim(MCMC_out[[il]])
        if (length(dims)==1){
          mat_out[,k] <- MCMC_out[[il]]
          k <- k + 1
        }else{
          mat <- matrix(MCMC_out[[il]], nrow = nr)
          for (im in 1:ncol(mat)){
            mat_out[,k] <- mat[,im]
            k <- k + 1
          }
        }
      }
      
      # get column names
      names <- c()
      for (il in 1:length(MCMC_out)){
        n <- names(MCMC_out)[il]
        dims <- dim(MCMC_out[[il]])[-1]
        if (length(dims) == 0){
          names <- c(names, n)
        }else{
          l = list()
          for (id in 1:length(dims)){
            l[[id]] <- 1:dims[id]
          }
          g <- expand.grid(l)
          for (ig in 1:nrow(g)){
            vn <- paste0(
              n,
              "[",
              paste0(g[ig,], collapse=","),
              "]"
            )
            names <- c(names, vn)
          }
        }
        
        # convert to data frame
        df_out <- as.data.frame(mat_out)
        names(df_out) <- names
        
      }
      
      # discard burnin, if applicable
      if (!is.na(burnin)){
        df_out <- df_out[(burnin+1):nrow(df_out),]
      }
      
      # add info to list
      list_out[[ic]] <- df_out
      
    }  
    return(list_out) 
  }
  
  # Save the results in a tabular format
  save_results <- function(fp_results, fp_out, burnin = 333){
    
    spost <- format_post(fp_results, burnin = burnin)
    results <- MCMCvis::MCMCsummary(spost)
    
    Nt <- results[substr(rownames(results),1,2)=="Nt",]
    Na <- results[substr(rownames(results),1,2)=="Na",]
    Pa <- results[substr(rownames(results),1,2)=="Pa",]
    Pr <- results[substr(rownames(results),1,2)=="Pr",]
    Ps <- results[substr(rownames(results),1,2)=="Ps",]
    Phi <- results[substr(rownames(results),1,3)=="Phi",]
    
    Nt <- data.frame(
      par = rownames(results)[substr(rownames(results),1,2)=="Nt"],
      mean = Nt$mean,
      med = Nt$`50%`,
      sd = Nt$sd
    )
    
    Na <- data.frame(
      par = rownames(results)[substr(rownames(results),1,2)=="Na"],
      mean = Na$mean,
      med = Na$`50%`,
      sd = Na$sd
    )
    
    Pa <- data.frame(
      par = rownames(results)[substr(rownames(results),1,2)=="Pa"],
      mean = Pa$mean,
      med = Pa$`50%`,
      sd = Pa$sd
    )
    
    Pr <- data.frame(
      par = rownames(results)[substr(rownames(results),1,2)=="Pr"],
      mean = Pr$mean,
      med = Pr$`50%`,
      sd = Pr$sd
    )
    
    Ps <- data.frame(
      par = rownames(results)[substr(rownames(results),1,2)=="Ps"],
      mean = Ps$mean,
      med = Ps$`50%`,
      sd = Ps$sd
    )
    
    Phi <- data.frame(
      par = rownames(results)[substr(rownames(results),1,3)=="Phi"],
      mean = Phi$mean,
      med = Phi$`50%`,
      sd = Phi$sd
    )
    
    out <- list(
      Nt = Nt,
      Na = Na,
      Pa = Pa,
      Pr = Pr,
      Ps = Ps,
      Phi = Phi
    )
    
    writexl::write_xlsx(out, path = fp_out)
    
    return(out)
    
  }
  
  # Some plotting functions
  {
    # Figure
    plot_abundance_raw <- function(data, ymin = 1994, ymax = 2022){
      data <- data[data$year >= ymin & data$year <= ymax, ]
      par(mfrow = c(1,1))
      plot(
        data$year,
        data$Nt_hat,
        ylim = c(0, 800),
        cex = 0.5,
        xlab = NA,
        ylab = "Abundance",
        xaxt = "n"
        #,main = "CIBW Abundance 1994 to 2022"
      )
      axis(1, at = ymin:ymax, labels = ymin:ymax)
      for (y in c(100, 200, 300, 400, 500, 600, 700)){
        lines(
          c(ymin-2, ymax+2),
          c(y, y),
          col = "lightgrey"
        )
      }
      lines(c(1999,1999), c(-1000, 1000), lty = 2)
      for (i in 1:nrow(data)){
        lines(
          c(data$year[i], data$year[i]),
          c(data$Nt_hat[i]-data$SD_Nt_hat[i], data$Nt_hat[i]+data$SD_Nt_hat[i]),
          col = "red"
        )
      }
      points(
        data$year,
        data$Nt_hat,
        pch = 19,
        cex = 0.75
      )
    }
    
    plot_age_comp <- function(fp_results){
      age_comp <- readxl::read_xlsx(fp_results, sheet = "Na")
      age_comp <- data.frame(
        year = 2005:2017,
        YOY = round(age_comp$med[1:13]),
        JUV = round(age_comp$med[15:27]),
        MAT = round(age_comp$med[29:41])
      )
      barplot(cbind(YOY, JUV, MAT) ~ year, data = age_comp, 
              ylim = c(0,500), 
              xlab = NA,
              legend.text = c("YOY", "Juvenile", "Mature"),
              col = c("#ADD8E6", "#4682B4", "#00008B"),
              ylab = "Abundance"
              #,main = "Beluga Abundance by Age Class"
      )
      lines(
        c(-4000, 4000),
        c(0, 0),
        col = "black"
      )
      for (n in c(100, 200, 300, 400, 500, 600, 700)){
        lines(
          c(-4000, 4000),
          c(n, n),
          col = "lightgrey"
        )
      }
      barplot(cbind(YOY, JUV, MAT) ~ year, data = age_comp, 
              ylim = c(0,500), 
              xlab = NA,
              legend.text = c("YOY", "Juvenile", "Mature"),
              col = c("#ADD8E6", "#4682B4", "#00008B"),
              add = T
      )
    }
  }
}

# Run MCMC routine
{
  
  out_dir <- "20250131"
  
  # Read in and format data
  {
    data_abundance <- readxl::read_xlsx("../Input Data/Beluga.xlsx", sheet = "Annual Abundance")
    data_harvest <- readxl::read_xlsx("../Input Data/Beluga.xlsx", sheet = "Annual Harvest") # Not used in report
    
    data_age_class_med <- readxl::read_xlsx("../Output Data/Beluga_Age_Class_Proportions.xlsx", sheet = "p_med")
    data_age_class_cv <- readxl::read_xlsx("../Output Data/Beluga_Age_Class_Proportions.xlsx",  sheet = "p_sd")
    
    ymin <- min(c(data_abundance$YR, data_harvest$YR, data_age_class_med$year))
    ymax <- max(c(data_abundance$YR, data_harvest$YR, data_age_class_med$year))
    ny <- length(ymin:ymax)
    data <- data.frame(
      year = ymin:ymax,
      y = 1:ny,
      Nt_hat = NA,
      SD_Nt_hat = NA,
      P1_hat = NA,
      SD_P1_hat = NA,
      P2_hat = NA,
      SD_P2_hat = NA,
      P3_hat = NA,
      SD_P3_hat = NA,
      H_min = NA,
      H_max = NA
    )
    for (y in data$y){
      year <- data$year[y]
      Nt_hat <- data_abundance$Nhat[data_abundance$YR == year]
      SD_Nt_hat <- data_abundance$SD_Nhat[data_abundance$YR == year] 
      P1_hat <- data_age_class_med$ac1[data_age_class_med$year == year]
      SD_P1_hat <- data_age_class_cv$ac1[data_age_class_cv$year == year]
      P2_hat <- data_age_class_med$ac2[data_age_class_med$year == year]
      SD_P2_hat <- data_age_class_cv$ac2[data_age_class_cv$year == year]
      P3_hat <- data_age_class_med$ac3[data_age_class_med$year == year]
      SD_P3_hat <- data_age_class_cv$ac3[data_age_class_cv$year == year]
      H_min <- data_harvest$Hmin[data_harvest$YR == year]
      H_max <- data_harvest$Hmax[data_harvest$YR == year]
      if (length(Nt_hat == 1))  data$Nt_hat[y] <- Nt_hat
      if (length(SD_Nt_hat)==1) data$SD_Nt_hat[y] <- SD_Nt_hat
      if (length(P1_hat)==1) data$P1_hat[y] <- P1_hat
      if (length(SD_P1_hat)==1) data$SD_P1_hat[y] <- SD_P1_hat
      if (length(P2_hat)==1) data$P2_hat[y] <- P2_hat
      if (length(SD_P2_hat)==1) data$SD_P2_hat[y] <- SD_P2_hat
      if (length(P3_hat)==1) data$P3_hat[y] <- P3_hat
      if (length(SD_P3_hat)==1) data$SD_P3_hat[y] <- SD_P3_hat
      if (length(H_min)==1) data$H_min[y] <- H_min
      if (length(H_max)==1) data$H_max[y] <- H_max
    }
    # Figure
    plot_abundance_raw(data)
    # writexl::write_xlsx(data, path = "../Xlsx/IPM_Input_Data_20240813.xlsx")
    data <- data[12:25,]
    data$y <- 1:nrow(data)
    ny <- nrow(data)
    print(data)
  }
  # Set initial values
  #    Some parameters were initialized using the results of Himes-Boor et al. 2020.
  #    Relevant results are outlined below:
  #       birth rate for established breeders -- 0.279
  #       YOY-1 yo calf survival -- 0.926
  #       adult nonbreeder survival -- 0.931
  #       adult breeder survival -- 0.962
  {
    inits <- list(
      Na = round(data$Nt_hat[1]*c(mean(data$P1_hat,na.rm = T), mean(data$P2_hat,na.rm = T), mean(data$P3_hat,na.rm = T))),
      Pr = rep(0.279*0.33*2, ny), # From Himes-Boor et al. 2020
      Ps = cbind( # From Himes-Boor et al. 2020
        rep(0.926, ny),
        rep((0.926+0.931+0.962)/3, ny),
        rep((0.931+0.962)/2, ny)
      ),
      H = round((data$H_min + data$H_max)/2),
      Phi = rep(0.5, ny)
    )
    print(inits)
  }
  
  # Run MCMC routine (should run for 5 days saving the results every 6 hours)
  a <- Sys.time()
  run_MCMC(
    data, 
    inits, 
    niter = 10000000,
    nthin = 10000,
    nchain = 4, 
    parallel = T, 
    dir_out = paste0("../Output Data/Vital Rates/", out_dir),
    nbatch = 1000000
  )
  print(Sys.time()-a)
  
  # Format posteriors
  MCMC_post <-  format_post(paste0("../Output Data/Vital Rates/", out_dir))
  
  # Save traceplots to disk
  MCMCvis::MCMCtrace(MCMC_post, wd = paste0("../Output Data/Vital Rates/", out_dir))
  
  # Save results in a tabular form
  fp_results <- paste0("../Output Data/Vital Rates/", out_dir)
  fp_out <- paste0("../Output Data/Vital Rates/", out_dir, "/Results.xlsx")
  save_results(fp_results, fp_out)
}
{
  # Figure 
  plot_age_comp(paste0("../Output Data/Vital Rates/", out_dir, "/Results.xlsx"))
}
# Create summary table for report
{
  Nt_y <- readxl::read_xlsx(paste0("../Output Data/Vital Rates/", out_dir, "/Results.xlsx"), sheet = "Nt")
  Na_ya <- readxl::read_xlsx(paste0("../Output Data/Vital Rates/", out_dir, "/Results.xlsx"), sheet = "Na")
  Pa_ya <- readxl::read_xlsx(paste0("../Output Data/Vital Rates/", out_dir, "/Results.xlsx"), sheet = "Pa")
  df1 <-  data.frame(
    year = 2005:2018,
    Nt_y = paste0(
      round(Nt_y$med),
      " (",
      round(Nt_y$sd),
      ")"
    ),
    Na_y1 = paste0(
      round(Na_ya$med[1:14]),
      " (",
      round(Na_ya$sd[1:14]),
      ")"
    ),
    Na_y2 = paste0(
      round(Na_ya$med[15:28]),
      " (",
      round(Na_ya$sd[15:28]),
      ")"
    ),
    Na_y3 = paste0(
      round(Na_ya$med[29:42]),
      " (",
      round(Na_ya$sd[29:42]),
      ")"
    ),
    Pa_y1 = paste0(
      round(Pa_ya$med[1:14],2),
      " (",
      round(Pa_ya$sd[1:14],2),
      ")"
    ),
    Pa_y2 = paste0(
      round(Pa_ya$med[15:28],2),
      " (",
      round(Pa_ya$sd[15:28],2),
      ")"
    ),
    Pa_y3 = paste0(
      round(Pa_ya$med[29:42],2),
      " (",
      round(Pa_ya$sd[29:42],2),
      ")"
    )
  )
  Pr_y <- readxl::read_xlsx(paste0("../Output Data/Vital Rates/", out_dir, "/Results.xlsx"), sheet = "Pr")
  Ps_ya <- readxl::read_xlsx(paste0("../Output Data/Vital Rates/", out_dir, "/Results.xlsx"), sheet = "Ps")
  Phi_y <- readxl::read_xlsx(paste0("../Output Data/Vital Rates/", out_dir, "/Results.xlsx"), sheet = "Phi")
  df2 <- data.frame(
    year = 2005:2018,
    Pr_y = paste0(
      round(Pr_y$med,2),
      " (",
      round(Pr_y$sd,2),
      ")"
    ),
    Ps_y1 = paste0(
      round(Ps_ya$med[1:14],2),
      " (",
      round(Ps_ya$sd[1:14],2),
      ")"
    ),
    Ps_y2 = paste0(
      round(Ps_ya$med[15:28],2),
      " (",
      round(Ps_ya$sd[15:28],2),
      ")"
    ),
    Ps_y3 = paste0(
      round(Ps_ya$med[29:42],2),
      " (",
      round(Ps_ya$sd[29:42],2),
      ")"
    ),
    Phi_y = paste0(
      round(Phi_y$med,2),
      " (",
      round(Phi_y$sd,2),
      ")"
    ) 
  )
  out <- list(
    df1,
    df2
  )
  names(out) <- c("df1", "df2")
  writexl::write_xlsx(out, path = "../Output Data/Population_Parameters.xlsx")
}
