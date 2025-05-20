# Clear workspace and set working directory
{
  rm(list=objects())
  this_fp <- function(){return(rstudioapi::getActiveDocumentContext()$path)}
  this_dir <- function(){return(dirname(this_fp()))}
  setwd(this_dir())
}
# Back out average vital rates from IPM output
get_IPM_averages <- function(fp_chain, burnin){
  load(fp_chain)
  Na1 <- MCMC_out$Na[(burnin+1):nrow(MCMC_out$Nt),-14,1]
  Na2 <- MCMC_out$Na[(burnin+1):nrow(MCMC_out$Nt),-14,2]
  Na3 <- MCMC_out$Na[(burnin+1):nrow(MCMC_out$Nt),-14,3]
  Pr <- MCMC_out$Pr[(burnin+1):nrow(MCMC_out$Nt),-14]
  Ps1 <- MCMC_out$Ps[(burnin+1):nrow(MCMC_out$Nt),-14,1]
  Ps2 <- MCMC_out$Ps[(burnin+1):nrow(MCMC_out$Nt),-14,2]
  Ps3 <- MCMC_out$Ps[(burnin+1):nrow(MCMC_out$Nt),-14,3]
  min(apply((Na1*Ps1+Na2*Ps2)/(Na1+Na2), MARGIN = 2, FUN = median))
  max(apply((Na1*Ps1+Na2*Ps2)/(Na1+Na2), MARGIN = 2, FUN = median))
  med = c(
    round(median(apply(Pr, MARGIN = 1, FUN = mean)), 2),
    round(median(apply((Na1*Ps1 + Na2*Ps2 + Na3*Ps3)/(Na1+Na2+Na3), MARGIN = 1, FUN = mean)),2)
  )
  sd = c(
    round(sd(apply(Pr, MARGIN = 1, FUN = mean)),3),
    round(sd(apply((Na1*Ps1 + Na2*Ps2 + Na3*Ps3)/(Na1+Na2+Na3), MARGIN = 1, FUN = mean)),3)
  )
  lcl_95 <- c(
    round(quantile(apply(Pr, MARGIN = 1, FUN = mean), probs = 0.025), 2),
    round(quantile(apply((Na1*Ps1 + Na2*Ps2 + Na3*Ps3)/(Na1+Na2+Na3), MARGIN = 1, FUN = mean), 0.025),2)
  )
  ucl_95 <- c(
    round(quantile(apply(Pr, MARGIN = 1, FUN = mean), probs = 0.975), 2),
    round(quantile(apply((Na1*Ps1 + Na2*Ps2 + Na3*Ps3)/(Na1+Na2+Na3), MARGIN = 1, FUN = mean), 0.975),2)
  )
  out <- data.frame(
    var = c("Pr", "Ps"),
    est = med,
    sd = sd,
    lcl_95 = lcl_95,
    ucl_95 = ucl_95
  )
  return(out)
}
get_IPM_averages("../Output Data/Vital Rates/20250114b/chain1.Rdata", burnin = 333)