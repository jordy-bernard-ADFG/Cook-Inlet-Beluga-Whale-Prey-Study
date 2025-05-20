# Clear workspace and set working directory
{
  rm(list=objects())
  this_fp <- function(){return(rstudioapi::getActiveDocumentContext()$path)}
  this_dir <- function(){return(dirname(this_fp()))}
  setwd(this_dir())
}
# Read in and restructure run size data
{
  Nt <- readxl::read_xlsx("../Input Data/Salmon.xlsx", sheet = "Annual_Run_Size")
  Kenai_Chinook_Early <- data.frame(
    Year = Nt$Year[Nt$Species=="Chinook" & Nt$Stock=="Kenai River Early"],
    Nt = Nt$Nt[Nt$Species=="Chinook" & Nt$Stock=="Kenai River Early"]
  )
  Kenai_Chinook_Late <- data.frame(
    Year = Nt$Year[Nt$Species=="Chinook" & Nt$Stock=="Kenai River Late"],
    Nt = Nt$Nt[Nt$Species=="Chinook" & Nt$Stock=="Kenai River Late"]
  )
  Kenai_Chinook <- data.frame(
    Year = Nt$Year[Nt$Species=="Chinook" & Nt$Stock=="Kenai River Early"],
    Nt = Nt$Nt[Nt$Species=="Chinook" & Nt$Stock=="Kenai River Early"] + Nt$Nt[Nt$Species=="Chinook" & Nt$Stock=="Kenai River Late"]
  )
  Susitna_Chinook <- data.frame(
    Year = Nt$Year[Nt$Species=="Chinook" & Nt$Stock == "Susitna River"],
    Nt = Nt$Nt[Nt$Species=="Chinook" & Nt$Stock == "Susitna River"]
  )  
  Kasilof_Sockeye <- data.frame(
    Year = Nt$Year[Nt$Species=="Sockeye" & Nt$Stock == "Kasilof River"],
    Nt = Nt$Nt[Nt$Species=="Sockeye" & Nt$Stock == "Kasilof River"]
  )
  Kenai_Sockeye <- data.frame(
    Year = Nt$Year[Nt$Species=="Sockeye" & Nt$Stock == "Kenai River"],
    Nt = Nt$Nt[Nt$Species=="Sockeye" & Nt$Stock == "Kenai River"]
  )
  Susitna_Sockeye <- data.frame(
    Year = Nt$Year[Nt$Species=="Sockeye" & Nt$Stock == "Susitna River"],
    Nt = Nt$Nt[Nt$Species=="Sockeye" & Nt$Stock == "Susitna River"]
  )
}
# Figure
{
  par(mfrow = c(2,1), mar = c(2, 4, 4, 2))
  plot(Kenai_Chinook$Year, Kenai_Chinook$Nt, xlim = c(1968, 2022), ylim = c(0, 250000), ty = "l", col = "red", lwd = 2, xlab = NA, ylab = "number of fish", main = "Chinook Annual Run Size")
  for (n in c(0, 50000, 100000, 150000, 200000, 250000)){
    lines(
      c(1950, 2050),
      c(n, n),
      col = "lightgrey"
    )
  }
  lines(Kenai_Chinook$Year, Kenai_Chinook$Nt,col = "red", lwd = 2)
  points(Kenai_Chinook$Year, Kenai_Chinook$Nt, pch = 19, col = "red")
  lines(Susitna_Chinook$Year, Susitna_Chinook$Nt, col = "blue", lwd = 2)
  points(Susitna_Chinook$Year, Susitna_Chinook$Nt, pch = 19, col = "blue")
  legend(
    x = "topright",
    legend = c(
      "Kenai",
      "Susitna"
    ),
    col = c("red", "blue"),
    lty = 1,
    lwd = 2
  )
  plot(Kasilof_Sockeye$Year, Kasilof_Sockeye$Nt/1000000, xlim = c(1968, 2022), ylim = c(0, 10000000/1000000), ty = "l", col = "green", lwd = 2, xlab = NA, ylab = "number of fish (in millions)", main = "Sockeye Annual Run Size")
  for (n in c(0, 2,4,6,8,10)){
    lines(
      c(1950, 2050),
      c(n, n),
      col = "lightgrey"
    )
  }
  lines(Kasilof_Sockeye$Year, Kasilof_Sockeye$Nt/1000000, col = "green", lwd = 2)
  points(Kasilof_Sockeye$Year, Kasilof_Sockeye$Nt/1000000, pch = 19, col = "green")
  lines(Kenai_Sockeye$Year, Kenai_Sockeye$Nt/1000000, ty = "l", col = "red", lwd = 2)
  points(Kenai_Sockeye$Year, Kenai_Sockeye$Nt/1000000, pch = 19, col = "red")
  lines(Susitna_Sockeye$Year, Susitna_Sockeye$Nt/1000000, col = "blue", lwd = 2)
  points(Susitna_Sockeye$Year, Susitna_Sockeye$Nt/1000000, pch = 19, col = "blue")
  legend(
    x = "topright",
    legend = c(
      "Kasilof",
      "Kenai",
      "Susitna"
    ),
    col = c("green", "red", "blue"),
    lty = 1,
    lwd = 2
  )
}
# Read in and restructure age-comp and age-length data
{
  #' Convert fish length to fish weight
  #' Function implements the predictive equations presented by Oke et al. 2020
  #' @param species Chinook, Chum, Coho, or Sockeye
  #' @param length mid-eye to fork length (in mm)
  #' @param return fish weight (in kg) 
  length_to_weight <- function(species, length){
    weight <- NA
    if (species == "Chinook"){
      weight <- 1.77*10^-9*length^3.30
    }
    if (species == "Chum"){
      weight <- 6.30*10^-9*length^3.14
    }
    if (species == "Coho"){
      weight <- 6.46*10^-9*length^3.14
    }
    if (species == "Sockeye"){
      weight <- 1.01*10^-8*length^3.10
    }
    return(weight)
  }
  #' Convert mass to energy at the individual level
  #' Function implements the predictive equations presented in O'Neil et al. 2014
  #' @param species Chinook, Chum, Coho, Pink, or Sockeye
  #' @param fish_weight fish weight in kg
  #' @return the energy content of the fish (in kcal)
  mass_to_energy <- function(species, fish_weight){
    energy <- NA
    if (species %in% c("Chinook", "Sockeye")){
      energy <- exp(0.94*log(fish_weight)+7.56)
    }
    if (species == "Coho"){
      energy <- exp(0.94*log(fish_weight)+7.31)
    }
    if (species %in% c("Pink", "Chum")){
      energy <- exp(0.94*log(fish_weight)+7.04)
    }
    return(energy)
  }
  # Read in Deshka Age, Sex, Length data
  deshka_asl <- readxl::read_xlsx("../Input Data/Salmon.xlsx", sheet = "ASL")
  # Compute Age, Length, Weight, and Energy Relationship
  alwe <- data.frame(
    age = c("1.1", "1.2", "1.3", "1.4"),
    age_numeric = c(3, 4, 5, 6),
    length = NA,
    weight = NA,
    energy = NA
  )
  for (i in 1:nrow(alwe)){
    age <- alwe$age[i]
    alwe$length[i] <- mean(deshka_asl$Length[deshka_asl$Age == age])
    alwe$weight[i] <- round(length_to_weight("Chinook", alwe$length[i]), 2)
    alwe$energy[i] <- round(mass_to_energy("Chinook", alwe$weight[i]))
  }
  # Read in and reformat age-comp data
  age_comp <- readxl::read_xlsx("../Input Data/Salmon.xlsx", sheet = "Age_Comp")
  # Get age comp for kenai chinook
  {
    age_comp_chinook_kenai_early <- age_comp[age_comp$Species=="Chinook" & age_comp$Stock=="Kenai Early",]
    age_comp_chinook_kenai_late <- age_comp[age_comp$Species=="Chinook" & age_comp$Stock=="Kenai Late",]
    age_comp_chinook_kenai <- data.frame(
      Year = Kenai_Chinook_Early$Year,
      Age_3 = (Kenai_Chinook_Early$Nt*age_comp_chinook_kenai_early$Age_3 + Kenai_Chinook_Late$Nt*age_comp_chinook_kenai_late$Age_3)/(Kenai_Chinook_Early$Nt+Kenai_Chinook_Late$Nt),
      Age_4 = (Kenai_Chinook_Early$Nt*age_comp_chinook_kenai_early$Age_4 + Kenai_Chinook_Late$Nt*age_comp_chinook_kenai_late$Age_4)/(Kenai_Chinook_Early$Nt+Kenai_Chinook_Late$Nt),
      Age_5 = (Kenai_Chinook_Early$Nt*age_comp_chinook_kenai_early$Age_5 + Kenai_Chinook_Late$Nt*age_comp_chinook_kenai_late$Age_5)/(Kenai_Chinook_Early$Nt+Kenai_Chinook_Late$Nt),
      Age_6 = (Kenai_Chinook_Early$Nt*age_comp_chinook_kenai_early$Age_6 + Kenai_Chinook_Late$Nt*age_comp_chinook_kenai_late$Age_6)/(Kenai_Chinook_Early$Nt+Kenai_Chinook_Late$Nt)
    )
  }
  # Age comp for susitna chinook
  age_comp_chinook_susitna <- age_comp[age_comp$Species=="Chinook" & age_comp$Stock=="Susitna",]
  # Estimate chinook biomass and energy content
  Kenai_Chinook$energy <- Susitna_Chinook$energy <- Kenai_Chinook$biomass <- Susitna_Chinook$biomass  <- NA
  for (iy in 1:nrow(Kenai_Chinook)){
    y <- Kenai_Chinook$Year[iy]
    Nt <- Kenai_Chinook$Nt[iy]
    ac <- c(
      age_comp_chinook_kenai$Age_3[age_comp_chinook_kenai$Year==y],
      age_comp_chinook_kenai$Age_4[age_comp_chinook_kenai$Year==y],
      age_comp_chinook_kenai$Age_5[age_comp_chinook_kenai$Year==y],
      age_comp_chinook_kenai$Age_6[age_comp_chinook_kenai$Year==y]
    )
    wt <- alwe$weight
    e <- mass_to_energy("Chinook", wt)
    Kenai_Chinook$biomass[iy] <- sum(Nt*ac*wt)
    Kenai_Chinook$energy[iy] <- sum(Nt*ac*e)
  }
  for (iy in 1:nrow(Susitna_Chinook)){
    y <- Susitna_Chinook$Year[iy]
    Nt <- Susitna_Chinook$Nt[iy]
    ac <- c(
      age_comp_chinook_susitna$Age_3[age_comp_chinook_susitna$Year==y],
      age_comp_chinook_susitna$Age_4[age_comp_chinook_susitna$Year==y],
      age_comp_chinook_susitna$Age_5[age_comp_chinook_susitna$Year==y],
      age_comp_chinook_susitna$Age_6[age_comp_chinook_susitna$Year==y]
    )
    wt <- alwe$weight
    e <- mass_to_energy("Chinook", wt)
    Susitna_Chinook$biomass[iy] <- sum(Nt*ac*wt)
    Susitna_Chinook$energy[iy] <- sum(Nt*ac*e)
  }
  # Estimate sockeye biomass and energy content
  annual_weight <- readxl::read_xlsx("../Input Data/Salmon.xlsx", "Annual_Weight")
  mean_sockeye_weight <- mean(annual_weight$KG[annual_weight$Species=="Sockeye"])
  Kasilof_Sockeye$biomass <- Kasilof_Sockeye$Nt*mean_sockeye_weight
  Kenai_Sockeye$biomass <- Kenai_Sockeye$Nt*mean_sockeye_weight
  Susitna_Sockeye$biomass <- Susitna_Sockeye$Nt*mean_sockeye_weight
  Kasilof_Sockeye$energy <- Kasilof_Sockeye$Nt*mass_to_energy("Sockeye", mean_sockeye_weight)
  Kenai_Sockeye$energy <- Kenai_Sockeye$Nt*mass_to_energy("Sockeye", mean_sockeye_weight)
  Susitna_Sockeye$energy <- Susitna_Sockeye$Nt*mass_to_energy("Sockeye", mean_sockeye_weight)
}
# See how stuff coorelates
{
  #  Kenai and Susitna River Chinook -- coorelation 1986-2021
  df1 <- data.frame(
    year = Kenai_Chinook$Year,
    biomass = Kenai_Chinook$biomass
  )
  df2 <- data.frame(
    year = Susitna_Chinook$Year,
    biomass = Susitna_Chinook$biomass
  )
  df1 <- df1[df1$year %in% 1986:2021,]
  df2 <- df2[df2$year %in% 1986:2021,]
  cor.test(df1$biomass, df2$biomass)
  # Kenai and Susitna River Sockeye  -- coorelation 2006-2015
  df1 <- data.frame(
    year = Kenai_Sockeye$Year,
    biomass = Kenai_Sockeye$biomass
  )
  df2 <- data.frame(
    year = Kasilof_Sockeye$Year,
    biomass = Kasilof_Sockeye$biomass
  )
  df3 <- data.frame(
    year = Susitna_Sockeye$Year,
    biomass = Susitna_Sockeye$biomass
  )
  df1 <- df1[df1$year %in% 2006:2015,]
  df2 <- df2[df2$year %in% 2006:2015,]
  df3 <- df3[df3$year %in% 2006:2015,]
  cor.test(df1$biomass+df2$biomass, df3$biomass)
  # Susitna Chinook and Sockeye -- coorelation
  df1 <- data.frame(
    year = Susitna_Sockeye$Year,
    biomass = Susitna_Sockeye$biomass
  )
  df2 <- data.frame(
    year = Susitna_Chinook$Year,
    biomass = Susitna_Chinook$biomass
  )
  df1 <- df1[df1$year %in% 2006:2015,]
  df2 <- df2[df2$year %in% 2006:2015,]
  mean(df1$biomass)/1000
  mean(df2$biomass)/1000
  cor.test(df1$biomass, df2$biomass)
  # Kenai Chinook and Sockeye -- coorelation
  df1 <- data.frame(
    year = Kenai_Sockeye$Year,
    biomass = Kenai_Sockeye$biomass
  )
  df2 <- data.frame(
    year = Kenai_Chinook$Year,
    biomass = Kenai_Chinook$biomass
  )
  df1 <- df1[df1$year %in% 1986:2015,]
  df2 <- df2[df2$year %in% 1986:2015,]
  mean(df1$biomass)/1000
  mean(df2$biomass)/1000
  cor.test(df1$biomass, df2$biomass)
}
# Figure
{
  par(mfrow = c(2,2), mar = c(2, 5, 2, 2))
  # panel 1
  plot(Susitna_Chinook$Year, Susitna_Chinook$biomass/1000, ty = "l", ylim = c(0, 1300), xlim = c(1968, 2022), col = "purple", lwd = 2, ylab = "biomass (metric tons)", main = "Northern District Chinook")
  for (m in c(0, 200, 400, 600, 800, 1000, 1200)){
    lines(
      c(1950, 2050),
      c(m, m),
      col = "lightgrey"
    )
  }
  lines(Susitna_Chinook$Year, Susitna_Chinook$biomass/1000, col = "purple", lwd = 2)
  points(Susitna_Chinook$Year, Susitna_Chinook$biomass/1000, pch = 19)
  # panel 2
  plot(Kenai_Chinook$Year, Kenai_Chinook$biomass/1000, ty = "l", ylim = c(0,650), xlim = c(1968, 2022), col = "purple", lwd = 2, ylab = "biomass (metric tons)", main = "Central District Chinook")
  for (m in c(0, 100, 200, 300, 400, 500, 600)){
    lines(
      c(1950, 2050),
      c(m, m),
      col = "lightgrey"
    )
  }
  lines(Kenai_Chinook$Year, Kenai_Chinook$biomass/1000, col = "purple", lwd = 2)
  points(Kenai_Chinook$Year, Kenai_Chinook$biomass/1000, pch = 19)
  # panel 3
  plot(Susitna_Sockeye$Year, Susitna_Sockeye$biomass/1000, ty = "l", ylim = c(0, 1700), xlim = c(1968, 2022), col = "blue", lwd = 2, ylab = "biomass (metric tons)", main = "Northern District Sockeye")
  for (m in c(0, 250, 500, 750, 1000, 1250, 1500, 1750)){
    lines(
      c(1950, 2050),
      c(m, m),
      col = "lightgrey"
    )
  }
  lines(Susitna_Sockeye$Year, Susitna_Sockeye$biomass/1000, col = "purple", lwd = 2)
  points(Susitna_Sockeye$Year, Susitna_Sockeye$biomass/1000, pch = 19)
  # panel 4
  plot(Kenai_Sockeye$Year, (Kenai_Sockeye$biomass + Kasilof_Sockeye$biomass[-49])/1000, ylim = c(0, 35000), ty = "l",  xlim = c(1968, 2022), col = "blue", lwd = 2, ylab = "biomass (metric tons)", main = "Central District Sockeye")
  for (m in c(0, 5000, 10000, 15000, 20000, 25000, 30000, 35000)){
    lines(
      c(1950, 2050),
      c(m, m),
      col = "lightgrey"
    )
  }
  lines(Kenai_Sockeye$Year, (Kenai_Sockeye$biomass + Kasilof_Sockeye$biomass[-49])/1000, col = "purple", lwd = 2)
  points(Kenai_Sockeye$Year, (Kenai_Sockeye$biomass + Kasilof_Sockeye$biomass[-49])/1000, pch = 19)
}
# Figure
{
  par(mfrow = c(1,1))
  plot(
    Kenai_Chinook$Year,
    Kenai_Chinook$biomass/1000,
    xlim = c(1968, 2022),
    ylim = c(0, 55000),
    ty = "l",
    col = "red",
    lwd = 2,
    ylab = "biomass (metric tons)",
    main = "Prey Biomass in Major Rivers"
  )
  for (m in c(0, 10000, 20000, 30000, 40000, 50000)){
    lines(
      c(1950, 2050),
      c(m, m),
      col = "lightgrey"
    )
  }
  lines(
    Susitna_Chinook$Year,
    Susitna_Chinook$biomass/1000,
    col = "blue",
    lwd = 2
  )
  lines(
    Kasilof_Sockeye$Year,
    Kasilof_Sockeye$biomass/1000,
    col = "green",
    lwd = 3,
    lty = 3
  )
  lines(
    Kenai_Sockeye$Year,
    Kenai_Sockeye$biomass/1000,
    col = "red",
    lwd = 3,
    lty = 3
  )
  lines(
    Susitna_Sockeye$Year,
    Susitna_Sockeye$biomass/1000,
    col = "blue",
    lwd = 3,
    lty = 3
  )
  points(
    2016,
    48000,
    pch = 19
  )
  legend(
    x = "topleft",
    legend = c(
      "Chinook",
      "Sockeye",
      "Kasilof",
      "Kenai",
      "Susitna"
    ),
    col = c("black", "black", "green", "red", "blue"),
    lty = c(1, 3, 1, 1, 1),
    lwd = c(2, 3, 2, 2, 2)
  )
  legend(
    x = "topright",
    legend = "Eulachon",
    col = "black",
    pch = 19
  )
}
# Compile and write to disk
writexl::write_xlsx(
  rbind(
    data.frame(
      Year = Kenai_Chinook$Year,
      Species = "Chinook",
      River = "Kenai",
      Run_Size = Kenai_Chinook$Nt,
      Biomass_MT = round(Kenai_Chinook$biomass/1000,0)
    ),
    data.frame(
      Year = Susitna_Chinook$Year,
      Species = "Chinook",
      River = "Susitna",
      Run_Size = Susitna_Chinook$Nt,
      Biomass_MT = round(Susitna_Chinook$biomass/1000,0)
    ),
    data.frame(
      Year = Kasilof_Sockeye$Year,
      Species = "Sockeye",
      River = "Kasilof",
      Run_Size = Kasilof_Sockeye$Nt,
      Biomass_MT = round(Kasilof_Sockeye$biomass/1000,0)
    ),
    data.frame(
      Year = Kenai_Sockeye$Year,
      Species = "Sockeye",
      River = "Kenai",
      Run_Size = Kenai_Sockeye$Nt,
      Biomass_MT = round(Kenai_Sockeye$biomass/1000,0)
    ),
    data.frame(
      Year = Susitna_Sockeye$Year,
      Species = "Sockeye",
      River = "Susitna",
      Run_Size = Susitna_Sockeye$Nt,
      Biomass_MT = round(Susitna_Sockeye$biomass/1000,0)
    )
  ),
  path = "../Output Data/Salmon_Availability.xlsx"
)
# Figure 4
{
  input <- readxl::read_xlsx("../Input Data/Salmon.xlsx", sheet = "Annual_Weight")
  chinook <- input[input$Species=="Chinook",]
  chum <- input[input$Species=="Chum",]
  coho <- input[input$Species=="Coho",]
  pink <- input[input$Species=="Pink",]
  sockeye <- input[input$Species=="Sockeye",]
  par(mfrow = c(1,1))
  plot(chinook$Year, chinook$LBS, ty = "l", ylim = c(0,35), col = "blue", lwd = 2, ylab = "mass (in pounds)", xlab = "")
  for (m in c(0, 5, 10, 15, 20, 25, 30, 35)){
    lines(
      c(1950, 2050),
      c(m, m),
      col = "lightgrey"
    )
  }
  lines(chinook$Year, chinook$LBS,col = "blue", lwd = 2)
  lines(chum$Year, chum$LBS, ty = "l", col = "green", lwd = 2)
  lines(coho$Year, coho$LBS, ty = "l", col = "grey", lwd = 2)
  lines(pink$Year, pink$LBS, ty = "l", col = "pink", lwd = 2)
  lines(sockeye$Year, sockeye$LBS, ty = "l", col = "red", lwd = 2)
  legend(
    x = "topright",
    legend = c(
      "Chinook",
      "Chum",
      "Coho",
      "Pink",
      "Sockeye"
    ),
    col = c("blue", "green", "grey", "pink", "red"),
    lty = 1,
    lwd = 2
  )
}
