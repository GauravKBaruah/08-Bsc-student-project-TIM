rm(list=ls())
source("01_functions.R")
library(statmod) ## for integrating ordinary differential equations
require(tidyverse) ## for efficient data manipulation & plotting
library(cowplot) ## for arranging plots in a grid
library(dplyr)
library(readr)
library(beepr)
library(deSolve)

#########################################
#         RStudio commands

# Ctrl + shift + c --> block comment

#########################################

g <- read.csv("web-of-life_2023-05-04_172346/M_PL_069_01.txt", header = T, sep = ",")
# if mutualistic interaction => setting strength to '1'
for (i in 1:nrow(g)) {        
  for (j in 1:ncol(g)) {
    if (g[i,j] > 1) {
      g[i,j] <- 1
    }
  }
}

view(g)

Aspecies<- dim(g)[2]
Pspecies<-dim(g)[1]

lowerTemp = 10
upperTemp = 25

Na <- rep(1, Aspecies)
Np <- rep(1, Pspecies)

bw  <- 0.1  #modulate tradeoff between maximum growth and tolerance range
aw  <- 0.1  #
gi <- 0.9   #
ki <-0.1    #mortality rate
w<- 5      #mutualism interaction width
Amatrix <- mat.comp(g)$Amatrix  #competition matrix , aij, for animals
Pmatrix <- mat.comp(g)$Pmatrix  #competition matrix , aij, for plants
mut.strength = 1 #average mutualistic strength

dat<-expand.grid(h2 = c(0, 0.4),
              temperature = seq(15,45, 1),
              trait.variation= c("high"),#, "low"),
              replicates= (1:1)) %>% 
  as_tibble() %>% 
  mutate(mean_animal_trait=0,
         mean_plant_trait=0,
         mean_community_trait=0,
         animal_biomass=0,
         plant_biomass=0,
         community_biomass=0,
         animal_richness=0,
         plant_richness=0,
         community_richness=0,
         temp_of_collapse=0,
         abruptness=0)

view(dat)

#----------------------------

#     trait values --> table for saving

#----------------------------

trait_values <- expand.grid(species = seq(1, Aspecies + Pspecies),
                            h2 = c(0,0.4),
                            temperature = seq(15, 45, 1)) %>% 
                as_tibble() %>% 
                mutate(trait_value = 0)

view(trait_values)

#----------------------------

for(i in 1:nrow(dat)) {
   
  if(dat$trait.variation[i] == "high") {
    sigma <- runif( (Aspecies+Pspecies), 0.01,0.05)
  } else {
    sigma <- runif( (Aspecies+Pspecies), 0.001,0.005)
  }
  
  optimumTemps <- runif(Aspecies + Pspecies, lowerTemp, upperTemp)
  muA <- optimumTemps[1:Aspecies]
  muP <- optimumTemps[(Aspecies + 1) : (Aspecies + Pspecies)]
  
  Temp<-dat$temperature[i]  
  
  h2 <- dat$h2[i]  

  #.........................................
  
  params<-list(matrix=g,bw=bw,aw=aw,h2=h2,w=w,Amatrix=Amatrix,Pmatrix=Pmatrix,
               gi=gi,ki=ki,Temp=Temp, sigma=sigma,
               mut.strength=mut.strength)
  
  ic<-c(Na,Np,muA,muP) ## initial conditions coerced into a vector
  
  tmax <- 1e5 ## time to integrate equations for
  sol<-ode(func=eqs, y=ic, parms=params, times=seq(0, tmax, by=tmax/1000)) %>% ## solve ODEs
    organize_results(params) #%>% ## put results in tidy table %>%## keep data for only the final point in time
  
  temp <- sol %>% filter(time == 100000)
  
  #extract animal & plant biomass 
  end_animal <- temp %>% filter(type == "N")
  end_plant <- temp %>% filter(type == "P")
  
  #sum over final animal & plant biomass
  dat$animal_biomass[i]<-sum(end_animal[,4])
  dat$plant_biomass[i]<-sum(end_plant[,4])
  
  dat$community_biomass[i] <- dat$animal_biomass[i] + dat$plant_biomass[i]
  
  #extract animal & plant mean trait 
  animal_trait <- temp %>% filter(type == "ma")
  plant_trait <- temp %>% filter(type == "mp")  
  
  
  #save final trait values
  
  if(h2 == 0) {
    for(a in 1:Aspecies) {
      trait_values$trait_value[2*(Temp - 15)*(Aspecies + Pspecies) + a] <- animal_trait$v[a]
    }
    for(p in 1:Pspecies) {
      trait_values$trait_value[2*(Temp - 15)*(Aspecies + Pspecies) + Aspecies + p] <- plant_trait$v[p]
    }
  }
  else {
    for(a in 1:Aspecies) {
      trait_values$trait_value[(2*(Temp - 15) + 1)*(Aspecies + Pspecies) + a] <- animal_trait$v[a]
    }
    for(p in 1:Pspecies) {
      trait_values$trait_value[(2*(Temp - 15) + 1)*(Aspecies + Pspecies) + Aspecies + p] <- plant_trait$v[p]
    }
  }
  
  #calculate animal & plant mean trait
  dat$mean_animal_trait[i] <- mean(animal_trait$v)
  dat$mean_plant_trait[i] <- mean(plant_trait$v)
  
  dat$mean_community_trait[i] <- mean((temp %>% filter (type == "ma" | type == "mp")) $v)
  
  #extract end community-richness

  animal_richness = length(which(end_animal$v > 1e-2))
  
  plant_richness = length(which(end_plant$v > 1e-2))
  
  dat$animal_richness[i]=animal_richness
  dat$plant_richness[i]=plant_richness
  dat$community_richness[i]=animal_richness + plant_richness
    
}

view(dat)
view(trait_values)

#----------------------------

#     saving data

#----------------------------
write.csv(sol, "230518_sol_evolution.csv")
write.csv(dat,"230530_dat_01.csv")
write.csv(trait_values, "230530_trait_values01.csv")

#----------------------------

#     loading data

#----------------------------

dat_01 = read.csv("230530_dat_01.csv", header = T)
trait_vals_01 = read.csv("230530_trait_values01.csv", header = T)

view(dat_01)
view(trait_vals_01)

#----------------------------

#     community-collapse & abruptness

#----------------------------

# collapse, toc = temperature of collapse
toc_noEvo = 0  
toc_evo = 0

collapse_noEvo = FALSE
collapse_evo = FALSE

for(i in 1:nrow(dat_01)) {
  
  if(dat_01$community_richness[i] == 0) {
    if(dat_01$h2[i] == 0 & !collapse_noEvo) {
      toc_noEvo = dat_01$temperature[i]
      collapse_noEvo = TRUE
    }
    else if(dat_01$h2[i] > 0 & !collapse_evo) {
      toc_evo = dat_01$temperature[i]
      collapse_evo = TRUE
    }
  }
  
}

print(toc_noEvo)
print(toc_evo)

#adding toc-value to dat_01 table --> is this needed?

# for(i in 1:nrow(dat_01)) {
#   if(dat_01$h2[i] == 0) {
#     dat_01$temp_of_collapse[i] = toc_noEvo
#   }
#   else {
#     dat_01$temp_of_collapse[i] = toc_evo
#   }
# 
# }

view(dat_01)


# abruptness --> (richness at temp_of_collapse - 1) - (richness at temp_of_collapse) <-- shouldn't this be 0?

abruptness_noEvo = (dat_01 %>% filter((temperature == toc_noEvo - 1) & h2 == 0))$community_richness -
                   (dat_01 %>% filter((temperature == toc_noEvo) & h2 == 0))$community_richness
abruptness_evo = (dat_01 %>% filter((temperature == toc_evo - 1) & h2 > 0))$community_richness -
                 (dat_01 %>% filter((temperature == toc_evo ) & h2 > 0))$community_richness 

for(i in 1:nrow(dat_01)) {
  if(dat_01$h2[i] == 0) {
    dat_01$abruptness[i] = abruptness_noEvo
  }
  else {
    dat_01$abruptness[i] = abruptness_evo
  }
}

view(dat_01)

#----------------------------

#     hysterisis

#----------------------------

dat<-expand.grid(h2 = c(0, 0.4),
                 temperature = seq(45,15, -1),
                 trait.variation= c("high"),#, "low"),
                 replicates= (1:1)) %>%
  as_tibble() %>%
  mutate(mean_animal_trait=0,
         mean_plant_trait=0,
         mean_community_trait=0,
         animal_biomass=0,
         plant_biomass=0,
         community_biomass=0,
         animal_richness=0,
         plant_richness=0,
         community_richness=0,
         temp_of_recovery=0,
         abruptness=0)

view(dat)

#trait_values for temperature with highest community richness
mu_temp <- (trait_vals_01 %>% filter(temperature == 18))$trait_value 
view(mu_temp)

for(i in 1:nrow(dat)) {

  if(dat$trait.variation[i] == "high") {
    sigma <- runif( (Aspecies+Pspecies), 0.01,0.05)
  } else {
    sigma <- runif( (Aspecies+Pspecies), 0.001,0.005)
  }

  #adjusting traits and biomasses
  Na <- runif(Aspecies, 0,0.005)
  Np <- runif(Pspecies, 0,0.005)

  h2 <- dat$h2[i]  
  
  if(h2 == 0) {
    muA <- mu_temp[1:Aspecies]
    muP <- mu_temp[(Aspecies + 1):(Aspecies + Pspecies)]
  }
  else {
    muA <- mu_temp[(Aspecies + Pspecies + 1):(2*Aspecies + Pspecies)]
    muP <- mu_temp[(2*Aspecies + Pspecies + 1):(2*(Aspecies + Pspecies))]
    
  }

  Temp<-dat$temperature[i]  

  #.........................................

  params<-list(matrix=g,bw=bw,aw=aw,h2=h2,w=w,Amatrix=Amatrix,Pmatrix=Pmatrix,
               gi=gi,ki=ki,Temp=Temp, sigma=sigma,
               mut.strength=mut.strength)

  ic<-c(Na,Np,muA,muP) ## initial conditions coerced into a vector

  tmax <- 1e5 ## time to integrate equations for
  sol<-ode(func=eqs, y=ic, parms=params, times=seq(0, tmax, by=tmax/1000)) %>% ## solve ODEs
    organize_results(params) #%>% ## put results in tidy table %>%## keep data for only the final point in time


  temp <- sol %>% filter(time == 100000)

  #extract animal & plant biomass
  end_animal <- temp %>% filter(type == "N")
  end_plant <- temp %>% filter(type == "P")

  #sum over final animal & plant biomass
  dat$animal_biomass[i]<-sum(end_animal[,4])
  dat$plant_biomass[i]<-sum(end_plant[,4])

  dat$community_biomass[i] <- dat$animal_biomass[i] + dat$plant_biomass[i]

  #extract animal & plant mean trait
  animal_mean_trait <- temp %>% filter(type == "ma")
  plant_mean_trait <- temp %>% filter(type == "mp")

  #calculate animal & plant mean trait
  dat$mean_animal_trait[i] <- mean(animal_mean_trait$v)
  dat$mean_plant_trait[i] <- mean(plant_mean_trait$v)

  dat$mean_community_trait[i] <- mean((temp %>% filter (type == "ma" | type == "mp")) $v)

  #extract end community-richness

  animal_richness = length(which(end_animal$v > 1e-2))

  plant_richness = length(which(end_plant$v > 1e-2))

  dat$animal_richness[i]=animal_richness
  dat$plant_richness[i]=plant_richness
  dat$community_richness[i]=animal_richness + plant_richness

}
view(dat)

#----------------------------

#     saving data hysterisis

#----------------------------

write.csv(dat,"230604_dat_01_hysterisis.csv")

#----------------------------

#     loading data hysterisis

#----------------------------

dat_01_hysterisis = read.csv("230604_dat_01_hysterisis.csv", header = T)

view(dat_01_hysterisis)

#----------------------------

#     plotting data

#----------------------------

#...template for plotting... 

#ggplot(data = <DATA>) + 
#  <GEOM_FUNCTION>(mapping = aes(<MAPPINGS>))

#--------------------------- normal


ggplot(data = dat_01) +
  geom_point(mapping = aes(x = temperature, y = community_richness, color = factor(h2)))

ggplot(data = dat_01) +
  geom_point(mapping = aes(x = temperature, y = plant_richness, color = factor(h2)))

ggplot(data = dat_01) +
  geom_point(mapping = aes(x = temperature, y = animal_richness, color = factor(h2)))

#--------------------------- hysterisis

ggplot(data = dat_01_hysterisis) +
  geom_point(mapping = aes(x = temperature, y = community_richness, color = factor(h2)))

ggplot(data = dat_01_hysterisis) +
  geom_point(mapping = aes(x = temperature, y = plant_richness, color = factor(h2)))

ggplot(data = dat_01_hysterisis) +
  geom_point(mapping = aes(x = temperature, y = animal_richness, color = factor(h2)))

#--------------------------- 

ggplot(data = dat_01) + 
  geom_point(mapping = aes(x = temperature, y = community_richness, color = factor(trait.variation))) #+
  #geom_abline(slope=1)

ggplot(dat_01) +
   geom_point(aes(x = temperature, y = community_richness, colour = h2)) +
   facet_wrap(~trait.variation) #+
   #geom_abline(slope=1, colour = "gray")
   