rm(list=ls())
source("01_functions.R")
library(statmod) ## for integrating ordinary differential equations
require(tidyverse) ## for efficient data manipulation & plotting
library(cowplot) ## for arranging plots in a grid
library(dplyr)
library(readr)
library(beepr)
library(deSolve)

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
#mut.strength = 1 #average mutualistic strength

dat<-expand.grid(h2 = c(0, 0.4),
                 mutualistic_strength = seq(0.5,5, 0.25),
                 temperature = c(15,20,28,40),
                 trait.variation= c("high"),
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
  mut.strength = dat$mutualistic_strength[i]
  
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

#----------------------------

#     saving data

#----------------------------

write.csv(dat,"230530_dat_01_mut_str.csv")

#----------------------------

#     loading data

#----------------------------

dat_01 = read.csv("230530_dat_01_mut_str.csv", header = T)

view(dat_01)

#----------------------------

#     plotting by temperature

#----------------------------

ggplot(data = dat_01 %>% filter(temperature == 15)) +
  geom_point(mapping = aes(x = mutualistic_strength, y = community_richness, color = factor(h2)))

ggplot(data = dat_01 %>% filter(temperature == 20)) +
  geom_point(mapping = aes(x = mutualistic_strength, y = community_richness, color = factor(h2)))

ggplot(data = dat_01 %>% filter(temperature == 28)) +
  geom_point(mapping = aes(x = mutualistic_strength, y = community_richness, color = factor(h2)))

ggplot(data = dat_01 %>% filter(temperature == 40)) +
  geom_point(mapping = aes(x = mutualistic_strength, y = community_richness, color = factor(h2)))

#----------------------------

#     plotting data, general

#----------------------------

#...template for plotting... 

#ggplot(data = <DATA>) + 
#  <GEOM_FUNCTION>(mapping = aes(<MAPPINGS>))

#--------------------------- 


ggplot(data = dat_01) +
  geom_point(mapping = aes(x = temperature, y = community_richness, color = (mutualistic_strength)))

ggplot(data = dat_01) +
  geom_point(mapping = aes(x = temperature, y = plant_richness, color = factor(mutualistic_strength)))




