rm(list=ls())
source("01_functions.R")
library(statmod) ## for integrating ordinary differential equations
require(tidyverse) ## for efficient data manipulation & plotting
library(cowplot) ## for arranging plots in a grid
library(dplyr)
library(readr)
library(beepr)
library(deSolve)



#example dynamics of three species plant-pollinator system 

#adjacency matrix
# g<-matrix(1,nrow=30,ncol=70)
# 
# Pspecies<-dim(g)[1]
# Aspecies<- dim(g)[2]

# #parameters for modelling: initialisation of the model
# sigma<-c(0.01,0.01, 0.01) #three species trait variance 
# Na<-1 #initial abundance abundance
# Np<-c(1,1) #initial abundance plants
# 
# muA<- 35  #initial mean phenotypic optimum trait values
# muP<- c(37,33)   #initial mean phenotypic optimum trait values


#parameters below are taken from another paper Akesson et al 2021 Nat Comms.
# bw  <- 0.1  #modulate tradeoff between maximum growth and tolerance range
# aw  <- 0.1  #
# gi <- 0.9   #
# ki <-0.1    #mortality rate
# w<- 10#mutualism interaction width
# Temp<-30  #seq(Tmin, Tmax, by=(Tmax-Tmin)/50)
# h2<-c(0.5,0.05,0.5)  #heritability of trait variance 
# Amatrix <- mat.comp(g)$Amatrix  #competition matrix , aij, for animals
# Pmatrix <- mat.comp(g)$Pmatrix  #competition matrix , aij, for plants
# mut.strength=1 #average mutualistic strength

#########################################
#         RStudio commands

# Ctrl + shift + c --> block comment

#########################################

# g<-matrix(1,nrow=3,ncol=10)
# g[1,3] <- 0
# g[1,5] <- 0
# g[1,7] <- 0
# g[1,8] <- 0
# g[2,1] <- 0
# g[2,2] <- 0
# g[2,4] <- 0
# g[2,6] <- 0
# g[3,2] <- 0
# g[3,8] <- 0

# Pspecies<-dim(g)[1]
# Aspecies<- dim(g)[2]
# 
# lowerTemp = 10
# upperTemp = 25

#sigma <- rep(c(0.01, 0.1), times = c(Pspecies, Aspecies))

#Np <- rep(1, Pspecies)
#Na <- rep(1, Aspecies)

#optimumTemps <- runif(Pspecies + Aspecies, lowerTemp, upperTemp)
#muP <- optimumTemps[1:Pspecies]
#muA <- optimumTemps[(Pspecies + 1) : (Pspecies + Aspecies)]

#bw  <- 0.1  #modulate tradeoff between maximum growth and tolerance range
#aw  <- 0.1  #
#gi <- 0.9   #
#ki <-0.1    #mortality rate
#w<- 10#mutualism interaction width
#Temp<-30  #seq(Tmin, Tmax, by=(Tmax-Tmin)/50)
#h2 <- rep(c(0.05, 0.5), times = c(Pspecies, Aspecies))
#Amatrix <- mat.comp(g)$Amatrix  #competition matrix , aij, for animals
#Pmatrix <- mat.comp(g)$Pmatrix  #competition matrix , aij, for plants
#mut.strength=1 #average mutualistic strength y_0
#.........................................

#params<-list(matrix=g,bw=bw,aw=aw,h2=h2,w=w,Amatrix=Amatrix,Pmatrix=Pmatrix,
#             gi=gi,ki=ki,Temp=Temp, sigma=sigma,
#             mut.strength=mut.strength)

#ic<-c(Na,Np,muA,muP) ## initial conditions coerced into a vector

#tmax <- 1e5 ## time to integrate equations for
#sol<-ode(func=eqs, y=ic, parms=params, times=seq(0, tmax, by=tmax/1000)) %>% ## solve ODEs
#  organize_results(params) #%>% ## put results in tidy table %>%## keep data for only the final point in time
  #plot_all()

#View(sol)

#temp <- sol %>% filter(time == 100000)
#View(temp)


#extract animal biomass
#animal_biomass <- temp %>% filter(type == "N")

#sum over final animal biomass
#sum(animal_biomass[,4])

#####################################################################################
#                                                                                   #
#         new                                                                       # 
#                                                                                   #
#####################################################################################
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

Pspecies<-dim(g)[1]
Aspecies<- dim(g)[2]

lowerTemp = 10
upperTemp = 25

Np <- rep(1, Pspecies)
Na <- rep(1, Aspecies)

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
              replicates= (1:3)) %>% 
  as_tibble() %>% 
  mutate(mean_animal_trait=0,
         mean_plant_trait=0,
         mean_community_trait=0,
         animal_biomass=0,
         plant_biomass=0,
         community_biomass=0,
         animal_richness=0,
         plant_richness=0,
         community_richness=0)

view(dat)

for(i in 1:nrow(dat)) {
   
  if(dat$trait.variation[i] == "high") {
    sigma <- runif( (Aspecies+Pspecies), 0.01,0.05)
  } else {
    sigma <- runif( (Aspecies+Pspecies), 0.001,0.005)
  }
  
  optimumTemps <- runif(Pspecies + Aspecies, lowerTemp, upperTemp)
  muP <- optimumTemps[1:Pspecies]
  muA <- optimumTemps[(Pspecies + 1) : (Pspecies + Aspecies)]
  
  Temp<-dat$temperature[i]  #seq(Tmin, Tmax, by=(Tmax-Tmin)/50)
  h2 <- dat$h2[i]  #rep(c(0.05, 0.5), times = c(Pspecies, Aspecies))

  #.........................................
  
  params<-list(matrix=g,bw=bw,aw=aw,h2=h2,w=w,Amatrix=Amatrix,Pmatrix=Pmatrix,
               gi=gi,ki=ki,Temp=Temp, sigma=sigma,
               mut.strength=1)
  
  ic<-c(Na,Np,muA,muP) ## initial conditions coerced into a vector
  
  tmax <- 1e5 ## time to integrate equations for
  sol<-ode(func=eqs, y=ic, parms=params, times=seq(0, tmax, by=tmax/1000)) %>% ## solve ODEs
    organize_results(params) #%>% ## put results in tidy table %>%## keep data for only the final point in time

  
  temp <- sol %>% filter(time == 100000)
  
  #view(temp)
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

#     test

#----------------------------
#view(sol_evolution)

#extracting replicates for given temperature and trait variation

# temp15 = dat_02 %>% filter(temperature == 15, trait_variation == "low")
# view(temp15)

#----------------------------

#     saving data

#----------------------------
write.csv(sol, "230518_sol_evolution.csv")
write.csv(dat,"230518_dat_01.csv")

#----------------------------

#     plotting data

#----------------------------

dat_01 = read.csv("230518_dat_01.csv", header = T)
sol_evolution = read.csv("230518_sol_evolution.csv", header = T)
view(dat_01)

#template for plotting-------------------

#ggplot(data = <DATA>) + 
#  <GEOM_FUNCTION>(mapping = aes(<MAPPINGS>))

#---------------------------


ggplot(data = dat_01) +
  geom_point(mapping = aes(x = temperature, y = community_richness, color = h2))


ggplot(data = dat_01) + 
  geom_point(mapping = aes(x = temperature, y = community_richness, color = trait.variation)) #+
  #geom_abline(slope=1)

ggplot(dat_01) +
   geom_point(aes(x = temperature, y = community_richness, colour = h2)) +
   facet_wrap(~trait.variation) #+
   #geom_abline(slope=1, colour = "gray")
   


