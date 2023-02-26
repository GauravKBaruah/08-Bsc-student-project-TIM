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
g<-matrix(1,nrow=2,ncol=1)

Aspecies<- dim(g)[2]
Pspecies<-dim(g)[1]


#parameters for modelling: intialisation of the model
sigma<-c(0.01,0.01, 0.01) #three species trait variance 
Na<-1 #intial abundance abundance
Np<-c(1,1) #initial abundance plants

muA<- 35  #initial mean phenotypic optimum trait values
muP<- c(37,33)   #intial mean phenotypic optimum trait values

#parameters below are taken from another paper Akesson et al 2021 Nat Comms.
bw  <- 0.1
aw  <- 0.1
gi <- 0.9
ki <-0.1 #mortality rate
w<- 5#mutualism interaction width
Temp<-30  #seq(Tmin, Tmax, by=(Tmax-Tmin)/50)
h2<-c(0.5,0.05,0.5)  #heritability of trait variance 
Amatrix <- mat.comp(g)$Amatrix  #competition matrix , aij, for animals
Pmatrix <- mat.comp(g)$Pmatrix  #competition matrix , aij, for plants
mut.strength=1 #average mutualistic strength

params<-list(matrix=g,bw=bw,aw=aw,h2=h2,w=w,Amatrix=Amatrix,Pmatrix=Pmatrix,
             gi=gi,ki=ki,Temp=Temp, sigma=sigma,
             mut.strength=1)

ic<-c(Na,Np,muA,muP) ## initial conditions coerced into a vector

tmax <- 1e5 ## time to integrate equations fors
sol<-ode(func=eqs, y=ic, parms=params, times=seq(0, tmax, by=tmax/1000)) %>% ## solve ODEs
  organize_results(params) %>% ## put results in tidy table %>%## keep data for only the final point in time
  plot_all()

sol




 