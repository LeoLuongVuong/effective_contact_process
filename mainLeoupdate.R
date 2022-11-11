
#Input Arguments
n<-1000                   # Population size
prop.immune <-0.8         # Proportion of immune individuals
rho <-0.8                 # Probability of symptomatic infection
q <-0.5                   # Transmission potential
alpha.as<- 0.1            # Relative infectiousness asymptomatic carriers
vacc.eff<- 0.4            # Vaccine effectiveness (against infectiousness and susceptibility to infection)
testing.prob<-1           # Probability that a symptomatic individual is tested
test.sens<-0.9            # Sensitivity of the test
test.delay<-1             # Delay from taking the test to the test result
contact.reduction<-0.1    # Reduction of contact rate
nSeeds<-3                 # Number initial infected individuals
lambda<-5                 # Number of daily contacts

#running simulations
source("scrLeo.R")
nSim<-20
set.seed(131714)
epi.outbreak<-list()
for (i in 1:nSim){
  print(i)
  epi.outbreak[[i]]<-sim.ekp(n=n,prop.immune = prop.immune,rho = rho,q=q,alpha.as = alpha.as,vacc.eff=vacc.eff,testing.prob = testing.prob,test.sens = test.sens,test.delay = test.delay,contact.reduction = contact.reduction,nSeeds = nSeeds, lambda = lambda)
}




# Plot Epidemiological quantities
#Look at the evolution for a single epidemic
setwd("D:/Hoc/Uantwerp/Infectious Disease Modelling Internship/scripts/Simulations and Rhistory")
single.ep<-epi.outbreak[[sample(1:nSim,1)]]

tiff('Prevalence_single_vacceff0.4.tiff', units="in", width=5, height=4, res=100, compression = 'lzw')
par(mai=c(0.8,0.8,0.1,0.1))
plot(single.ep$epi.evo$Days, single.ep$epi.evo$Prevalence, type="l", xlab = "Days", ylab = "Prevalence")
dev.off()

tiff('Incidence_vacceff0.4.tiff', units="in", width=5, height=4, res=100, compression = 'lzw')
par(mai=c(0.8,0.8,0.1,0.1))
plot(single.ep$epi.evo$Days, single.ep$epi.evo$Incidence, col="red", xlab = "Days", ylab = " Incidence")
dev.off()

tiff('Eff_Rep_Num_vacceff0.4.tiff', units="in", width=5, height=4, res=100, compression = 'lzw')
par(mai=c(0.8,0.8,0.1,0.1))
plot(single.ep$epi.evo$Days, single.ep$epi.evo$Rt, col="blue", ylab = "Effective Reproductive Number", xlab = "Days" )
dev.off()

#look at the evolution of different epidemics
n.epidemics<-20 #Number of epidemics you want to consider
epidemics<-sample(1:nSim,n.epidemics) #Make sure you are selecting a number of epidemics smaller than the number of simulations


tiff('Prevalence_multiple_vacceff0.4.tiff', units="in", width=5, height=4, res=100, compression = 'lzw')
par(mai=c(0.8,0.8,0.1,0.1))
plot(epi.outbreak[[epidemics[1]]]$epi.evo$Days, epi.outbreak[[epidemics[1]]]$epi.evo$Prevalence, type="l", xlab = "Days", ylab = "Prevalence", xlim = c(0,110),ylim = c(0,25)) # you need to adjust the x and y axis, to see all the epidemics. To do so adjust xlim and ylim
for (i in 1:n.epidemics){
  lines(epi.outbreak[[epidemics[i]]]$epi.evo$Days, epi.outbreak[[epidemics[i]]]$epi.evo$Prevalence, col=i)
}
dev.off()



#Look at summary measures among different epidemics
FinalSize<-NULL
PeakIncidence<-NULL
PeakPrevalence<-NULL

for (i in 1:nSim){
  FinalSize<-c(FinalSize,epi.outbreak[[i]]$FinalSize$FinalSize)
  PeakIncidence<-c(PeakIncidence, epi.outbreak[[i]]$PeakIncidence$PeakIncidence)
  PeakPrevalence<-c(PeakPrevalence,epi.outbreak[[i]]$PeakPrevalence$PeakPrevalence)
}

tiff('Finalsize_vacceff0.4.tiff', units="in", width=5, height=4, res=100, compression = 'lzw')
par(mai=c(0.8,0.8,0.1,0.1))
boxplot(FinalSize, ylab="Final Size")
dev.off()

tiff('PeakIncidence_vacceff0.4.tiff', units="in", width=5, height=4, res=100, compression = 'lzw')
par(mai=c(0.8,0.8,0.1,0.1))
boxplot(PeakIncidence, ylab="Peak Incidence")
dev.off()

tiff('PeakPrevalence_vacceff0.4.tiff', units="in", width=5, height=4, res=100, compression = 'lzw')
par(mai=c(0.8,0.8,0.1,0.1))
boxplot(PeakPrevalence, ylab="Peak Prevalence")
dev.off()

# Save Simulations

name<-paste("EpiOutbreak", "_N",n,"_nSeeds",nSeeds,"_PropImm",prop.immune,"_rho",rho,"_q",q,"_alpha",alpha.as,"_vacc.eff",vacc.eff,"_testingProb",testing.prob,"_testSens",test.sens,"_testdelay",test.delay,"_contact reduction",contact.reduction,sep = "")
save(epi.outbreak, file = paste(name,".RData",sep = ""))
setwd("D:/Hoc/Uantwerp/Infectious Disease Modelling Internship/scripts/effective_contact_process")

