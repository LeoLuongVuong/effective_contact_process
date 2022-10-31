###################################################
# Functions used to assign:
# - incubation period 
# - infectious period
# - infectivity Measure
###################################################

# The following functions are set to describe the spreading of SARS-CoV-2

#infectious period length

infectious.period.length<-function(){
return(rgamma(1,shape = 32.14,scale = 0.47))#better to move to Gamma, e.g. Gamma with shape 26.01 scale 0.392 -> mean 10.2 and 2 sd
}

# infectivity measures for vaccinated and unvaccinated individuals

nCov.InfMeasure.unvacc<-function(t){
  vload.comp<-NULL
  for (j in 1:length(t)){
    if (t[j]>=1.7 & t[j]<=5.2){
      vload.comp<-c(vload.comp,(5.514286*t[j]-49.37429+40)*0.00942063)
    }
    if (t[j]>5.2 & t[j]<=12.7){
      vload.comp<-c(vload.comp,(-2.573333*t[j]-7.318667+40)*0.00942063)
    }
    if (t[j]<1.7 | t[j]>12.7){
      vload.comp<-c(vload.comp,0)
    }
  }
  return(vload.comp)
}

nCov.InfMeasure.vacc<-function(t){
  vload.comp<-NULL
  for (j in 1:length(t)){
    if (t[j]>=2 & t[j]<=5.2){
      vload.comp<-c(vload.comp,(6.09375*t[j]-52.1875+40)*0.01243202)
    }
    if (t[j]>5.2 & t[j]<=10.7){
      vload.comp<-c(vload.comp,(-3.545455*t[j]-2.063636+40)*0.01243202)
    }
    if (t[j]<2 | t[j]>10.7){
      vload.comp<-c(vload.comp,0)
    }
  }
  return(vload.comp)
}

# We compute the day of symptom onset as the peak of the infectivity measure fpr vaccinated and unvaccinated individuals

incubation.period.vacc<-function(lengthIP){ 
  time.points<-seq(0,lengthIP,0.01)
  infectiousnessmeasure.values<-nCov.InfMeasure.vacc(time.points) #viral-loads
  time.max<-time.points[which(infectiousnessmeasure.values==max(infectiousnessmeasure.values))] 
  return(time.max) 
}

incubation.period.unvacc<-function(lengthIP){ 
  time.points<-seq(0,lengthIP,0.01)
  infectiousnessmeasure.values<-nCov.InfMeasure.unvacc(time.points) #viral-loads
  time.max<-time.points[which(infectiousnessmeasure.values==max(infectiousnessmeasure.values))] 
  return(time.max) 
}

#########################3
# Simulating Function
#########################3

sim.ekp<-function(n,prop.immune, rho,q, alpha.as,vacc.eff,testing.prob,test.sens,test.delay,contact.reduction,lambda,nSeeds){
  
  
  
  status.matrix <- data.frame(infected          = rep(0,n), # status variable: indicate if individuals are infectious (1), susceptible (0) or recovered (-1)
                              time.of.infection = NA,       # time of infection
                              infector          = NA,       # infector
                              severity          = -1,       # severity: 1 symptomatic infection, 2 asymptomatic infection
                              TimeSymptomOnset  = Inf,      # symptom onset date
                              IPLength          = NA,       # length of infectious period
                              Recovery          = Inf,      # recovery time 
                              TimePosTest       = Inf,      # time at which the test gives positive result
                              Vaccinated        = 0)        # Vaccination status: 1 vaccinated, 0 unvaccinated            
  
  
  
  events<-data.frame("NextCtc"=Inf,"TestPositive"=Inf, "Recovery"=Inf)
  
  infectives<-rep(0,n) # vector that indicates who is infectious at the current time: 1 infectious 0 non infectious
  current.time<-0
  index.contact<-rep(0,n) # vector that selects the individuals that have to propose a new social contact (global) - 1 yes 0 no
  time.events<-matrix(NA,1,3)
  
  #transmission parameter dataframe: each line is an individual, the first column is the transmission coefficient and the second is the contact rate
  transmission.parameters<-data.frame("id"=1:n,"q"=rep(q,n),"contact_rate"=rep(lambda,n))   
  
  #Proportion of immune/vaccination
  if (prop.immune>0){
    #status.matrix$infected[sample(1:n,round(prop.immune*n))]<--2 #-2 means they are immune # We drop this such that there's no immune individual
    status.matrix$Vaccinated[sample(1:n,round(prop.immune*n))]<-1
    }
    
  testpositive.day<-rep(Inf,n)
  contact.time<-data.frame("id"=1:n,"pr.ctc"=rep(NA,n),"pr.infectee"=rep(NA,n))   #matrix containing the proposed time of the next contact (first column) and the contact individual (second column)
  
  # first infected: randomly chosen in the population (among the susceptibilities)
  first.cases<-sample(which(status.matrix$infected==0),nSeeds) #those who are immune can get infected or not?
    
  # update the infectious status and epidemiological characteristics of the seeded cases  
  for (j in first.cases){
    first<-j
    status.matrix$infected[first] <- 1 
    status.matrix$time.of.infection[first] <- 0
  
  #if (status.matrix$Vaccinated[first]<-1){ #if a person is vaccinated  
    #status.matrix$IPLength[first]<-infectious.period.length()*vaccine.effectiveness #need to be distinguished between vaccinated and unvaccinated individuals (Leo add)
  #}else{
    status.matrix$IPLength[first]<-infectious.period.length()
  #}
    status.matrix$Recovery[first]<-current.time+status.matrix$IPLength[first]
    if (runif(1)<rho){ #if symptomatic
      status.matrix$severity[first]<-1
      if(status.matrix$Vaccinated[first]==1){ #if vaccinated
      status.matrix$TimeSymptomOnset[first]<-current.time+incubation.period.vacc(status.matrix$IPLength[first])
      }else{
        status.matrix$TimeSymptomOnset[first]<-current.time+incubation.period.unvacc(status.matrix$IPLength[first])
      }
      if (runif(1)<testing.prob){ #a patient is tested (Leo add)
        if (runif(1)<test.sens){ #the test is positive (Leo add)
          status.matrix$TimePosTest[first]<-status.matrix$TimeSymptomOnset[first]+test.delay #assume only symptomatic individuals are tested (Leo add)
        }
      }
      time.events<-rbind(time.events,c(current.time,1.1,first)) #1.1 represents a symptomatic infection event
    }else{
      status.matrix$severity[first]<-2
      transmission.parameters$q[first]<-transmission.parameters$q[first]*alpha.as 
      time.events<-rbind(time.events,c(current.time,1.2,first)) #1.2 represents an asymptomatic infection event
    }
    infectives[first]<-1
    contact.time$pr.ctc[first]<-rexp(1,transmission.parameters$contact_rate[first])+current.time       # Generate the next inter-arrival time for individual i 
  }
  
  proposed.individual<-0
  temp.contact.time<-0
  indiv.prop.ctc<-0
  recovered<-0
  err<-0
  
    while((sum(infectives))>0){ #while there are still infectives
    #Phase 1: individuals that has to, propose a new social contact
    for (i in which(index.contact==1) ){ # for all the individuals that has to propose a global contact    #why 1 because previously we assign index.contact to 0 (Leo ask)
      contact.time$pr.ctc[i]<-rexp(1,transmission.parameters$contact_rate[i])+current.time# I generate the next interarrival time for individual i
      index.contact[i]<-0 #i only propose 1 contact?
    }
    
    next.contact<-contact.time$pr.ctc #overall contact times    
    recovery.vector<-status.matrix$Recovery
    events$Recovery<-min(recovery.vector)
    testpositive.day<-status.matrix$TimePosTest
    #Identify the next event (the one that is realized first in time)
    ifelse(length(which(is.na(next.contact)==FALSE))>0,events$NextCtc<-min(next.contact, na.rm = T),events$NextCtc<-Inf) # among all the proposed social contact between households we select the minimum
    ifelse(length(which(!is.infinite(testpositive.day)))>0,events$TestPositive<-min(testpositive.day),events$TestPositive<-Inf ) # among all the test positive day we select the minimum
    
    next.evts<-colnames(events)[which(min(events)==events)] # if two at the same time we pick one random 
    if (length(next.evts)>1){ 
      next.evts<-sample(colnames(events)[which(min(events)==events)],1)
    }
    
    if (next.evts=="NextCtc"){
      current.time<-events$NextCtc
      if (length(min(contact.time, na.rm = T))>1){ #when two contacts happen at the same time we select one at random 
        infector<-sample(which(contact.time==current.time),1) 
        infectee<-sample(setdiff(1:n,infector),1)
        index.contact[infector]<-1 #infectors are thet one who have to propose a new contact
        contact.time$pr.ctc[infector]<-NA
      }else{
        infector<-which(contact.time$pr.ctc ==events$NextCtc)
        infectee<-sample(setdiff(1:n,infector),1)
        index.contact[infector]<-1
        contact.time$pr.ctc[infector]<-NA
      }
      inf.ctc<-transmission.parameters$q[infector]
      if (status.matrix$Vaccinated[infector]==1){
      acc.rate<-nCov.InfMeasure.vacc(t=current.time-status.matrix$time.of.infection[infector])*inf.ctc*(1-vacc.eff)
      }else{
        acc.rate<-nCov.InfMeasure.unvacc(t=current.time-status.matrix$time.of.infection[infector])*inf.ctc
      }
      if (status.matrix$infected[infector]!=1){acc.rate<-0}
      if (acc.rate>1){err<-err+1}
      if (status.matrix$infected[infectee]==0 & runif(1)<acc.rate){
        status.matrix$infected[infectee] <- 1 
        status.matrix$time.of.infection[infectee] <-current.time
        status.matrix$infector[infectee]<-infector
        status.matrix$IPLength[infectee]<-infectious.period.length()
        status.matrix$Recovery[infectee]<-current.time+status.matrix$IPLength[infectee]
        if (runif(1)<rho){ #if symptomatic
          status.matrix$severity[infectee]<-1
          if(status.matrix$Vaccinated[infectee]==1){ #if vaccinated
           status.matrix$TimeSymptomOnset[infectee]<-current.time+incubation.period.vacc(status.matrix$IPLength[infectee])
          }else{
            status.matrix$TimeSymptomOnset[infectee]<-current.time+incubation.period.unvacc(status.matrix$IPLength[infectee])
          }
          if (runif(1)<testing.prob){
            if (runif(1)<test.sens){
               status.matrix$TimePosTest[infectee]<-status.matrix$TimeSymptomOnset[infectee]+test.delay
            }
          }
          time.events<-rbind(time.events,c(current.time,1.1,infectee))
        }else{
          status.matrix$severity[infectee]<-2
          transmission.parameters$q[infectee]<-transmission.parameters$q[infectee]*alpha.as #A single q parameter for everyone
          time.events<-rbind(time.events,c(current.time,1.2,infectee))
        }
        infectives[infectee]<-1
        contact.time$pr.ctc[infectee]<-rexp(1,transmission.parameters$contact_rate[infectee])+current.time      # I generate the next interarrival time for individual i
      }
    }
    
    
    if (next.evts=="TestPositive"){
      current.time<-events$TestPositive
      testpositive.individuals<-which(testpositive.day==current.time)
      status.matrix$TimePosTest[testpositive.individuals]<-Inf
      for (k in testpositive.individuals){
        transmission.parameters$contact_rate[k]<-transmission.parameters$contact_rate[k]*contact.reduction
      }
    }
    
    
    if (next.evts=="Recovery"){
      current.time<-events$Recovery
      temp.recovered<-which(recovery.vector==events$Recovery)
      for (recovered in temp.recovered){
        status.matrix$infected[recovered]<--1
        recovery.vector[recovered]<-Inf
        status.matrix$Recovery[recovered]<-Inf
        infectives[recovered]<-0
        time.events<-rbind(time.events,c(current.time,-1,recovered))
      }
    }
  }
  
  time.events<-time.events[-1,]
  timev.name<-c("time","event","who")
  dimnames(time.events)<-list(NULL,timev.name)
  

  # Compute Epidemic Meausure
  Rt<-NULL
  first.cases<-which(status.matrix$time.of.infection==0)
  temp.sec.cases<-NULL
  for (j in first.cases){
    ifelse(length(which(status.matrix$infector==j)>0),temp.sec.cases<-c(temp.sec.cases,length(which(status.matrix$infector==j))),temp.sec.cases<-c(temp.sec.cases,0))
  }
  Rt<-mean(temp.sec.cases)
  
  NewCases<-nSeeds
  Prevalence<-nSeeds
  last.day<-ceiling(max(time.events[,1]))
  
  for (i in 1:last.day){
    temp.time<-setdiff(which(time.events[,1]>(i-1)),which(time.events[,1]>i))
    if (length(temp.time)>0){
      new.cases.ev<-c(which(time.events[temp.time,2]==1.1),which(time.events[temp.time,2]==1.2))
      if (length(new.cases.ev)>0){
        NewCases<-c(NewCases,(length(new.cases.ev)))
        newly.infected.id<-time.events[temp.time[new.cases.ev],3]
        temp.sec.cases<-NULL
        for (k in newly.infected.id) {
          ifelse(length(which(status.matrix$infector ==k)>0),temp.sec.cases<-c(temp.sec.cases,length(which(status.matrix$infector==k))),temp.sec.cases<-c(temp.sec.cases,0))
        }
        Rt<-c(Rt,mean(temp.sec.cases))
      }else{
        NewCases<-c(NewCases,0)
        Rt<-c(Rt,NA)
      }
      ctime<-which(time.events[,1]<i)
      Prevalence<-c(Prevalence, (length(which(time.events[ctime,2]==1.1))+length(which(time.events[ctime,2]==1.2))-length(which(time.events[ctime,2]==-1))))
    }else{
      ctime<-which(time.events[,1]<i)
      Prevalence<-c(Prevalence, (length(which(time.events[ctime,2]==1.1))+length(which(time.events[ctime,2]==1.2))-length(which(time.events[ctime,2]==-1))))
      NewCases<-c(NewCases,0)
      Rt<-c(Rt,NA)
    }
  }
  
  epi.evo<-data.frame("Days"=0:last.day, "Incidence"=NewCases,"Prevalence"=Prevalence,"Rt"=Rt)
  FinalSize<-data.frame("FinalSize"=length(which(status.matrix$infected==-1)))
  PeakIncidence<-data.frame("PeakIncidence"=max(epi.evo$Incidence),"TimePeakIncidence"=which(epi.evo$Incidence==max(epi.evo$Incidence))[1])
  PeakPrevalence<-data.frame("PeakPrevalence"=max(epi.evo$Prevalence),"TimePeakPrevalence"=which(epi.evo$Prevalence==max(epi.evo$Prevalence))[1])
  
  return(list(time.events=time.events, status.matrix=status.matrix,epi.evo=epi.evo, FinalSize=FinalSize, PeakIncidence=PeakIncidence, PeakPrevalence=PeakPrevalence))
  }














