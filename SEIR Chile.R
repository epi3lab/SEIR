

# Population structure: (S, E, E_t, I, I_t, I_s, H, R, D)
#   S = susceptible
#   E = latent
#   E_t = latent but test-posI_tive
#   I = infectious
#   I_t = infectious and test-posI_tive
#   I_s = symptomatic
#   H = hospI_talized
#   R = recovered
#   D = dead
#####################
# Assumptions:
# 1. latent and infectious individuals can be identified by testing
#   1a. test-positive infectious individuals will reduce contacts 90% (modifiable)
# 2. latent individuals will become infectious, then symptomatic before seeking hospital care
#   2a. only a portion (age-specific) of symptomatic individuals will require hospitalization
# 3. recovery I_s based on age and disease status
#   3a. a portion of latent individuals will recover before being symptomatic
#   3b. symptomatic individuals who are not hospitalized will recover
#   3c. hospitalized individuals will either recover or die

params=list(
  # most from Harvard Document
  l=2.52,   #length of latency
  p_s=0.75, #0.96 probability of symptoms #https://www.nytimes.com/2020/03/31/health/coronavirus-asymptomatic-transmission.html
  i=2.7,   #length of infectiousness before symptomatic (subclinical)
  r=11,   #time to recovery if mild/asymptomatic
  p_h=0.14, #probability of hospitalization given symptoms, from https://pediatrics.aappublications.org/content/early/2020/03/16/peds.2020-0702
  s=5,   #length of symptoms before hospitalization
  h=10,  #length of hospitalization
  p_d=0.04 #probability of death given hospitalization
)
paramlist=list( #set up separate parameters for each age group
  params,params,params,params,params,params,params,params,params
)
paramlist[[1]]$p_s=0.1 #from https://pediatrics.aappublications.org/content/early/2020/03/16/peds.2020-0702
#other parameters can be made age-specific as needed

lambda_params=list(
  # from http://gabgoh.github.io/COVID/index.html
  p_i=.45,  #probability of infection given contact
  p_hcw=.25, #probability of infection of Health Care Workers given contact (could change after Personal Protective Equipment decrease)
  c_hcw=10, #number of contacts by Health Care Workers per hospitalized patient per day
  # modifiable 
  q=0.9, #proportional decrease in c due to quarantine after positive test
  d=0.9, #proportional decrease in c due to social distancing
  w=0.9 #proportional decrease in p_i due to hygiene
)

# from https://journals.plos.org/plosmedicine/article/file?id=10.1371/journal.pmed.0050074&type=printable
waifw=read.csv("contact_matrix.csv") 

cfr_byage=data.frame(
  #from https://www.cdc.gov/mmwr/volumes/69/wr/mm6912e2.htm
  #assume HCW similar to 22-54
  agegroup=colnames(waifw),
  p_h=c(2,.175,.175,.2475,.253,.36,.446,.508,.21),
  cfr=c(0,.0015,.0015,.0065,.02,.038,.074,.19,.004)
)
cfr_byage$p_d=cfr_byage$cfr/cfr_byage$p_h #convert p(d) to p(d|h)

#pops is a matrix with the current Chilean population in each age category, census data 2020 used by Minsal
#(S, E, E_t, I, I_t, I_s, H, R, D)

Lambda=function(lambda_params,waifw,pops){
  lambda=rep(NA,nrow(pops))
  beta=waifw*lambda_params$p_i
  I_mat=apply(cbind(pops[1:8,4]*(1-lambda_params$d), #untested infectious
                    pops[1:8,5:6]*(1-lambda_params$q)), #quarantined
              1,sum)/apply(pops[1:8,1:8],1,sum) #add up pressure from infectious groups, normalizes by group size
  for(i in 1:ncol(pops)){
    lambda[i]=sum(beta[1:8,i]*I_mat) #sum up infectious pressure from each age group
  }
  lambda[9]=lambda_params$p_hcw*sum(pops[,7])/sum(pops[9,])*lambda_params$c_hcw #average all hospitalized patients across HCWs
  return(lambda)
}


DZstep=function(pop,params,lambda){
  S=pop[1]; E=pop[2]; E_t=pop[3]; I=pop[4]; I_t=pop[5]; I_s=pop[6]; H=pop[7]; R=pop[8]; D=pop[9]
  # hospitalized
  newdeaths=min(H,rpois(1,params$p_d*(1/params$h)*H)) #how many will die
  D=D+newdeaths;H=H-newdeaths
  recoverH=min(H,rpois(1,(1-params$p_d)*(1/params$h)*H)) #hospitalized recover
  R=R+recoverH;H=H-recoverH
  # symptomatic
  hospitalize=min(I_s,rpois(1,params$p_h*(1/params$s)*I_s)) #symptomatic become hospitalized
  H=H+hospitalize;I_s=I_s-hospitalize
  recoverI_s=min(I_s,rpois(1,(1-params$p_h)*(1/params$r)*I_s)) #symptomatic recover at home
  R=R+recoverI_s;I_s=I_s-recoverI_s
  # infectious
  symptomsI=min(I,rpois(1,I*(1/params$i)*params$p_s)) #infectious become symptomatic
  symptomsI_t=min(I_t,rpois(1,I_t*(1/params$i)*params$p_s)) #tested infectious become symptomatic
  I_s=I_s+symptomsI+symptomsI_t;I=I-symptomsI;I_t=I_t-symptomsI_t
  recoverI=min(I,rpois(1,I*(1/params$i)*(1-params$p_s))) #infectious recover
  recoverI_t=min(I_t,rpois(1,I_t*(1/params$i)*(1-params$p_s))) #tested infectious recover
  R=R+recoverI+recoverI_t;I=I-recoverI;I_t=I_t-recoverI_t
  # latent
  infectious=min(E,rpois(1,E/params$l)) #latent become infectious
  I=I+infectious;E=E-infectious
  infectious_t=min(E_t,rpois(1,E_t/params$l)) #tested latent become infectious
  I_t=I_t+infectious_t;E_t=E_t-infectious_t
  # susceptible
  infection=min(S,rpois(1,S*lambda)) #susceptible become infected
  S=S-infection;E=E+infection
  #recover population
  newpop=c(S,E,E_t,I,I_t,I_s,H,R,D)
  return(list(pop=newpop,deaths=newdeaths,infections=infection))
}


#####Running a simulation
#source("SEIR_functions.R")

###############No control, Chile
N=19500000 #total population of Chile per proyection census 2017
agedist=c(0.261533392, 0.037854616,0.342484541,0.129616582,0.109695186,0.070296303,0.034908441,0.013610939, #proportion in each age category, per census #https://www.populationpyramid.net/es/chile/2019/
          109706) #number of HCW, per https://www.minsal.cl/wp-content/uploads/2015/08/Informe-Brechas-RHS-en-Sector-P%C3%BAblico_Abril2017.pdf 
agenums=c(round((N-agedist[9])*agedist[1:8],0),agedist[9]) #distribute numbers across age groups

#https://www.minsal.cl/wp-content/uploads/2015/08/Informe-Brechas-RHS-en-Sector-P%C3%BAblico_Abril2017.pdf
#Grafico 2
#41622 MD in Chile (2017)
#44473 Nurses (2017)
#12166 Tech med in (2017)
#10000 (bluff) medical assistants, by guessing

timestop=30*6 #run for 6 months

####one iteration
RunIter=function(N,agenums,timestop,waifw,params,lambda_params){
  poparray=array(NA,dim=c(9,9,timestop),dimnames = list(colnames(waifw),c("S","E","E_t","I","I_t","I_s","H","R","D"),paste("t",1:timestop,sep="")))
  #sets up an array for the population at each timestep in each age and disease category
  poparray[,,1]=0
  poparray[,1,1]=agenums #set up the starting population as fully susceptible
  startdz=5 #how many diseased to start with
  startdist=rmultinom(1,startdz,agedist[4:8]) #distribute the diseased across the older age categories
  poparray[4:8,1,1]=poparray[4:8,1,1]-startdist #take diseased out of S
  poparray[4:8,4,1]=startdist #put diseased in I
  #initialize saving
  deaths=infections=lambda=matrix(NA,nrow=timestop,ncol=9)
  colnames(deaths)=colnames(infections)=colnames(lambda)=colnames(waifw)
  seir=matrix(NA,nrow=timestop,ncol=9);colnames(seir)=c("S","E","E_t","I","I_t","I_s","H","R","D")
  #run the simulation
  for(t in 2:timestop){ #daily
    lambda[t,]=Lambda(lambda_params,waifw,poparray[,,t-1])
    #step each agegroup through infections
    for(g in 1:9){
      tstep=DZstep(poparray[g,,t-1],paramlist[[g]],lambda[t,g])
      deaths[t,g]=tstep$deaths
      infections[t,g]=tstep$infections
      poparray[g,,t]=tstep$pop
    }
    seir[t,]=apply(poparray[,,t],2,sum)
  }
  return(list(poparray=poparray,seir=seir,deaths=deaths,infections=infections))
}

#######multiple iterations
n_iter=1000 #number of iterations
iter_out=list(NULL) #initialize
seir_array=array(NA,dim=c(timestop,9,n_iter),
                 dimnames = list(paste("t",1:timestop,sep=""),c("S","E","E_t","I","I_t","I_s","H","R","D"),1:n_iter))
death_array=infection_array=array(NA,dim=c(timestop,9,n_iter),
                                  dimnames = list(paste("t",1:timestop,sep=""),colnames(waifw),1:n_iter))

for(i in 1:n_iter){
  iter_out[[i]]=RunIter(N,agenums,timestop,waifw,params,lambda_params)
  seir_array[,,i]=iter_out[[i]]$seir
  death_array[,,i]=iter_out[[i]]$deaths
  infection_array[,,i]=iter_out[[i]]$infections
}

matplot(infection_array[,7,],type = "l",col = "gray") #plot hospitalized
matplot(seir_array[,9,],type = "l",col = "gray") #plot deaths
matplot(seir_array[,6,],type = "l",col = "gray") #plot Infected
