#length/age model process test

##SECTION 1: Setup basic values for an example model

#Generic growth data to create length at age this could be replaced
#by empirical length at age in practice if available but provides
#an easy way to parameterize a generic example population
max_age <- 15
Linf <- 100
K <- .3
a0 <- 0
cv <- .2
alpha <- 0.001
beta <- 3
rzero <- 10000
steep <- 1
phi_0 <- 560.049
age_50_mat <- 6
mat_slope <- 2
M <- 0.2

#Set up length/age composition data ranges
comp_ages <- 0:max_age
comp_lengths <- 0:ceiling((1+3*cv)*Linf)
n_length_comps <- length(comp_lengths)

#Set up number of growth morphs to track
n_growth_morphs <- 7
length_props <- 1+(-3:3)*cv

#Vector to specify the ages tracked in the model
#For simplicity start with assumption that ages in
#comps and model are the same
ages <- comp_ages
n_ages <- length(ages)

#Vector of model years
start_year <- 1
end_year <- 10
years <- start_year:end_year
n_years <- length(years)

#Vector of annual recruitment values we need n_years + n_ages length vector
#so that initial depletion comps have a recruitment to track back to.
#pre start year recruitments and initial depletion values will be confounded
#so one of these needs to be fixed (i.e. pre start year recruitment set to R0
# or initial depletion set at an equilibrium F)
recruits <- rep(rzero,(n_years+n_ages))
recruit_year <- c(rev(years[1]-1-ages),years)
#Initialize ssb
ssb <- rep(0,n_years)
unfished_ssb <- rep(0,n_years)

##SECTION 2: Model sub functions used throughout script

#Function to get expected length at age from growth params
AtoL <- function(a_f,Linf_f,K_f,a_0){
  L <- Linf_f*(1-exp(-K_f*(a_f-a_0)))+0.001
}

#Function to get expected weight at length from growth params
LtoW <- function(L_f,alpha_f,beta_f){
  W <- alpha_f*L_f^beta_f
}

logistic <- function(x_f, inflection_point_f, slope_f){
  prop <- 1 / (1 + exp(-1 * slope_f * (x_f - inflection_point_f)))
}

get_recruit <- function(ssb_f,rzero_f,steep_f,phi_0_f){
  recruits <- (0.8*rzero_f*steep_f*ssb_f)/(0.2*phi_0_f*rzero_f*(1.0-steep_f)+ssb_f*(steep_f-0.2))
}

##SECTION 3a: Setup population model vectors to track population bin values these
# should all be length n_years*n_ages*n_growth_morphs

#For simplicity start with the assumption that these are not time varying
#Assume that an equal number of lengths are tracked
#for each age. This simplifies setting things up but I
#don't think will be stricktly necessary

year <- sort(rep(years,n_ages*n_growth_morphs))

#Assume that an equal number of lengths are tracked
#for each age. This simplifies setting things up but I
#don't think will be strictly necessary
ages_x_lengths <- sort(rep(ages,n_growth_morphs))

age <- rep(ages_x_lengths,n_years)

#For this example I'm setting the lengths tracked for each
#age based on a VB growth curve and 3 standard deviations
#either side of the mean in 1 sd steps. This should mean that
#almost all individuals (99.7%) are inside these bounds
#assuming a normal distribution.
lengths_x_ages <- AtoL(ages_x_lengths,
                       rep(Linf,length(ages_x_lengths)),
                       rep(K,length(ages_x_lengths)),
                       rep(a0,length(ages_x_lengths))
                       )*rep(length_props,length(ages))

len <- rep(lengths_x_ages,n_years)

#For this example set a single natural mortality for all years, ages, and lengths
nat_mort <- rep(M,n_years*n_ages*n_growth_morphs)

#Total mortality will be calculated as a sum of nat_mort and fleet specific
#fishing mortalities
total_mort <- rep(0,n_years*n_ages*n_growth_morphs)

#Vectors to hold depletion data specific to population year, age, and length bins
#Abundance will be interpolated across all comp lengths later using this
#depletion

depletion <- rep(0,n_years*n_ages*n_growth_morphs)
unfished_depletion <- rep(0,n_years*n_ages*n_growth_morphs)

#SECTION 3b: Setup fleet and survey vectors to be tracked by the population
#model

#Setup fleet and survey structure
n_fleets <- 2
n_surveys <- 1

if((n_fleets+n_surveys)>0){
  fleets <- 1:(n_fleets+n_surveys)
}else{
  fleets <- NULL
}

Fleet_F <- c(0.15,0.05,0)

fleet_fleet <- sort(rep(fleets,n_years*n_ages*n_growth_morphs))
fleet_year <- rep(year,(n_fleets+n_surveys))
fleet_age <- rep(age,(n_fleets+n_surveys))
fleet_length <- rep(len,(n_fleets+n_surveys))
fleet_type <- sort(rep(c(rep(1,n_fleets),rep(2,n_surveys)),n_years*n_ages*n_growth_morphs))

fleet_fish_mort <- rep(0,n_years*(n_fleets+n_surveys)*n_ages*n_growth_morphs)
for(i in seq_along(Fleet_F)){
  fleet_fish_mort[fleet_fleet==i] <- Fleet_F[i]
}

fleet_catch_proportion <- rep(0,n_years*(n_fleets+n_surveys)*n_ages*n_growth_morphs)

fleet_abund_index <- rep(0,n_years*(n_fleets+n_surveys)*n_ages*n_growth_morphs)

fleet_q <- rep(0,n_years*(n_fleets+n_surveys)*n_ages*n_growth_morphs)

##SECTION 3c: Setup initial population structure vectors to specify fished
#and unfished initial depletion

#Initialize the fished and unfished depletion at age and length in first year
#of the population model
init_age <- ages_x_lengths
init_length <- lengths_x_ages
init_F <- M
init_DAAL <- exp(-(M+init_F)*init_age) #rep(0,n_ages*n_growth_morphs)
init_unfished_DAAL <- exp(-(M)*init_age) #rep(0,n_ages*n_growth_morphs)


##SECTION 3d: Setup population composition vectors that will track length and
#age comp calculations of actual abundance used to fit observed data and
#calculate ssb and expected recruitment.

#First set up vectors to store indexing values of year, age, composition length,
#the model length for lower interpolation bound, model length for upper interpolation
#bound, and the abundance for the composition component.
pop_length_comp_year <- sort(rep(years,(n_ages*n_length_comps)))
pop_length_comp_age <- rep(sort(rep(ages,n_length_comps)),n_years)
pop_length_comp_length <- (rep(comp_lengths,n_ages*n_years))
#For this example I'm setting the weights as fixed weight at length.
pop_length_comp_weight <- LtoW(pop_length_comp_length,
                               rep(alpha,length(pop_length_comp_length)),
                               rep(beta,length(pop_length_comp_length)))
#For this example set constant 50% female T
pop_length_comp_prop_female <- rep(0.5,n_years*n_ages*n_length_comps)
#For this example set proportion mature to logistic maturity by age
pop_length_comp_prop_mature<-logistic(pop_length_comp_age,age_50_mat,mat_slope)

#Setup fleet composition matrix using selectivity by age,length
fleet_comp <- matrix(0,nrow=length(pop_length_comp_year),ncol=(n_fleets+n_surveys))
fleet_select_prop <- matrix(0,nrow=length(pop_length_comp_year),ncol=(n_fleets+n_surveys))
fleet_age_50_select <- c(2,3,4)
fleet_age_slope_select <- c(2,2,2)
fleet_length_50_select <- c(40,60,80)
fleet_length_slope_select <- c(10,10,10)
for(i in seq_along(fleet_select_prop[1,])){
  fleet_select_prop[,i] <- logistic(pop_length_comp_age,
                                    fleet_age_50_select[i],
                                    fleet_age_slope_select[i]) *
                           logistic(pop_length_comp_length,
                                    fleet_length_50_select[i],
                                    fleet_length_slope_select[i])
}
#Setup temporary storage vectors to calculate interpolation bounds and abundance
temp_abun_per_rec<-NULL
temp_comp_lower_interp<-rep(NA,n_length_comps)
temp_comp_upper_interp<-rep(NA,n_length_comps)
temp_comp_weight_upper_interp<-rep(NA,n_length_comps)
comp_lower_interp <- NULL
comp_upper_interp <- NULL
comp_weight_upper_interp <- NULL

pop_length_comp_lower_interp <- rep(NA,length(pop_length_comp_length))
pop_length_comp_upper_interp <- rep(NA,length(pop_length_comp_length))
pop_length_comp_abund <- rep(NA,length(pop_length_comp_length))
pop_length_comp_abund_unfished <- rep(NA,length(pop_length_comp_length))

##SECTION 3e: Setup model abundance transition function to identify the source
#and sink for calculating depletion at year/age/length for this growth morph
#setup these just end up being 100% transfers along a cohort line.
growth_source <- NULL
growth_sink <- NULL
for(i in seq_along(years)){
  for(j in rev(seq_along(ages))){
    for(k in 1:n_growth_morphs){
      if(i==1){
        #In first year no transfer so just leave in same year
        #an initial population function with fill these.
        growth_source<-c(growth_source,which(year==years[i] &
                                             age==ages[j])[k])

        growth_sink<-c(growth_sink,which(year==years[i] &
                                         age==ages[j])[k])
      }else if(j==1){
        #The first age class will be filled using a recruitment function.
        growth_source<-c(growth_source,which(year==years[i] &
                                               age==ages[j])[k])

        growth_sink<-c(growth_sink,which(year==years[i] &
                                           age==ages[j])[k])
      }else{
        #For all other years this will define the proportional transfer
        #between age/length bins each year. Currently simplified to 100%
        #along a range of growth morphs. Abundance will be adjusted by
        #natural and fishing mortality impacts.
        growth_source<-c(growth_source,which(year==years[i-1] &
                                               age==ages[j-1])[k])

        growth_sink<-c(growth_sink,which(year==years[i] &
                                           age==ages[j])[k])
      }
    }
  }
}
#Specify the transfer of abundance between year/age/length bins
transition_proportion <- rep(1,length(growth_source))

##SECTION 3f: Setup composition interpolation tracking values including proportion
#of recruits by length and link model bins to comp values for interpolation.

#Sequence along all the model ages to calculate the mean length at that age
#and a normal distribution of lengths around the mean as well as the closest
#model lengths at age to identify interpolation bounds.
for(i in seq_along(ages)){

  #Calculate mean length at age to spread lengths around
  mean_length <- AtoL(ages[i],Linf,K,a0)

  #Calculate the cumulative proportion shorter than each composition length
  temp_len_probs<-pnorm(q=comp_lengths,mean=mean_length,sd=mean_length*cv)
  #Reset the first length proportion to zero so the first bin includes all
  #density smaller than that bin
  temp_len_probs[1]<-0
  #subtract the offset length propabilities to calculate the proportion in each
  #bin. For each length bin the proportion is how many fish are larger than this
  #length but shorter than the next bin length.
  temp_len_probs<- c(temp_len_probs[-1],1)-temp_len_probs

  #Store the probability of length by bin as a zero mortality abundance
  #distribution this can be multiplied by the interpolated depletion estimates
  #from the population grid. This way the smooth normal shape is retained and
  #just the relative depletion is interpolated which should be much smoother.
  temp_abun_per_rec <- c(temp_abun_per_rec,temp_len_probs)

  #Calculate the population bin lengths. This could also be pulled from the
  #model data but this works for this example.
  pop_lengths <- mean_length*length_props


  #For the first bin set the lower bound at 1 for comp bins that fall outside
  #of the smallest population bin we will add a length zero bin with zero abundance
  temp_comp_lower_interp[comp_lengths<pop_lengths[1]] <- 0
  temp_comp_upper_interp[comp_lengths<pop_lengths[1]] <- pop_lengths[1]
  #Now loop over population length bins to determine which comp bins they should
  #be applied to for interpolation.
  for(j in seq_along(pop_lengths)[-1]){
    temp_comp_lower_interp[comp_lengths<pop_lengths[j] &
                             comp_lengths>=pop_lengths[j-1] ] <- pop_lengths[j-1]

    temp_comp_upper_interp[comp_lengths<pop_lengths[j] &
                             comp_lengths>=pop_lengths[j-1] ] <- pop_lengths[j]
  }
  #For length bins larger than the largest population bin set the lower bin to
  #the largest population bin and another bin will be added at twice the maximum
  #of the largest population/comp length and set at zero abundance.
  temp_comp_lower_interp[comp_lengths>=pop_lengths[length(pop_lengths)]] <- pop_lengths[length(pop_lengths)]
  temp_comp_upper_interp[comp_lengths>=pop_lengths[length(pop_lengths)]] <- 2*max(pop_lengths,comp_lengths)

  #For linear interpolation calculate the proportion weighting of the upper bound
  #depletion. depletion=upper_depletion*upper_weight+lower_depletion*(1-upper_weight)
  temp_comp_weight_upper_interp <- (comp_lengths-temp_comp_lower_interp)/
    (temp_comp_upper_interp-temp_comp_lower_interp)
  #compound the interpolation bounds across all ages
  comp_lower_interp <- c(comp_lower_interp,temp_comp_lower_interp)
  comp_upper_interp <- c(comp_upper_interp,temp_comp_upper_interp)
  comp_weight_upper_interp <- c(comp_weight_upper_interp,temp_comp_weight_upper_interp)
}
#replicate abundance and interpolation values across years
pop_length_comp_no_mort_abund_per_rec <- rep(temp_abun_per_rec,length(years))
pop_length_comp_lower_interp <- rep(comp_lower_interp,length(years))
pop_length_comp_upper_interp <- rep(comp_upper_interp,length(years))
pop_length_comp_weight_upper_interp <- rep(comp_weight_upper_interp,length(years))



##SECTION 4a: Calculate population depletion by year/age/length

#Loop over all population bins to calculate depletion, relative
#landings, and relative abundance indicies.
for(i in seq_along(transition_proportion))
{
  #Sink is a row reference to the depletion bin being calculated year, age, length
  sink <- growth_sink[i]
  #Source is a row reference to the bin used for starting depletion (i.e. For
  #this example case last year, one year younger, same length cohort)
  source <- growth_source[i]
  #The the proportion of the source bin to start with (i.e. 1 for this example)
  #In a scenario with multiple source bins this would be a fraction.
  trans_prop <- transition_proportion[i]
  #Same as above but for referencing the removals and abundance index of each
  #fleet and/or survey.
  fleet_sink <- which(fleet_year==year[sink] &
                      fleet_age==age[sink] &
                      fleet_length==len[sink])
  fleet_source <- which(fleet_year==year[source] &
                        fleet_age==age[source] &
                        fleet_length==len[source])

  if(year[sink]==start_year){ #This sets up initial conditions in the first model year
    #initial depletion at age/length is an input or estimated parameter.
    depletion[sink] <- init_DAAL[which(init_age==age[sink] & init_length==len[sink])]
    unfished_depletion[sink] <- init_unfished_DAAL[which(init_age==age[sink] & init_length==len[sink])]
  }else if(age[sink]==0){ #This sets initial depletion at age 0 recruitment to the
    #annual recruitment as a fraction of rzero
    depletion[sink] <- recruits[recruit_year==year[sink]]/rzero
    unfished_depletion[sink] <- recruits[recruit_year==year[sink]]/rzero
  }else if(age[sink]==max_age){ #This calculates plus group depletion by applying
    #mortality to the existing plus group depletion then adding incoming abundance
    #from the previous age/year

    #Total mortality is just the sum of natural and fishing mortality specified
    #for this year/age/length/fleet
    total_mort[sink] <- sum(nat_mort[sink],fleet_fish_mort[fleet_sink])

    #Calculate what fraction of abundance is removed by each fleet using fleet
    #specific F values for the existing plus group abundance
    fleet_catch_proportion[fleet_sink] <- fleet_catch_proportion[fleet_sink] + (fleet_fish_mort[fleet_sink]/total_mort[sink])*depletion[sink]*(1-exp(-total_mort[sink]))

    #Calculate population abundance index estimated from fleet specific q
    # for the existing plus group abundance
    #this allows surveys to still get an index estimate with F=0
    fleet_abund_index[fleet_sink] <- fleet_abund_index[fleet_sink] + fleet_q[fleet_sink]*depletion[sink]*(1-exp(-total_mort[sink]))

    #This applies mortality to the existing plus group abundance
    depletion[sink] <- depletion[sink]*exp(-total_mort[sink])

    #This applies natural mortality to the existing unfished plus group abundance
    unfished_depletion[sink] <- unfished_depletion[sink]*exp(-nat_mort[sink])

    #This calculates total mortality on the age group moving into the plus group
    total_mort[source] <- sum(nat_mort[source],fleet_fish_mort[fleet_source])

    #This applies mortality to the age group moving into the plus group and adds
    #that abundance to the plus group
    depletion[sink] <- depletion[sink] + trans_prop*depletion[source]*exp(-total_mort[source])

    #This applies mortality to the age group moving into the unfished plus group
    #and adds that abundance to the plus group
    unfished_depletion[sink] <- unfished_depletion[sink] + trans_prop*unfished_depletion[source]*exp(-nat_mort[source])

    #This calculates fleet specific removals for incoming age group and adds to
    #the existing plus group landings
    fleet_catch_proportion[fleet_sink] <- fleet_catch_proportion[fleet_sink] + (fleet_fish_mort[fleet_source]/total_mort[fleet_source])*trans_prop*depletion[source]*(1-exp(-total_mort[source]))

    #This calculates fleet specific abundance index for incoming age group and adds to
    #the existing plus group index
    fleet_abund_index[fleet_sink] <- fleet_abund_index[fleet_sink] + fleet_q[fleet_source]*trans_prop*depletion[source]*(1-exp(-total_mort[source]))

  }else{
    total_mort[source] <- sum(nat_mort[source],fleet_fish_mort[fleet_source])

    depletion[sink] <- trans_prop*depletion[source]*exp(-total_mort[source])

    unfished_depletion[sink] <- trans_prop*unfished_depletion[source]*exp(-nat_mort[source])

    fleet_catch_proportion[fleet_sink] <- fleet_catch_proportion[fleet_sink] + (fleet_fish_mort[fleet_sink]/total_mort[source])*trans_prop*depletion[source]*(1-exp(-total_mort[source]))

    fleet_abund_index[fleet_sink] <- fleet_abund_index[fleet_sink] + fleet_q[fleet_source]*trans_prop*depletion[source]*(1-exp(-total_mort[source]))
  }
}

##SECTION 4b: Calculate population composition abundance, ssb, and expected recruitment

#Here we calculate interpolated relative abundance at length
#param len #is the length that relative abundance is being calculated for
#param year #is the time at which this is being calculated (fixed at whole years for now)
#param pop_val_id #is a vector of values identifying the row id of model estimates
#                  to include in the interpolation

#Loop over all length components to calculate abundance
for(i in seq_along(pop_length_comp_year)){

  #Absolute abundance in compostion bins is calculated as
  #rzero*proportion_at_length_by_age*interpolated_depletion
  pop_length_comp_abund[i] <- rzero *
                              pop_length_comp_no_mort_abund_per_rec[i] * (

                              (1-pop_length_comp_weight_upper_interp[i]) *
                                max(0,depletion[year==pop_length_comp_year[i] &
                                          age==pop_length_comp_age[i] &
                                          len==pop_length_comp_lower_interp[i]])
                              +
                              pop_length_comp_weight_upper_interp[i] *
                                max(0,depletion[year==pop_length_comp_year[i] &
                                          age==pop_length_comp_age[i] &
                                          len==pop_length_comp_upper_interp[i]])
                          )

  pop_length_comp_abund_unfished[i] <- rzero *
                                    pop_length_comp_no_mort_abund_per_rec[i] * (

                                    (1-pop_length_comp_weight_upper_interp[i]) *
                                      max(0,unfished_depletion[year==pop_length_comp_year[i] &
                                                  age==pop_length_comp_age[i] &
                                                  len==pop_length_comp_lower_interp[i]])
                                    +
                                      pop_length_comp_weight_upper_interp[i] *
                                      max(0,unfished_depletion[year==pop_length_comp_year[i] &
                                                  age==pop_length_comp_age[i] &
                                                  len==pop_length_comp_upper_interp[i]])
                                  )

  #Identify which year this comp bin is part of
  temp_year_index <- which(years==pop_length_comp_year[i])

  #Add this comp bin to spawning stock biomass based on abundance, weight, sex, and maturity
  ssb[temp_year_index] <- ssb[temp_year_index] + pop_length_comp_abund[i]*
                                                 pop_length_comp_weight[i]*
                                                 pop_length_comp_prop_mature[i]*
                                                 pop_length_comp_prop_female[i]

  unfished_ssb[temp_year_index] <- unfished_ssb[temp_year_index] +
                                   pop_length_comp_abund_unfished[i]*
                                   pop_length_comp_weight[i]*
                                   pop_length_comp_prop_mature[i]*
                                   pop_length_comp_prop_female[i]
}

#Calculate expected compositions for fleets based on fleet selectivity at length
#and age. This is not currently integrated with the bin specific F's which are
#independent parameters. Not sure of the best approach to blend a state space
#based estimation of F by year/age/length/fleet vs the parametric selectivity
#used to specify selectivity proportions. Maybe some sort of penalty
#such as comparing fleet_select_prop to F/max(F) or something like that?
for(i in seq_along(fleet_comp[1,])){
    fleet_comp[,i] <- pop_length_comp_abund*
                      fleet_select_prop[,i]
}

#Predicted recruitment from a BH recruit function using SSB.
pred_unfished_recruitment <- get_recruit(unfished_ssb,rzero,steep,phi_0)
pred_recruitment <- get_recruit(ssb,rzero,steep,phi_0)

##SECTION 5: Likelihood calculations.
#TODO: Setup likelihood components including predicted recruitment then loop
#over priors and observations to sum likelihood components

NLL <- 0
parameters <- vector()
priors <- vector()

for(i in seq_along(parameters)){

}

for(i in seq_along(data_sources)){

}



