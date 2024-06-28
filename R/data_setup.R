example_data <- function(){
  load("data/am2022.RData")
  load("data/sizeage_matrix.RData")
  input$sizeage <- sizeage; rm(sizeage)
  source("R/helper.R")

  head(input$obsdf, 5) # long format with all observations
  # obs_type # 0=catch, 1=index, 2=agecom, 3=lencomp
  # nll_type # 0=dnorm, 1=dmultinom
  # fit_data # 1/0=TRUE/FALSE
  # fleet    # 1=fishery, 2=survey
  # obs      # transformed appropriately for nll_type (becomes keep vec)
  # obserror # if nll_type obs error is an input (note this is Neff for dmultinom)

  #Remove length comps for now till we implement them
  input$obsdf <- input$obsdf[input$obsdf$obs_type!=3,]


  return()
}

example_dat<-function(input){
  # data list ----
  dat <- list()
  dat$obs <- input$obsdf$obs
  dat$aux <- input$obsdf
  dat$aux <- get_id(dat$aux)
  dat$aux <- get_likelihood_index(dat$aux)
  dat$year <- input$years
  dat$minYear <- min(dat$year)
  dat$age <-  input$ages
  dat$len <-  input$lens
  dat$minAge <- min(dat$age)
  dat$sampleTimes <- input$srv_frac
  dat$spawnTimes <- input$sp_frac
  dat$waa <- input$waa
  dat$mature <- input$maturity
  dat$sizeage <- input$sizeage
  dat$fleetTypes <- unique(input$obsdf$fleet)
  dat$srmode <- 0 #
  dat$logN_mode <- 0 # 0 = deterministic SCAA, 1 = sigR estimated, logR estimated as ranef, 2 = full state-space with/ shared sigN for a > 1 (same as n_NAA_sigma in WHAM)

  # prediction data frame
  dat$aux <- get_pred(dat$aux, input)
}

example_par<-function(input){
  # parameter ----
  par <- list()
  par$logsigR <- log(input$sigr)
  par$logsigN <- if(dat$logN_mode==2){log(0.5)}else{numeric(0)}
  par$logQ <- 0
  # is M a constant in FIMS or by year/age?
  par$logM <- matrix(log(input$natmort), nrow=length(dat$year), ncol=length(dat$age))
  par$rickerpar <- if(dat$srmode==1){c(1,1)}else{numeric(0)}
  par$bhpar <- if(dat$srmode==2){c(1,1)}else{numeric(0)}
  par$logN <- matrix(10, nrow=length(dat$year), ncol=length(dat$age)) # Tim suggested initializing at 10 rather than zero
  par$logFmort <- matrix(0, nrow=length(dat$year), ncol=1)
  par$logfshslx <- log(input$fsh_slx) # need parametric selectivity
  par$logsrvslx <- log(input$srv_slx)

  return(par)
}

example_map <- function(){
  fill_vals <- function(x,vals){rep(as.factor(vals), length(x))}
  map <- list()
  map$logsigR <- if(dat$logN_mode==0){fill_vals(par$logsigR, NA)}else{factor(1)}
  # map$logQ <- fill_vals(par$logQ, NA)
  map$logM <- fill_vals(par$logM, NA)
  map$logfshslx <- fill_vals(par$logfshslx, NA)
  map$logsrvslx <- fill_vals(par$logsrvslx, NA)

  nyr <- length(dat$year)
  nage <- length(dat$age)
  tmp <- matrix(data = NA, ncol = nage, nrow = nyr)
  tmp[,1] <- 1:nyr
  tmp[1,2:nage] <- (nyr+1):(nyr+nage-1)
  map$logN <- as.factor(as.vector(tmp))
  return(map)
}
