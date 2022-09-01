rm(list = ls())

library(ggplot2)
library(grid)
library(gridExtra)

workdir=getwd()

# number of R(t) values (1-4)
nb = 4

# Load the fortran codes

dyn.load("gsrc/detsirh.so")
is.loaded("detsirh")

dyn.load("gsrc/stosirh.so")
is.loaded("stosirh")

dyn.load("gsrc/mcmc.so")
is.loaded("mcmc")

# change this as needed to the path with your data 

data_path = "data/DICE_state_hosp.rds"

# load the data
state_hosp = readRDS(data_path)
# examine the structure
#str(state_hosp, 1)

# get a list of state abbvs
all_states = as.character(state_hosp$states$state)

# I kept all 4 influenza-specific columns, but the one we want for the 
# FluSight challenge is 'conf_inf_admin' confirmed influenza admissions
#str(state_hosp$conf_inf_admin)

# to be consistent with the DICE format, this is a data.frame of incidence
# with each column corresponding to a location. The dates for this incidence
# are a separate vector state_hosp$dates

# If you want to remove a location (American Samoa)
remove_loc = "AS"
remove_index = state_hosp$states$state==remove_loc
state_hosp$states = state_hosp$states[!remove_index, ]
state_hosp$conf_inf_admin = state_hosp$conf_inf_admin[, !remove_index]
state_hosp$conf_inf_admin_cvrg = state_hosp$conf_inf_admin_cvrg[, !remove_index]
state_hosp$inf_death = state_hosp$inf_death[, !remove_index]
state_hosp$inf_death_cvrg = state_hosp$inf_death_cvrg[, !remove_index]

all_states = as.character(state_hosp$states$state)

# to loop through individual locations:
dates = state_hosp$dates
for (ii in 1:length(state_hosp$states$state)) {
  # retrieve state name, population, and daily hospitalization incidence
  state_abbv = as.character(state_hosp$states$state[ii])
  state_pop = state_hosp$states$pop[ii]
  state_inc = state_hosp$conf_inf_admin[, ii]
}

# get a list of last NA
min_date = rep(Sys.Date(), length(all_states))
for (ii in 1:length(all_states)) {
  state_inc = state_hosp$conf_inf_admin[, ii]
  min_date_temp = dates[max(which(is.na(state_inc)))]
  min_date[ii] = min_date_temp
}


# removing everything before 2020-11-01 will remove spotty data
dates = state_hosp$dates
n_dates = length(dates)
#min_date_ind = which(dates==as.Date(max(min_date)))
min_date_ind = which(dates==as.Date("2021-11-01"))
keep_ind = min_date_ind:n_dates
state_hosp$dates = state_hosp$dates[keep_ind]
state_hosp$conf_inf_admin = state_hosp$conf_inf_admin[keep_ind, ]
state_hosp$conf_inf_admin_cvrg = state_hosp$conf_inf_admin_cvrg[keep_ind, ]
state_hosp$inf_death = state_hosp$inf_death[keep_ind, ]
state_hosp$inf_death_cvrg = state_hosp$inf_death_cvrg[keep_ind, ]

# -----------------------------------

nstate = length(state_hosp$states$state)

# calculate baseline for each state

basetime = 7

baseline = rep(NA, nstate)
for (ii in 1:nstate) {
  state_inc = state_hosp$conf_inf_admin[, ii]
  baseline[ii] = mean(state_inc[1:basetime], na.rm = TRUE)
  baseline[ii] = max(1/basetime, baseline[ii])
}

# -----------------------------------

# dates array 
dates = state_hosp$dates

# create a list for the results

all_state_hosp_frcst =list()

today = Sys.Date()

results_dir = "~/Dropbox/CSMB01/CDC2021-2022/weekly_profiles/"

results_dir = ''

# adjust filenames based on one or two value R(t)

plot_filename = paste0(results_dir, 'gstate',nb,'_forecast_', today,'.pdf')
data_filename = paste0(results_dir, 'gstate',nb,'_forecast_', today,'.rds')
cat("\n\n Fitting influenza daily admission using ", nb, " values for R(t)\n\n")

pdf(file = plot_filename)
par(mfrow = c(3, 1))

for (istate in 1:nstate) { 

  # retrieve state name, population, and daily hospitalization incidence
  state_abbv = as.character(state_hosp$states$state[istate])
  state_pop = state_hosp$states$pop[istate]
  state_inc = state_hosp$conf_inf_admin[, istate]
  
  population = Ntotal = state_pop
 
  obs = state_inc
  
  # in case there is an NA replace with zero
  obs[is.na(obs)] = 0
  
  ndays = length(dates)
  ndates = length(dates)

  # set number of forecast days to four weeks plus two days  (for safety)
  ndays_forecast = 7*4 +2
  ndays_total = ndays + ndays_forecast 
  dates_total = seq(from=dates[1], by ='day', length=ndays_total)


  # Check if to skip the fitting
  if (sum(state_inc, na.rm = TRUE) < 10) {
    data = data.frame(dates = dates, obs = obs)
    state_hosp_frcst = list(state = state_abbv, data = data, profiles = list(dates=dates_total, profiles = NULL), fit = FALSE)    
    all_state_hosp_frcst[[istate]] = state_hosp_frcst
    next 
  }
  
  cat("Processing: ", state_abbv, '\n')
  I0total = 1000
  
  dt = 0.05
    
  xvals = seq(0, ndays, 1)
  vecTcalc=seq(0, ndays,dt)
  nTimes = length(vecTcalc)
  
  vecN=Ntotal
  vecI0=I0total
  
  
  # pC - percent clinical - will be fitted 

  pC = 0.5
  
  # probability of going from I2_clinical to Hospital - do not fit this if fitting pC
  
  vecPH=0.03 
  
  # inverse rate of recovery for mild and clinical cases ( 3 days)
  D_IM = D_IC = 3.0
  
  # inverse rate of going from clinical to hospital 
  D_IH = 3.0
  
  # inverse rate of going from hospital to either recovered or death - needed only for book keeping 
  
  D_HOS=10.
  
  
  ## reporting rate - not fitted 
    
  rho = 1.0 
  
  ## fraction susceptible - set to 1 if not fitting
  ##
  eta = 1.0

  # transmission rates - will be fitted 
  
  R0=1.1

  ##
  # array for beta values - will be fitted
  beta_vec = rep(R0/D_IM, nb)

  ## time of change in beta - will be fitted 
  ##
  tcng_vec = rep(round(ndays/(2*nb)), nb)
    
  ## pack the parameter in the par array
  ##
  
  parnames = c("pC",'pH', "D_IC", "D_IM", "D_IH", "D_HOS", "rho", "eta", 'baseline')

  nparam = length(parnames)


  betanames = paste0('beta', 1:nb)
  tcngnames = paste0('tcng', 1:nb)


  parnames = c(parnames, betanames, tcngnames)
  
   
  par =rep(NA, (nparam + nb * 2))
  

  par[1] = pC
  par[2] = vecPH
  par[3] = D_IC
  par[4] = D_IM
  par[5] = D_IH
  par[6] = D_HOS
  par[7] = rho
  par[8] = eta
  par[9] = round(baseline[istate] /2)


  for (ii in (nparam+1):(nparam+nb)) {
    par[ii] = beta_vec[(ii-nparam)]
  }

  for (ii in (nparam+nb+1):(nparam+nb*2)) {
    par[ii] = tcng_vec[(ii-nparam-nb)]
  }
  

  names(par) <- parnames

  nTimes = length(vecTcalc)
  new.seed <- as.integer(runif(1)*2e4)
    
  parmin = parmax = par
  
  #pC
  parmin[1] = 0.05
  parmax[1] = 1.0
  
  #pH
  parmin[2] = 0.01
  parmax[2] = 0.1
  
  #D_IC - rate of recovery in Ic
  parmin[3] = 3.
  parmax[3] = 3.
  
  #D_IM - rate of recovery in Im
  parmin[4] = 3.
  parmax[4] = 3.
  
  #D_IH - rate from Ic to hospital
  parmin[5] = 3.
  parmax[5] = 3.
  
  #D_HOS - rate of exiting hospital- does not affect the fit
  parmin[6] = D_HOS
  parmax[6] = D_HOS
  
  # rho - reporting rate
  
  parmin[7] = 0.1
  parmax[7] = 1.
  
  # eta - fraction susceptible
  
  parmin[8] = 0.1
  parmax[8] = 1.0
  
  # baseline - will be fitted
  parmin[9] = 0
  parmax[9] = baseline[istate]

  # beta_vec
  parmin[(nparam+1):(nparam+nb)] = 0.9/3.
  parmax[(nparam+1):(nparam+nb)] = 1.30/3.
  
  # time of change
  
  parmin[(nparam+nb+1):(nparam+nb*2)] = 1
  parmax[(nparam+nb+1):(nparam+nb*2)] = ndays-5

  # 
  nparamtot = nparam + nb * 2

  # step size of parameters 
  
  step = rep(1e-2, nparamtot)
  
  
  # logflag: 0=uniform 1=log integer
  
  ilogflag = rep(0, nparamtot)
  
  # imask - 0 fixed parameter, 1 fit
  
  imask = rep(1, nparamtot)
  imask[1:nparam] = 0
  
  # fit only pC and baseline (along with beta)
  imask[c(1, 9)] = 1

  # a weight array for all the days - start with weigths of 1 for all days
  
  wght = rep(1., length(obs))
  
  arr_daily = rep(0, ndays)
  
  gamaObs = lgamma((round(obs) + 1))
  
  nMCMC = 1e6
  
  nlines = 1e4
  
  nlines = min(nMCMC, nlines)
  
  ithin = nMCMC/ nlines
  
  tab = array(0, c(nMCMC/ithin, (nparamtot+1)))
  
  # add noise to all parameters that are being optimized
  savepar = par
  
  for (ii in 1:nparamtot) {
    if (imask[ii] == 0) next
    par[ii] = parmin[ii] + (parmax[ii] - parmin[ii]) * runif(1, 0, 1)
    if (par[ii] >= parmax[ii] || par[ii] <= parmin [ii]) par[ii] = savepar[ii] 
  }


  cat("\n Initial Parameters: ", par, '\n\n')
    out <- .Fortran("mcmc", vecN = as.double(vecN), vecI0 = as.double(vecI0), nparam = as.integer(nparam), nb = as.integer(nb), par = as.double(par), parmin = as.double(parmin),
                    parmax = as.double(parmax), step = as.double(step), ilogflag = as.integer(ilogflag), imask = as.integer(imask), iseed = as.integer(new.seed), wght = as.double(wght), obs = as.double(obs),
                    gamaObs = as.double(gamaObs), trickle = as.double(0), vecTcalc = as.double(vecTcalc), nTimes = as.integer(nTimes), ndays = as.integer(ndays) , 
                    rtn_daily = as.double(arr_daily), rec_daily = as.double(arr_daily), 
                    hos_daily = as.double(arr_daily), nMCMC = as.integer(nMCMC), ithin = as.integer(ithin), tab = as.single(tab))
    
    tab = array(out$tab, c(nMCMC/ithin, nparamtot+1))
    colnames(tab) = c(parnames, 'LLK')
    parBest = out$par
  
  nReals = 1000
  
  cat('\n', state_abbv, 'parBest', parBest, '\n')  
  
  par_mcmc = tab[1:nReals, 1:nparamtot]
  for(ii in 1:nReals) par_mcmc[ii,1:nparamtot] = parBest
  
  
  # here ndays will be replaced by ndays + nforecast 
  # tps will be updated accordignly as will vecTcals, nTimes 
  
 
  arr_daily = arr_daily = array(0, c(ndays_total, nReals))
  vecTfrcst = seq(0, ndays_total,dt)
  nFrcst = length(vecTfrcst)
  
  y2 <- .Fortran("stosirh", vecN = as.integer(vecN), vecI0 = as.integer(vecI0), nparam = as.integer(nparam), nb = as.integer(nb), par_tab = as.double(par_mcmc), 
                 trickle = as.double(0), vecTfrcst = as.double(vecTfrcst), nReals = as.integer(nReals), nFrcst = as.integer(nFrcst), ndays_total = as.integer(ndays_total) , 
                 rtn_daily = as.integer(arr_daily), rec_daily = as.integer(arr_daily), 
                 hos_daily = as.integer(arr_daily))
  
  
  
  tps <- seq(from = 0, to = round(max(vecTfrcst)))
  
  noTPts <- length(tps)
  
  rtn_daily = array(y2$rtn_daily, c(ndays_total,nReals))
  rec_daily = array(y2$rec_daily, c(ndays_total,nReals))
  hos_daily = array(y2$hos_daily, c(ndays_total,nReals))
  
  hos_mean = hos_med =  hos_ll = hos_ul = hos_l80 = hos_u80 = rep(0, ndays_total)
  for (ii in 1:ndays_total) {
    hos_mean[ii] = mean(hos_daily[ii, ])
    hos_ll[ii] = quantile(hos_daily[ii, ], c(0.025))
    hos_ul[ii] = quantile(hos_daily[ii, ], c(0.975))
    hos_med[ii] = quantile(hos_daily[ii, ], c(0.5))
    hos_l80[ii] = quantile(hos_daily[ii, ], c(0.1))
    hos_u80[ii] = quantile(hos_daily[ii, ], c(0.9))  
  }
  
  data = data.frame(dates = dates, obs = obs)
  state_hosp_frcst = list(state = state_abbv, data = data, profiles = list(dates=dates_total, profiles = hos_daily), fit = TRUE)
    
  all_state_hosp_frcst[[istate]] = state_hosp_frcst
  
  # update the RData file with the information we have 
  
  ymax = max(hos_ul, obs)
  plot(dates_total, hos_mean, ylim = c(0, ymax), type = "l", lty = 2, col = "blue", lwd = 2, xlab = "", ylab = "Daily Hos. Admin.", xaxt='n')
  polygon(c(dates_total, rev(dates_total)), c(hos_ll, rev(hos_ul)), col = '#EBECF0', border = NA)
  polygon(c(dates_total, rev(dates_total)), c(hos_l80, rev(hos_u80)), col = '#D3D3D3', border = NA)
  lines(dates_total, hos_mean, type = "l", col = "blue", lwd = 2, lty = 2)
  lines(dates, hos_mean[1:ndays], type = 'l', lwd = 2, col = 'blue')
  lines(dates, obs, type = "h", col = "cornflowerblue", lwd = 2)
  abline(v=dates[ndays], lty = 2, col = 'coral')
  legend("topleft", c(state_abbv,"Reported", "Fit/forecast"), text.col = c("coral", "cornflowerblue",'blue'), bty = "n")
  legend_beta = round(parBest[(nparam+1)] * D_IC, digits = 2)
  for (jj in 2:nb) {
    legend_beta = c(legend_beta, round(parBest[(nparam+jj)] * D_IC, digits = 2))
  }
  legend('topright', legend = legend_beta, bty = 'n')
  axis.Date(1, at = seq(dates_total[1], dates_total[ndays_total]+6, 'weeks'), format = '%y-%m-%d') 
  axis.Date(1, at = seq(dates_total[1], dates_total[ndays_total]+6, 'days'), tcl = -0.2, labels = FALSE)
}

dev.off()

saveRDS(all_state_hosp_frcst, file = data_filename)
cat('\n\nSaved results to: ', data_filename, '\n')
cat('\n\nSaved plots to: ', plot_filename, '\n')

