#--------------------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------------------
#devtools::install_github("PalEON-Project/stepps-cal")
library(stepps)
library(dplyr)
library(DT)
library(neotoma)
library(sp)
library(fields)
library(rgdal)
library(abind)
library(rstan)

#-----------------------------------------------------------------------------------------------
setwd('~/workflow_stepps_calibration/prediction/')
help.fun.loc <- 'utils/'
data.loc <- 'data/'
plot.loc <- 'plots/'

source(paste(help.fun.loc,'pred_helper_funs.r',sep='')) #make sure I use medians of a
source(paste(help.fun.loc,'process_funs.r',sep='')) #make sure I use medians of a



#load the data used for the vegetation model (we have to use the same knots)
load('~/workflow_stepps_calibration/vegetation/data/veg_data_15_taxa_6796_cells_260_knots.rdata')

fit <- read_stan_csv('~/workflow_stepps_calibration/vegetation/output/veg_paciorek_neus_nb.csv')
post <- rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
param.names <-colnames(post[,1,]) #find parameter names of posterior
param.greek <- sapply(param.names,function(x) strsplit(x,'[[]')[[1]][1]) #only take greek letter
param.e.bayes <- param.greek%in%c('eta','rho') #only take eta and rho
est.e.bayes <- colMeans(post[,1,param.e.bayes])

eta <- est.e.bayes[grep('eta',names(est.e.bayes))]
rho <- est.e.bayes[grep('rho',names(est.e.bayes))]

#-----------------------------------------------------------------------------------------------
#produce 
#-----------------------------------------------------------------------------------------------

for(num_sites in c(88,62,73)){
#load pollen data
#88 sites including sites from the edge of the domain 

  if(num_sites==88) {
    load('~/workflow_stepps_prediction/prediction/data/elicitation_neus_certainty_median_only_abies_88_sites.RData')
  }
  #-----------
  if(num_sites==62) {
    load('~/workflow_stepps_calibration/calibration/data/elicitation_neus_certainty_median_only_abies.RData')
  }
  #73 sites no sites on the boundaries and not sites 
  if(num_sites==73) {
    load('~/workflow_stepps_prediction/prediction/data/elicitation_neus_certainty_median_only_abies_73_sites.RData')
  }
  
  

#pollen_coords
# K
# idx_cores
# N_cores
# d
# y
  N <- N_cells

#------------------------------------------------------------------------------------------------
#produce knots
#------------------------------------------------------------------------------------------------
#not needed here because knots and distance matrices are loaded from the data used for STEPPS-vegetation
# clust_n <- 260
# knot_coords = kmeans(veg_coords, clust_n)$centers 
# knot_coords = unname(knot_coords)
# N_knots     = dim(knot_coords)[1]
# 
# #-------------------------------------------------------------------------------------------------------------------------------
# # new d_inter
# #-------------------------------------------------------------------------------------------------------------------------------
# d_inter <- rdist(veg_coords,knot_coords)/10^6
# d_knots <- rdist(knot_coords,knot_coords)/10^6
#--------------------------------------------------------------------------------------------
# have to load parameters from calibration model
#--------------------------------------------------------------------------------------------
  

  source('R/build_cal_main.r') # this is strange....
  for (run in runs) {
    kernel <- run$kernel
    num_a <- run$one_a
    one_psi <- run$one_psi
    handle <- strsplit(run$suff_fit,'_A')[[1]][1]


# loading results of the calibration model is a bit messy. I ran STEPPS calibration for 200 samples 
# for 10 parallel chains
    
    for (i in 0:9) { 
      if(kernel=='pl'){
        if(num_a==FALSE){
          fname = paste('~/stepps_median/output_cal_pl_Ka_Kgamma_EPs_',i,'.csv',sep='') #sprintf('%s/%s.csv', path_out, suff_fit)
        }
        if((num_a==TRUE)&(i<5)){
          fname = paste('~/workflow_stepps_calibration/results/data/stepps_median/output_cal_pl_',i,'.csv',sep='') #sprintf('%s/%s.csv', path_out, suff_fit)
        }
        if(run$a_const==TRUE){
          fname = paste('~/stepps_neus_abies/output/output_cal_pl_Ka_Kgamma_EPs_',i,'.csv',sep='') #sprintf('%s/%s.csv', path_out, suff_fit)
        }  
      }
      if(kernel=='gaussian'){
        if(one_psi==FALSE){
          fname = paste('~/stepps_median/output_cal_g_Kpsi_Kgamma_EPs_',i,'.csv',sep='') #sprintf('%s/%s.csv', path_out, suff_fit)
        }
        if(one_psi==TRUE){
          fname = paste('~/stepps_median/output_neus_',i,'.csv',sep='') #sprintf('%s/%s.csv', path_out, suff_fit)
        }
      }
      fit <- read_stan_csv(fname)
      if (i==0)  {post = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)}
      else{post_new <-  rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
      post <- abind(post,post_new,along=1)
    }
  }

  param.names.cal <-colnames(post[,1,]) #find parameter names of posterior
  param.greek.cal <- sapply(param.names.cal,function(x) strsplit(x,'[[]')[[1]][1]) #only take greek letter
  param.e.bayes.cal <- param.greek.cal%in%c('gamma','phi') #only take gamma and phi
  est.e.bayes.cal <- colMeans(post[,1,param.e.bayes])
  gamma <- est.e.bayes.cal[grep('gamma',names(est.e.bayes.cal))]
  phi <- est.e.bayes.cal[grep('phi',names(est.e.bayes.cal))]

  w <- build_weight_matrix(post = post,d = t(d),idx_cores = idx_cores,
                         N = N_cells,N_cores =  N_cores,run = run)


  sum_w_pot <- build_sumw_pot(post = post, K=K,N_pot = N_pot,
                            d_pot =  d_pot,run = run)

#-------------------------------------------------------------------------------------------------------------
# prediction model requires 
  if(length(gamma)==1) {
    gamma <- rep(gamma,K)
    w1 <- array(dim=c(K,N_cores,N_cells))
    for(i in 1:K) {w1[i,,]<-w}
    w <- w1
  }
#-------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
  res <- 1 # 
  T <- 1
#-------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#save data in .dump file
  stan_rdump(
    c('K', 'N','T', 'N_knots', 'N_cores',
      'y', 'res',
      'sum_w_pot',
      'rho','eta','gamma','phi',
      'idx_cores',
      'd','d_knots','d_inter','w'
      ), 
  file=paste('~/workflow_stepps_calibration/vegetation/data_nb/prediction_',K,'_taxa_',N_cells,'_cells_',N_knots,'_knots_',handle,'_',num_sites,'_sites_only_abies.dump',sep=""))


  save(K, N,T, N_knots, N_cores,
       y, res,
       sum_w_pot,
       rho,eta,gamma,phi,
       idx_cores,
       d,d_knots,d_inter,w,
      pollen_coords,r,veg_coords, 
  file=paste('~/workflow_stepps_calibration/vegetation/data_nb/prediction_',K,'_taxa_',N_cells,'_cells_',N_knots,'_knots_',handle,'_',num_sites,'_sites_only_abies.rdata',sep=""))
  }
}
