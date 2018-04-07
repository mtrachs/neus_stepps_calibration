library(rstan)
library(RColorBrewer)
library(fields)
library(sp)
library(maptools)
#-----------------------------------------------------------------------------------------------
setwd('~/workflow_stepps_calibration/prediction/')
data.loc <- 'output_low_resolution/'
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
taxa <- sort(c('Maple','Birch','Hickory','Chestnut','Other conifer','Other hardwood',
               'Beech','Ash','Tamarack','Spruce','Pine','Poplar','Oak','Hemlock','Elm'))





prediction.files <- list.files(data.loc)
x <- prediction.files

K = 15
N= 754
T=1



fit <- read_stan_csv(paste(data.loc,x,sep=''))
post <- rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
#dim(post)
var_names<- colnames(post[,1,])
par_names <- sapply(var_names, function(x) unlist(strsplit(x,'[[]'))[1])
#-----------------
gauss.process <- grep('g',par_names)
gauss.process <- gauss.process[-1]

post.exp <- exp(post[,1,gauss.process])
#post.exp <- post.exp[(nrow(post.exp)-199):nrow(post.exp),] 
#find indexes that start the fun
ash.index <- seq(1,ncol(post.exp),K)

summary.proportion <- 
  sapply(ash.index,function(z){
    exp.site <- cbind(post.exp[,z:(z+K-1)]) # last one is adding exp(0) = 1 baseline taxon
    prop.site <- exp.site/rowSums(exp.site)
    prop.mean.site <- colMeans(prop.site)
    uncertainty <- apply(prop.site,2,function(x) quantile(x,probs=c(0.025,0.5,0.975)))
    data.frame(mean.proportion = prop.mean.site,median = uncertainty['50%',],lb = uncertainty['2.5%',],ub = uncertainty['97.5%',])
  })
summary.proportion <- t(summary.proportion)

mean.proportion <- matrix(ncol=K,unlist(summary.proportion[,'mean.proportion']),byrow=TRUE)
colnames(mean.proportion) <- taxa

#------------------------------------------------------------------------------------------------------------------
#use Andrias code 
setwd('~/github_clones/stepps-prediction/')
source('r/utils/pred_plot_funs.r')
source('r/utils/pred_helper_funs.r')
source('r/utils/build_mu_g.r')
source('r/read_stanbin.r')
source('r/')




post_dat = list(post = post, par_names = par_names)


process_out = build_r_nb(post_dat, N, T, K)

dim(process_out$r)
site_means <- 
  sapply(1:N,function(x){
  rowMeans(process_out$r[x,,])
})    
site_means <- t(site_means)
colnames(site_means) <- colnames(mean.proportion)

#test if equal
all.equal(mean.proportion,site_means)
#there we go!