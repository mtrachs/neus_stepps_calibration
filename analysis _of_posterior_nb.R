library(rstan)
library(RColorBrewer)
library(fields)
library(sp)
library(maptools)
#-----------------------------------------------------------------------------------------------
setwd('~/workflow_stepps_calibration/prediction/')
help.fun.loc <- 'utils/'
data.loc <- 'output_nb//'
plot.loc <- 'plots_nb/'#'plots_w_full/'
#-----------------------------------------------------------------------------------------------
us.shp <- readShapeLines('data/map_data/us_alb.shp',proj4string=CRS('+init=epsg:3175'))
#-----------------------------------------------------------------------------------------------
taxa <- sort(c('Maple','Birch','Hickory','Chestnut','Other conifer','Other hardwood',
               'Beech','Ash','Tamarack','Spruce','Pine','Poplar','Oak','Hemlock','Elm'))
K <- length(taxa)
prediction.files <- list.files(data.loc)

sapply(prediction.files,function(x){
  
  
  fit <- read_stan_csv(paste(data.loc,x,sep=''))
  post <- rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
  #dim(post)
  var_names<- colnames(post[,1,])
  par_names <- sapply(var_names, function(x) unlist(strsplit(x,'[[]'))[1])
  gauss.process <- grep('g',par_names)
  gauss.process <- gauss.process[-1]
  
  #---------------------------------------------------------------------------------------------------------------------
  post.exp <- exp(post[,1,gauss.process])
  #find indexes that start the fun
  ash.index <- seq(1,ncol(post.exp),K)
  
  summary.proportion <- 
    sapply(ash.index,function(z){
      exp.site <- cbind(post.exp[,z:(z+(K-1))]) # last one is adding exp(0) = 1 baseline taxon
      prop.site <- exp.site/rowSums(exp.site)
      prop.mean.site <- colMeans(prop.site)
      uncertainty <- apply(prop.site,2,function(x) quantile(x,probs=c(0.025,0.5,0.975)))
      data.frame(mean.proportion = prop.mean.site,median = uncertainty['50%',],lb = uncertainty['2.5%',],ub = uncertainty['97.5%',])
    })
  summary.proportion <- t(summary.proportion)
  
  mean.proportion <- matrix(ncol=K,unlist(summary.proportion[,'mean.proportion']),byrow=TRUE)
  colnames(mean.proportion) <- taxa
  
  
  # mean.proportion <- matrix(ncol=15,unlist(summary.proportion[,'mean.proportion']),byrow=TRUE)
  # colnames(mean.proportion) <- taxa
  
  
  median.proportion <- matrix(ncol=K,unlist(summary.proportion[,'median']),byrow=TRUE)
  colnames(median.proportion) <- taxa
  
  lb.proportion <- matrix(ncol=K,unlist(summary.proportion[,'lb']),byrow=TRUE)
  colnames(lb.proportion) <- taxa
  
  ub.proportion <- matrix(ncol=K,unlist(summary.proportion[,'ub']),byrow=TRUE)
  colnames(ub.proportion) <- taxa
  
  uncertainty <- ub.proportion - lb.proportion
  
  
  figure.name <- strsplit(x,'.csv')[[1]][1]
  
  
  #transfrom pollen coordinates to us coordinates
  pollen_input <- load('~/workflow_stepps_calibration/calibration/data/elicitation_neus_certainty_median_only_abies.RData')
  sputm <- SpatialPoints(pollen_coords, proj4string=CRS("+init=epsg:4326"))
  spgeo <- spTransform(sputm, CRS("+init=epsg:3175"))
  pollen_coord_us <- spgeo@coords
  
  x1 <- strsplit(x,'.csv')[[1]][1]
  
  pdf(paste(plot.loc,x1,'.pdf',sep=''),width = 10,height = 5)
  par(mfrow=c(1,2),oma=c(1,1,1,2),cex = 0.8,cex.axis = 0.8)
  sapply(taxa, function(y){
    
    
    breaks <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,1)
    categories <- cut(mean.proportion[,y],breaks,include.lowest = TRUE,labels = FALSE)
    
    colours <- rev(brewer.pal(10,'RdYlBu'))
    colours.plot <- colours[categories]
    
    east <- sort(unique(veg_coords[,'meters.east']))
    north <- sort(unique(veg_coords[,'meters.north']))
    
    breaks1 <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7)
    breaks2 <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,NA)
    
    plot(us.shp,xlim=range(east),ylim=range(north),main = paste(y,'STEPPS'))
    image.plot(east, north,z = matrix(ncol= length(north),rep(-1,length(east)*length(north))),zlim=c(0,0.7),col = colours,
               breaks= breaks1,lab.breaks = breaks2,cex.axis = 0.8,add=TRUE)
    points(veg_coords,col=colours.plot,pch = 15)
    points(pollen_coord_us,col=1,pch = 16,cex = 0.75)
    plot(us.shp,add=TRUE)
    
    #----------------------------------------------------------------------------------------------------------------------
    veg <- as.matrix(r)
    categories <- cut(veg[,y],breaks,include.lowest = TRUE,labels = FALSE)
    colours.plot <- colours[categories]
    
    
    plot(us.shp,xlim=range(east),ylim=range(north),main = paste(y,'Paciorek et al. (2016)'))  
    image.plot(east, north,z = matrix(ncol= length(north),rep(-1,length(east)*length(north))),zlim=c(0,0.7),col = colours,
               breaks= breaks1,lab.breaks = breaks2,cex.axis = 0.8,add = TRUE)
    points(veg_coords,col=colours.plot,pch = 15)
    points(pollen_coord_us,col=1,pch = 16,cex = 0.75)
    plot(us.shp,add=TRUE)
  })
  dev.off()
})


