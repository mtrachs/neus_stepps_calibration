# library(sp)
# library(rgdal)
library(maptools)

#-------------------------------------------------------------------------------------------------------------------
setwd('~/workflow_stepps_calibration/calibration/')
help.fun.loc <- 'calibration_helper_funs/'
data.loc <- 'data/'
plot.loc <- '~/workflow_stepps_calibration/vegetation/plots/'
source('utils/dataPlotFuns.r')
source('utils/simDataFuns.r')

#have to change that one
#load('data/elicitation_neus_certainty_median.RData')
for(num_sites in c(88,62,73)){
  #load pollen data
  #88 sites including sites from the edge of the domain 
  
load(paste('~/workflow_stepps_calibration/vegetation/data_nb/prediction_15_taxa_6796_cells_260_knots_cal_g_',num_sites,'_sites_only_abies.rdata',sep=''))
  
########################################################################################################
## pls pie map
########################################################################################################
# reorder data by pollen abundance
new.order = order(colSums(y), decreasing=TRUE)
taxa <- colnames(y)

taxa = taxa[new.order]
y = y[,new.order]


r = r[,new.order]
colnames(r) = taxa

col_list = c("#1F78B4", "#33A02C", "#E31A1C", "#6A3D9A", "#B15928", "#FF7F00", 
             "#FFFF99", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6",'black','grey','white')

stepps_cols = data.frame(taxa=taxa, cols=col_list)

########################################################################################################
## pls pie map
########################################################################################################

ntaxa = length(taxa)
#what is this
centers   = veg_coords#centers_veg
colnames(centers) = c('x', 'y')

xlo = min(centers[,1])
xhi = max(centers[,1])
ylo = min(centers[,2])
yhi = max(centers[,2])

# subgrid = regular_subgrid(centers, dx=34000,dy=34000)
# knots_in = knots_in_domain4(subgrid, centers, cell_width = 8000)
# 
# plot(centers[,1], centers[,2], asp=1)
# points(knots_in[,1], knots_in[,2], col='red')
# 
# d = rdist(knots_in, centers)
# 
# pls_coarse = matrix(0, nrow=nrow(d), ncol=ntaxa)
# colnames(pls_coarse) = taxa
# 
# for (i in 1: nrow(centers)){
#   #print(min(d[,i]))
#   close_knot = which.min(d[,i])
#   
#   #   print(as.double(pls_cut[i,]))
#   #   print(pls_cut[i,])
#   
#   
#   pls_coarse[close_knot,] <- pls_coarse[close_knot,] + as.matrix(r[i,])
# }
# 
# pls_props  = t(apply(pls_coarse, 1, function(x) if (sum(x) != 0){x/sum(x)} else {x}))

shift=30000


#setwd('C:/Work/Paleon/STEPPS/stepps-calibration-master/')
# postscript('r/data/figs/pie_plot_pls_UMW_v0.2.eps', width=8, height=6)
# pdf(paste(plot.loc,'pie_plot_veg_NEUS.pdf',sep=''), width=12, height=10)
# par(mfrow=c(1,1))
# pieMap(proportions = pls_props, 
#        centers  = knots_in,
#        restrict = FALSE,
#        inputRestricted = FALSE,
#        xlim   = c(xlo+shift, xhi-shift),
#        ylim   = c(ylo+shift, yhi-shift),
#        radius = 14000,
#        scale  = 1,
#        xlab   = 'x',
#        ylab   = 'y',
#        add_legend = TRUE, 
#        main_title='',
#        col_list= col_list)
# dev.off()


########################################################################################################
## pollen pie map
########################################################################################################

pollen = y
pollen_props  = t(apply(pollen, 1, function(x) x/sum(x)))
colnames(pollen_props) = taxa


sputm <- SpatialPoints(pollen_coords, proj4string=CRS('+init=epsg:4326'))
spgeo <- spTransform(sputm, CRS("+init=epsg:3175"))


centers   = spgeo@coords#centers_polA
colnames(centers) = c('x', 'y')

pdf(paste(plot.loc,'pie_plot_pollen_NEUS_',num_sites,'_sites_only_abies.pdf',sep=''), width=12, height=10)
par(mfrow=c(1,1))
pieMap(proportions = pollen_props, 
       centers  = centers,
       restrict = FALSE,
       inputRestricted = FALSE,
       xlim   = c(xlo+shift, xhi-shift),
       ylim   = c(ylo+shift, yhi-shift),
       radius = 18000,
       scale  = 1,
       xlab   = 'x',
       ylab   = 'y', 
       add_legend=TRUE,
       main_title='',
       col_list=col_list)
dev.off()
}
# pdf('C:/Work/Paleon/STEPPS/NEUS/pie_plot_pollen_NEUS_ensemble_new.pdf', width=12, height=10)
# 
#   sapply(1:100,function(zzz) {
# 
#     load(paste('C:/Work/Paleon/STEPPS/NEUS/data_for_stan/elicitation_neus_',zzz,'.RData',sep=''))
# 
#     y = y[,new.order]
#     pollen = y
#     pollen_props  = t(apply(pollen, 1, function(x) x/sum(x)))
#     colnames(pollen_props) = taxa
#     sputm <- SpatialPoints(pollen_coords, proj4string=CRS('+init=epsg:4326'))
#     spgeo <- spTransform(sputm, CRS("+init=epsg:3175"))
#     centers   = spgeo@coords#centers_polA
#     colnames(centers) = c('x', 'y')
# 
# 
# 
#   par(mfrow=c(1,1))
#   pieMap(proportions = pollen_props, 
#        centers  = centers,
#        restrict = FALSE,
#        inputRestricted = FALSE,
#        xlim   = c(xlo+shift, xhi-shift),
#        ylim   = c(ylo+shift, yhi-shift),
#        radius = 18000,
#        scale  = 1,
#        xlab   = 'x',
#        ylab   = 'y', 
#        add_legend=TRUE,
#        main_title='',
#        col_list=col_list)
#   })
# dev.off()
