# clear R's brain
rm(list = ls())


# packages

library(ggplot2)
library(dplyr)
library(reshape2)
library(maps)
library(ape)
library(phytools)
library(viridisLite)
library(TreeTools)



setwd('/home/habich/AH112/')

tree_all <- read.newick('forEff_with_outgroup_Ppar_HKY.treefile')
strains_in_tree <- c(TipLabels(tree_all))

ES <- read.csv('all_effector_sets_single_eff.csv', header = T)
rownames(ES) <- ES$Sequence_ID
strains_in_ES_file <- c(ES$Sequence_ID)


test_all <- midpoint.root(test_all)

all_acc_eff <- c('PA0260','PALA36_00282','PA0262','PA0822','PA3487','PA14-43100', 'PAKAF_03765','PALA52_00191',
                 'PA0093_tse6', 'PA0093_tas1', 'PA0093_tse6a',
                 'PA0099_tse7', 'PA0099_tse7a', 'PA0099_tse7b', 'PA0099_tse7c', 'PA0099_tse7d', 'PA0099_tse7e',
                 'PA5265_tspE1b', 'PA5265_tspE1c',
                 'PA14_33970_tepB','PA14_33970_tepBa','PA14_33970_tepBb')


for(locus in all_acc_eff) {
  #locus <- 'PALA52_00191'
  
  eff <- ES[colnames(ES) %in% locus]
  eff <- setNames(eff[,1], rownames(ES))
  eff <- as.factor(eff)
  levels(eff)
  
  eff.ER_model <- fitMk(test_tree, eff, model = 'ER')
  eff.ARD_model <- fitMk(test_tree, eff, model = 'ARD')
  eff.Irr1_model <- fitMk(test_tree, eff, model = matrix(c(0,1,0,0),2,2,byrow=TRUE))
  eff.Irr2_model <- fitMk(test_tree, eff, model = matrix(c(0,0,1,0),2,2,byrow=TRUE))
  
  eff.aov <- anova(eff.ER_model, eff.ARD_model, eff.Irr1_model, eff.Irr2_model)
  capture.output(eff.aov, file = paste0('results/', locus, '_aov.txt'), row.names = F)
 
  eff.simmap <- simmap(eff.aov, nsim=1000)
  summary(eff.simmap)
  capture.output(summary(eff.simmap), file = paste0('results/', locus, '_simmap.txt'), row.names = F)

  cols <- setNames(viridisLite::viridis(n=length(levels(eff))), levels(eff))
  
  #plot posterior probabilites: pie charts
  pdf(paste0('results/', locus, '_piecharts.pdf'), width = 25, height = 400)
  
  plot(summary(eff.simmap), ftype='i', fsize=0.75,
       colors=cols, cex=c(0.15,0.05))
  legend('topright', levels(eff),pt.cex=1.5,pt.bg=cols, bty='n',cex=0.8,fill=cols)
  dev.off()
  
  
  eff.densityMap <- densityMap(eff.simmap,plot=FALSE, res=1000)
  
  #plot posterior probabilites: gradient
  pdf(paste0('results/', locus, '_gradient.pdf'), width = 25, height = 400)
  eff.densityMap <- setMap(eff.densityMap, viridisLite::viridis(n=10))
  plot(eff.densityMap,lwd=3,outline=TRUE,fsize=c(0.6,0.7),legend=0.1)
  dev.off()
  
  eff.density <- density(eff.simmap)
  eff.density
  capture.output(eff.density, file = paste0('results/', locus, '_density.txt'), row.names = F)
  pdf(paste0('results/', locus, '_density.pdf'), width = 10, height = 5) 

  par(mfrow=c(1,2),las=1,cex.axis=0.7,cex.lab=0.8)
  COLS<-setNames(cols[2:1],eff.density$trans)
  plot(eff.density,ylim=c(0,0.6),
       transition=names(COLS)[1],colors=COLS[1],main="")
  mtext('A',line=1,adj=0,
        cex=0.8)
  plot(eff.density,ylim=c(0,0.6),
       transition=names(COLS)[2],colors=COLS[2],main="")
  mtext('B',line=1,adj=0,
        cex=0.8)
  dev.off()

}
