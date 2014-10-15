## This script will do some plots to show how the simulation deviates
## 		from the real value of the parameters.
## Need to do some frequency histogramns to show the parameter values
## 		however, all the distributions need to be in the same plot
## 		for adequate comparison. Use the 'add' tool to help with this.

library(geiger)
library(diversitree)
load("./sim_mk.RData")
source("./sim_mk_functions.R")

real <- list("lambda0" = 1,"lambda1" = 1,"mu0" = 0,"mu1" = 0, "q01" = c(0.2,0.5), "q10" = c(0.5,0.2))
index <- 1:10

mcmc.unif <- lapply(index, FUN = function(x) out[[1]][[x]][[1]])
mcmc.jump <- lapply(index, FUN = function(x) out[[2]][[x]][[1]])

nm <- names(mcmc.unif[[1]][,-c(1,dim(mcmc.unif[[1]])[2])])

for(i in nm){
    pos <- which(names(mcmc.unif[[1]]) == i)

    dd.unif <- lapply(index, FUN = function(x) density(mcmc.unif[[x]][,pos]) )
    dd.jump <- lapply(index, FUN = function(x) density(mcmc.jump[[x]][,pos]) )
    
    limits <- to.range(dd.unif, dd.jump)
    
    jpeg(paste("res_",i,".jpeg",sep = ""), width = 600, height = 600, quality = 90)
    plot(dd.unif[[1]], xlim = limits[[1]], ylim = limits[[2]], col = "blue", main = i)
    for(i in index){
        lines(dd.unif[[i]], col = "blue")
        lines(dd.jump[[i]], col = "red")
    }
    abline(v = real[[pos-1]], col = c("blue","red"), lty = 3, lwd = 1.5)
    dev.off()
}
