## This is the simulation of a Mk model with a jump and check for BiSSE behaviour on parameter estimates.
## This script create the data, runs the simulations and estimate the parameters in BiSSE.
## Use the plot results script for plotting.

library(geiger)
library(diversitree)
library(multicore)

source("./sim_mk_functios.R")

## Get data for simulation:
tr <- lapply(1:10, function(x) sim.bdtree(stop = "taxa", n = 100))
node <- sapply(tr, function(x) select.node(x, percent = c(0.2, 0.5)))

## Set values for the single and jump MK models:
q1 <- rbind(c(-.2, .2), c(.5, -.5)); q1
q2 <- rbind(c(-.5, .5), c(.2, -.2)); q2
mm <- list(q1,q2)

## Make the jumpMK model simulation:
mk.jump <- lapply(1:10, function(x) sim.mk.jump(tr[[x]], matrix = mm, node = node[x], root = 1)$tip.state[,,1] )
mk.unif <- lapply(1:10, function(x) sim.char(tr[[x]], par = mm[[1]], model = "discrete", root = 1)[,,1] )

run.sim <- function(phy, state){
    state[which(state == "1")] <- 0
    state[which(state == "2")] <- 1
    res <- run.bisse(phy, state, unres = NULL, tun.steps = 100, chain.steps = 10000, constrain = FALSE)
    return(res)
}

## Run the parameter estimates in parallel:
index <- 1:10

tasks <- list(
    job1 <- function() lapply(index, FUN = function(x) run.sim(phy = tr[[x]], state = mk.jump[[x]])),
    job2 <- function() lapply(index, FUN = function(x) run.sim(phy = tr[[x]], state = mk.unif[[x]]))
    )

out <- mclapply(tasks, function(f) f(), mc.cores = 2)

save.image("sim_mk.RData")
