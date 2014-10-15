library(geiger)

## Simulations rvMK: Heterogeneous variation of the MK in the tree.
## I am using the 'geiger' sim.char as a base to do the simulation with two mk models.

source("./sim_mk_functions.R")

## Get a tree and a list with Q matrices:
tr <- sim.bdtree(stop = "taxa", n = 100, seed = 1212)
q1 <- rbind(c(-.0, .0), c(1, -1)); q1
q2 <- rbind(c(-2, 2), c(.0, -.0)); q2
mm <- list(q1,q2)
mm

## Select a node to change the matrix:
plot(tr, direction="upward"); nodelabels(c("Selected \n 142","Root state \n 1"), c(142,101))

## Run the mk.jump simulation:
mk.sim <- sim.mk.jump(tr, matrix = mm, node = 142, root = 1)

## What are the simulated states?
tra <- mk.sim[[3]]
tra$tip.label <- mk.sim[[1]][,,1]
plot(tra, direction="upward")
txt <- as.character(mk.sim[[2]][,,1])
nds <- as.numeric(names(mk.sim[[2]][,,1]))
col <- rep("black", times = length(mk.sim[[2]][,,1]) )
col[which(nds == "142")] <- "red"
nodelabels(text = txt, node = nds, col = col)
