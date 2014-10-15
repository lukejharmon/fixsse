sim.mk.jump <- function(phy, matrix, node, nsim=1, root=1)
{
	## This function simulates a mk model with a jump.
	
    ## node = node where to transition matrix[[1]] -> matrix[[2]]
    ## matrix = list with two Q matrices, or list of lists for several traits
	## nsim = number of character simulations. Will return output in geiger2 format.
	## root = value for the root state (values here are 1 & 2).

    model <- "discrete"
    ## if only one char
    if(is.matrix(matrix[[1]])){
        nchar <- 1
        model.matrix <- lapply(matrix, function(x) geiger:::.make.modelmatrix(x, model)[[1]])
        model.matrix <- list(model.matrix)
    }
    ## if more than one char (two level list object)
    if(is.list(matrix[[1]])){
        nchar <- length(matrix)
        model.matrix <- lapply(1:nchar, function(y) lapply(matrix[[y]], function(x) geiger:::.make.modelmatrix(x, model)[[1]]) )
    }
    
    if(length(root) > 1) stop("'root' should be a single value")

    ## Find the vector for 'jump.node'
    jump.node <- geiger:::.get.descendants.of.node(node, phy)
    jump.node <- append(jump.node, node)

    nbranches <- nrow(phy$edge)
    nspecies <- Ntip(phy)
    nnode <- Nnode(phy)
    rt <- nspecies+1
    zphy <- reorder.phylo(phy, "postorder")
    el <- zphy$edge.length
    result.tip <- array(0, dim=c(nspecies, nchar, nsim))
    result.node <- array(0, dim=c(nnode, nchar, nsim))

    .get.state <- function(s, p){
        pp <- cumsum(p[s,])
        min(which(runif(1)<pp))
    }
    
    for(j in 1:nchar) {
        mA <- model.matrix[[j]][[1]]
        mB <- model.matrix[[j]][[2]]
        if(!root%in%c(1:nrow(mA))) stop(paste("'root' must be a character state from 1 to ", nrow(m), sep=""))
        pA <- lapply(el, function(l) matexpo(mA*l))
        pB <- lapply(el, function(l) matexpo(mB*l))
        for(k in 1:nsim) {
            node.value <- numeric(nspecies+nnode)
            node.value[rt] <- root
            for(i in nbranches:1) {
                cur <- zphy$edge[i,2]
                anc <- zphy$edge[i,1]
                if(cur %in% jump.node){
                    curp <- pB[[i]] ## matrix B
                } else {
                    curp <- pA[[i]] ## matrix A
                }
                s <- node.value[anc]
                node.value[cur] <- .get.state(s, curp)
            }
            result.tip[,j,k] <- node.value[1:nspecies]
            result.node[,j,k] <- node.value[nspecies+1:nnode]
        }
    }
	## Return the tip states, node states and the reordered tree for plot.
    rownames(result.tip) <- zphy$tip.label
    rownames(result.node) <- nspecies+1:nnode
    out <- list(tip.state = result.tip, node.state = result.node, reorder.phylo = zphy)
    return(out)
}

select.node <- function(phy, percent){
    t.nd <- Nnode(phy)
    t.sp <- Ntip(phy)
    nodes <- t.sp+1:t.nd
    repeat
    {
        nds <- sample(nodes, 1)
        tps <- get.descendants(nds, phy, tips.only = TRUE)
        if( (t.sp * percent[2]) > length(tps) && length(tps) > (t.sp * percent[1]) )
            {
                break
            }
    }
    return(nds)
}

run.bisse <- function(tree, st, unres, tun.steps, chain.steps, constrain = "TRUE"){
	## Function to prepare prior, make MLE starting point estimate, tunning and running
	##		BiSSE mcmc.
    ## The returning object of this function will be a list with the mcmc run, the lik function and the
    ## 		prior distribution.
	
    ## tree = phylo
	## st = vector of states
	## unres = unresolved matrix (as in make.bisse)
	## steps = number of steps of the mcmc chain
    ## save = save every x steps
	## file = name of the file to save.
    
    tree <- multi2di(tree)
    mmm <- match(tree$tip.label, names(st))
    lik <- make.bisse(tree, st[mmm], unresolved = unres)
	
	if(constrain == "TRUE"){
		w.init <- rep(1,4)
		lik <- constrain(lik, lambda1~lambda0, mu1~mu0)
		print("Constrained model.")
	} else {
		w.init <- rep(1,6) 
		print("Full model.")
	}

	## Flag:
	print("Start analysis...")

    start <- starting.point.bisse(tree)
    prior <- make.prior.exponential(1 / 2 * (start[1] - start[3]))

    fit <- find.mle(lik, start[argnames(lik)])

	## Flag:
	print("ML estimate finished...")

    tun <- mcmc(lik, fit$par, nsteps = tun.steps, w = w.init, lower = 0, print.every = 0, prior = prior)

	if(constrain == "TRUE"){
		w <- diff(sapply(tun[2:5], range))
	} else {
		w <- diff(sapply(tun[2:7], range))
	}

	## Flag:
	print("MCMC tunning finished...")
	print("Starting MCMC chain...")

    run <- mcmc(lik, fit$par, nsteps = chain.steps, w = w, lower = 0, print.every = 0, prior = prior)

	## Flag:
	print("Done!")

    return(list(run, lik, prior))
}

to.range <- function(dens1, dens2){
	## Function to calculate range of densities for plotting results.
	## Only to be called by the 'plot_results.R' script.

    x.rg.unif <- lapply(index, FUN = function(x) range(dens1[[x]]$x) )
    x.rg.jump <- lapply(index, FUN = function(x) range(dens2[[x]]$x) )
    x.rg.unif <- do.call(rbind, x.rg.unif)
    x.rg.jump <- do.call(rbind, x.rg.jump)
    x.limit.unif <- c( min(x.rg.unif[,1]), max(x.rg.unif[,2]) )
    x.limit.jump <- c( min(x.rg.jump[,1]), max(x.rg.jump[,2]) )
    x.limit <- c( min(x.limit.unif[1],x.limit.jump[1]), max(x.limit.unif[2],x.limit.jump[2]) )

    y.rg.unif <- lapply(index, FUN = function(x) range(dens1[[x]]$y) )
    y.rg.jump <- lapply(index, FUN = function(x) range(dens2[[x]]$y) )
    y.rg.unif <- do.call(rbind, y.rg.unif)
    y.rg.jump <- do.call(rbind, y.rg.jump)
    y.limit.unif <- c( min(y.rg.unif[,1]), max(y.rg.unif[,2]) )
    y.limit.jump <- c( min(y.rg.jump[,1]), max(y.rg.jump[,2]) )
    y.limit <- c( min(y.limit.unif[1],y.limit.jump[1]), max(y.limit.unif[2],y.limit.jump[2]) )

    res <- list(xlim=x.limit, ylim=y.limit)
    return(res)
}
