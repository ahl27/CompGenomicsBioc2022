to.dendrogram <- function(object, states=NULL, p=NULL, s=NULL) {
	z <- list()
	oHgts <- object$lengths
	oHgt <- object$height
	nMerge <- length(oHgt)
	hMax <- oHgt[nMerge]
	
	one <- 1L
	two <- 2L
	for (k in seq_len(nMerge)) {
		x <- as.integer(object$merge[k, ])
		neg <- x < 0
		if (all(neg)) { # two leaves
			zk <- as.list(-x)
			attr(zk, "members") <- two
			objlabels <- object$labels[-x]
			attr(zk[[1L]], "label") <- objlabels[1L]
			attr(zk[[2L]], "label") <- objlabels[2L]
			attr(zk[[1L]], "members") <- attr(zk[[2L]], "members") <- one
			attr(zk[[1L]], "height") <- oHgt[k] - oHgts[k, 1]
			attr(zk[[2L]], "height") <- oHgt[k] - oHgts[k, 2]
			attr(zk[[1L]], "leaf") <- attr(zk[[2L]], "leaf") <- TRUE
		} else if (any(neg)) { # one leaf, one node
			X <- as.character(x)
			isL <- x[1L] < 0 # is leaf left?
			if (isL) {
				zk <- list(-x[1L], z[[X[2L]]])
			} else {
				zk <- list(z[[X[1L]]], -x[2L])
			}
			attr(zk, "members") <- attr(z[[X[1 + isL]]], "members") + one
			attr(zk[[2 - isL]], "members") <- one
			attr(zk[[2 - isL]], "height") <- oHgt[k] - oHgts[k, 2 - isL]
			attr(zk[[2 - isL]], "label") <- object$labels[-x[2 - isL]]
			attr(zk[[2 - isL]], "leaf") <- TRUE
		} else { # two nodes
			x <- as.character(x)
			zk <- list(z[[x[1L]]], z[[x[2L]]])
			attr(zk, "members") <- attr(z[[x[1L]]], "members") + attr(z[[x[2L]]], "members")
		}
		attr(zk, "height") <- oHgt[k]
		attr(zk, "state") <- states[k]
		attr(zk, "probability") <- p[k]
		attr(zk, "support") <- s[k]
		k <- as.character(k)
		z[[k]] <- zk
	}
	z <- z[[k]]
	attr(z, "state") <- states[nMerge]
	attr(z, "probability") <- NULL
	attr(z, "support") <- NULL
	class(z) <- "dendrogram"
	z
}

.collapse <- function(dend, collapse, dim) {
	# initialize a stack of maximum length (dim)
	stack <- vector("list", dim)
	visit <- logical(dim) # node already visited
	parent <- integer(dim) # index of parent node
	index <- integer(dim) # index in parent node
	pos <- 1L # current position in the stack
	stack[[pos]] <- dend
	Score <- attr(dend, "score")
	LnLs <- attr(dend, "siteLnLs")
	while (pos > 0L) { # more nodes to visit
		if (visit[pos]) { # ascending tree
			visit[pos] <- FALSE # reset visit
			
			if (!is.leaf(stack[[pos]][[1]])) {
				h1 <- attr(stack[[pos]][[1]], "height")
			} else {
				h1 <- -Inf
			}
			if (!is.leaf(stack[[pos]][[2]])) {
				h2 <- attr(stack[[pos]][[2]], "height")
			} else {
				h2 <- -Inf
			}
			
			h <- attr(stack[[pos]], "height")
			
			if ((h - h1) <= collapse || (h - h2) <= collapse) {
				# make multifurcating
				m1 <- attr(stack[[pos]][[1L]], "members")
				m2 <- attr(stack[[pos]][[2L]], "members")
				states <- c(attr(stack[[pos]][[1L]], "state"),
					attr(stack[[pos]][[2L]], "state"))
				probs <- c(attr(stack[[pos]][[1L]], "probability"),
					attr(stack[[pos]][[2L]], "probability"))
				support <- c(attr(stack[[pos]][[1L]], "support"),
					attr(stack[[pos]][[2L]], "support"))
				m <- m1 + m2
				if ((h - h1) <= collapse && (h - h2) <= collapse) {
					l1 <- length(stack[[pos]][[1L]])
					l2 <- length(stack[[pos]][[2L]])
					x <- vector("list", l1 + l2)
					x[seq_len(l1)] <- stack[[pos]][[1L]][seq_len(l1)]
					x[seq_len(l2) + l1] <- stack[[pos]][[2L]][seq_len(l2)]
				} else if ((h - h1) <= collapse) {
					l <- length(stack[[pos]][[1L]])
					x <- vector("list", l + 1L)
					x[seq_len(l)] <- stack[[pos]][[1L]][seq_len(l)]
					x[l + 1L] <- stack[[pos]][-1L]
				} else if ((h - h2) <= collapse) {
					l <- length(stack[[pos]][[2L]])
					x <- vector("list", l + 1L)
					x[1L] <- stack[[pos]][-2L]
					x[seq_len(l) + 1L] <- stack[[pos]][[2L]][seq_len(l)]
				}
				stack[[pos]] <- x
				
				attr(stack[[pos]], "height") <- h
				attr(stack[[pos]], "members") <- m
				attr(stack[[pos]], "state") <- unique(states)
				if (length(probs) > 0)
					attr(stack[[pos]], "probability") <- min(probs)
				if (length(support) > 0)
					attr(stack[[pos]], "support") <- min(support)
				
				class(stack[[pos]]) <- "dendrogram"
			}
			
			# replace self in parent
			if (parent[pos] > 0)
				stack[[parent[pos]]][[index[pos]]] <- stack[[pos]]
			pos <- pos - 1L # pop off of stack
		} else { # descending tree
			visit[pos] <- TRUE
			p <- pos
			for (i in seq_along(stack[[p]])) {
				if (!is.leaf(stack[[p]][[i]])) {
					# push subtree onto stack
					pos <- pos + 1L
					stack[[pos]] <- stack[[p]][[i]]
					parent[[pos]] <- p
					index[[pos]] <- i
				}
			}
		}
	}
	
	dend <- stack[[1L]]
	attr(dend, "score") <- Score
	attr(dend, "siteLnLs") <- LnLs
	
	return(dend)
}

.organizeClusters <- function(myClusters,
	dNames,
	o) {
	l <- length(dNames)
	clusters <- data.frame(cluster=integer(l),
		row.names=dNames)
	w <- which(myClusters[, 7] < 0)
	if (length(w) > 0)
		clusters$cluster[-1*myClusters[w, 7]] <- as.integer(myClusters[w, 9])
	w <- which(myClusters[, 8] < 0)
	if (length(w) > 0)
		clusters$cluster[-1*myClusters[w, 8]] <- as.integer(myClusters[w, 10])
	
	# order the cluster numbers to match
	# the order of the dendrogram
	temp <- 0
	l <- max(clusters$cluster)
	clustersTemp <- clusters
	v <- vector(mode="numeric",length=l)
	j <- 0
	for (i in 1:length(o)) {
		if (clusters$cluster[o[i]] != temp &
			length(which(v==clusters$cluster[o[i]]))==0) {
			temp <- clusters$cluster[o[i]]
			j <- j + 1
			v[j] <- temp
		}
	}
	for (k in 1:l) {
		w <- which(clusters$cluster == v[k])
		clustersTemp$cluster[w] <- k
	}
	clusters <- clustersTemp
	
	return(clusters)
}

.organizeClustersFast <- function(myClusters,
	dNames) {
	l <- length(dNames)
	clusters <- data.frame(cluster=integer(l),
		row.names=dNames)
	w <- which(myClusters[, 7] < 0)
	if (length(w) > 0)
		clusters$cluster[-1*myClusters[w, 7]] <- as.integer(myClusters[w, 9])
	w <- which(myClusters[, 8] < 0)
	if (length(w) > 0)
		clusters$cluster[-1*myClusters[w, 8]] <- as.integer(myClusters[w, 10])
	
	return(clusters)
}

.rates1 <- function(alpha, nBins) {
	# Determine rates based on alpha and the number of bins
	# bins roots normalized to 1 of the Laguerre quadrature
	# first nBins elements are rates with mean 1
	# second nBins elements are probabilities with sum 1
	
	findRoots <- function(alpha, nBins) {
		
		# Determine rates based on Gamma's alpha and the number of bins
		# bins roots normalized to 1 of the General Laguerre Polynomial (GLP)
		
		coeff  <- integer(nBins + 1)
		for (i in 0:nBins) {
			n <- nBins + alpha
			k <- nBins - i
			coeff[i + 1] <- (-1)^i*choose(nBins + alpha, nBins - i)/factorial(i)
		}
		
		return(sort(Re(polyroot(coeff))))
	}
	
	roots <- findRoots(alpha - 1, nBins)
	
	Laguerre <- function(x, alpha, degree) {
		y <- 0
		for (i in 0:degree) {
			y <- y + (-1)^i*choose(degree + alpha, degree - i)*x^i/factorial(i)
		}
		return(y)
	}
	
	weights <- numeric(nBins)
	f <- prod(1 + (alpha - 1)/(1:nBins))
	
	for (i in 1:nBins) {
		weights[i] <- f*roots[i]/((nBins + 1)^2*Laguerre(roots[i],
			alpha - 1,
			nBins + 1)^2)
	}
	
	roots <- roots/alpha
	
	return(c(roots, weights))
}

.rates2 <- function(alpha, nBins) {
	if (nBins == 1)
		return(c(1, 1))
	quants <- qgamma((1:(nBins - 1))/nBins, shape = alpha, rate = alpha)
	c(diff(c(0, pgamma(quants * alpha, alpha + 1), 1)) * nBins,
		rep(1/nBins, nBins))
}

.optimizeModelDNA <- function(myClusters,
	model,
	myDNAStringSet,
	N,
	scaleTree,
	.rates,
	defaults,
	weights,
	factr,
	processors=1) {
	defaults <- defaults[1:12]
	if (factr < 1) # below machine precision
		factr <- 1
	indels <- grepl("\\+Indels?", model, ignore.case=TRUE)
	if (indels)
		model <- sub("\\+Indels?", "", model, ignore.case=TRUE)
	empirical <- grepl("\\+F", model, ignore.case=TRUE)
	if (empirical)
		model <- sub("\\+F", "", model, ignore.case=TRUE)
	rates <- as.integer(sub("([^+]*)(\\+G(\\d+))?", "\\3", model, ignore.case=TRUE))
	model <- sub("(.*)\\+G\\d+", "\\1", model, ignore.case=TRUE)
	
	top <- c(0.49, 0.49, 0.49, NA_real_, 1, 1e2, 1e2, 1e2, 1e2, 1e2, 1e2, 500)
	mid <- c(0.25, 0.25, 0.25, 0.25, 0.1, 1, 1, 1, 1, 1, 1, 1)
	low <- c(0.01, 0.01, 0.01, NA_real_, 1e-5, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 0.01)
	var <- mid
	var[4] <- NA_real_
	
	if (model == "JC69") {
		var[c(1:3, 6:10)] <- NA_real_
	} else if (model == "K80") {
		var[c(1:3, 7:10)] <- NA_real_
		mid[6] <- defaults[6]
	} else if (model == "T92") {
		var[c(2:3, 7:10)] <- NA_real_
		mid[c(1, 6)] <- defaults[c(1, 6)]
	} else if (model == "F81") {
		var[6:10] <- NA_real_
		mid[1:4] <- defaults[1:4]
	} else if (model == "HKY85") {
		var[7:10] <- NA_real_
		mid[c(1:4, 6)] <- defaults[c(1:4, 6)]
	} else if (model == "TN93") {
		var[8:10] <- NA_real_
		mid[c(1:4, 6:7)] <- defaults[c(1:4, 6:7)]
	} else if (model == "SYM") {
		var[1:3] <- NA_real_
		mid[6:10] <- defaults[6:10]
	} else { # model == "GTR"
		mid[c(1:4, 6:10)] <- defaults[c(1:4, 6:10)]
	}
	if (indels) {
		mid[c(5, 11)] <- defaults[c(5, 11)]
		if (empirical)
			var[5] <- NA_real_
	} else {
		var[c(5, 11)] <- NA_real_
		mid[c(5, 11)] <- 0
	}
	if (is.na(rates))
		var[12] <- NA_real_
	
	if (empirical) {
		if (model == "T92") {
			defaults[c(1, 4)] <- mean(defaults[c(1, 4)], na.rm=TRUE)
			defaults[2:3] <- (1 - 2*defaults[1])/2
		}
		w <- which(!is.na(var[1:5]))
		if (length(w) > 0) {
			mid[w] <- defaults[w]
			var[w] <- NA_real_
		}
	}
	free <- !is.na(var)
	tot <- sum(free)
	
	f <- function(params) {
		count <<- count + 1L
		if (scaleTree) {
			myClusters[, 4:5] <- params[1]*myClusters[, 4:5]
			params <- params[-1]
		}
		
		vals <- mid
		vals[free] <- params
		if (empirical) {
			vals[4] <- 1 - sum(vals[1:3])
		} else if (free[1]) {
			if (free[2]) { # all bases free
				vals[4] <- 1 - sum(vals[1:3])
			} else {
				vals[1:4] <- c(vals[1],
					rep((1 - 2*vals[1])/2, 2),
					vals[1])
			}
			if (vals[4] <= 0)
				return(highest)
		}
		if (free[5])
			vals[1:5] <- vals[1:5]/sum(vals[1:5])
		if (!free[7])
			vals[7] <- vals[6]
		if (free[12]) {
			vals <- c(vals[-12],
				.rates(vals[12], rates))
		} else {
			vals <- c(vals[-12], 1, 1)
		}
		
		LnL <- .Call("clusterML",
			myClusters,
			myDNAStringSet,
			vals,
			integer(),
			numeric(),
			0,
			1L,
			weights,
			processors,
			PACKAGE="DECIPHER")
		
		if (LnL > highest)
			highest <<- LnL
		
		LnL
	}
	
	highest <- -Inf
	count <- 0L
	if (scaleTree) {
		o <- optim(c(1, mid[free]),
			f,
			method="L-BFGS-B",
			lower=c(0.5, low[free]),
			upper=c(2, top[free]),
			control=list(factr=factr,
				parscale=c(0.1, mid[free]*0.1)))
		var[free] <- o$par[-1]
	} else {
		o <- optim(mid[free],
			f,
			method="L-BFGS-B",
			lower=low[free],
			upper=top[free],
			control=list(factr=factr,
				parscale=mid[free]*0.1))
		var[free] <- o$par
	}
	
	if (empirical && length(w) > 0)
		var[w] <- defaults[w]
	
	LnL <- o$value
	K <- 2*dim(myClusters)[1] - 1 + tot
	AICc <- 2*K + 2*LnL + 2*K*(K + 1)/(N - K - 1)
	BIC <- 2*LnL + K*log(N)
	
	if (scaleTree) {
		c(o$par[1], count, var, LnL, AICc, BIC)
	} else {
		c(1, count, var, LnL, AICc, BIC)
	}
}

.optimizeModelAA <- function(myClusters,
	model,
	myAAStringSet,
	N,
	scaleTree,
	.rates,
	defaults,
	weights,
	factr,
	processors=1) {
	defaults <- defaults[1:213]
	if (factr < 1) # below machine precision
		factr <- 1
	indels <- grepl("\\+Indels?", model, ignore.case=TRUE)
	if (indels)
		model <- sub("\\+Indels?", "", model, ignore.case=TRUE)
	empirical <- grepl("\\+F", model, ignore.case=TRUE)
	if (empirical)
		model <- sub("\\+F", "", model, ignore.case=TRUE)
	rates <- as.integer(sub("([^+]*)(\\+G(\\d+))?", "\\3", model, ignore.case=TRUE))
	model <- sub("(.*)\\+G\\d+", "\\1", model, ignore.case=TRUE)
	
	top <- c(defaults[1:210], 1, 1e2, 500)
	mid <- c(defaults[1:210], 0.1, 1, 1)
	low <- c(defaults[1:210], 1e-5, 1e-2, 0.01)
	var <- mid
	
	if (indels) {
		if (defaults[211] < low[211]) {
			mid[211] <- low[211]
		} else {
			mid[211] <- defaults[211]
		}
		mid[212] <- defaults[212]
		if (empirical)
			var[211] <- NA_real_
	} else {
		var[211:212] <- NA_real_
		mid[211:212] <- 0
	}
	if (is.na(rates))
		var[213] <- NA_real_
	
	free <- c(rep(FALSE, 210), !is.na(var[211:213]))
	tot <- sum(free)
	
	f <- function(params) {
		count <<- count + 1L
		if (scaleTree) {
			myClusters[, 4:5] <- params[1]*myClusters[, 4:5]
			params <- params[-1]
		}
		
		vals <- mid
		vals[free] <- params
		if (free[211])
			vals[191:211] <- vals[191:211]/sum(vals[191:211])
		if (free[213]) {
			vals <- c(vals[-213],
				.rates(vals[213], rates))
		} else {
			vals <- c(vals[-213], 1, 1)
		}
		
		LnL <- .Call("clusterML",
			myClusters,
			myAAStringSet,
			vals,
			integer(),
			numeric(),
			0,
			3L,
			weights,
			processors,
			PACKAGE="DECIPHER")
		
		if (LnL > highest)
			highest <<- LnL
		
		LnL
	}
	
	highest <- -Inf
	count <- 0L
	if (scaleTree) {
		o <- optim(c(1, mid[free]),
			f,
			method="L-BFGS-B",
			lower=c(0.5, low[free]),
			upper=c(2, top[free]),
			control=list(factr=factr,
				parscale=c(0.1, mid[free]*0.1)))
		var[free] <- o$par[-1]
	} else {
		o <- optim(mid[free],
			f,
			method="L-BFGS-B",
			lower=low[free],
			upper=top[free],
			control=list(factr=factr,
				parscale=mid[free]*0.1))
		var[free] <- o$par
	}
	
	if (empirical)
		var[211] <- defaults[211]
	
	LnL <- o$value
	K <- 2*dim(myClusters)[1] - 1 + tot
	AICc <- 2*K + 2*LnL + 2*K*(K + 1)/(N - K - 1)
	BIC <- 2*LnL + K*log(N)
	
	if (scaleTree) {
		c(o$par[1], count, var, LnL, AICc, BIC)
	} else {
		c(1, count, var, LnL, AICc, BIC)
	}
}

.globalBranches <- function(f, # function to minimize
	x, # initial guesses
	h=0.5, # starting neighborhood around center point
	epsilon=1e-4, # convergence precision
	tol=0.001, # accuracy of x values
	coef=0.7, # damping on x to prevent oscillations (0, 1]
	minIterations=5,
	maxIterations=10,
	lower=1e-8) {
	# optimize global branch lengths on whole tree
	x[x < lower] <- lower
	x <- X <- log(x) # enforce positive lengths during optimization
	log_lower <- log(lower)
	prev <- best <- Inf
	W <- seq_along(x)
	prev_delta <- numeric(length(x))
	count <- 0L
	while (length(W) > 0) {
		count <- count + 1L
		H <- h/count
		
		v <- x[W]
		u <- matrix(exp(c(v + H, v - H)), nrow=2L, byrow=TRUE)
		fu <- f(exp(x), rep(W, each=2), u)
		
		f_right <- fu[-1L][c(TRUE, FALSE)]
		f_left <- fu[-1L][c(FALSE, TRUE)]
		fu <- fu[1L]
		change <- prev - fu
		prev <- fu
		
		# only record likeihood improvements
		if (fu < best) {
			X <- x
			best <- prev
		}
		
		# calculate the gradient
		df <- (f_right - f_left)/(2*H)
		ddf <- (f_right - 2*fu + f_left)/H^2
		
		delta <- ifelse(ddf == 0, H, abs(df/ddf))
		delta[delta > 1] <- 1
		delta <- ifelse(df <= 0, delta, -delta)
		# dampen dependency among x with exponential smoothing
		delta <- coef*delta + (1 - coef)*prev_delta[W]
		prev_delta[W] <- delta
		x[W] <- x[W] + delta
		x[W] <- ifelse(x[W] < log_lower, log_lower, x[W])
		
		# `change` is positive when converging
		if (count < minIterations) {
			cont <- df <= -epsilon | df >= epsilon
			W <- W[cont | abs(change) > tol]
		} else if (count < maxIterations) {
			W <- W[abs(change) > tol] # prevent oscillations
		} else { # must continue converging
			W <- W[change > tol] # prevent oscillations
		}
	}
	
	exp(X)
}

.reorderClusters <- function(myClusters, all=FALSE) {
	# order clusters by branching pattern
	repeat {
		a <- pmax(myClusters[, 7], myClusters[, 8])
		w <- which(a > seq_len(nrow(myClusters)))
		if (length(w) == 0)
			break
		w <- w[1]
		temp <- myClusters[w, c(4, 5, 7, 8)]
		myClusters[w:(a[w] - 1), c(4, 5, 7, 8)] <- myClusters[(w + 1):a[w], c(4, 5, 7, 8)]
		myClusters[a[w], c(4, 5, 7, 8)] <- temp
		w1 <- which(myClusters[w:dim(myClusters)[1], 7:8] %in% (w + 1):a[w])
		w2 <- which(myClusters[w:dim(myClusters)[1], 7:8] %in% w)
		if (length(w1) > 0)
			myClusters[w:dim(myClusters)[1], 7:8][w1] <- myClusters[w:dim(myClusters)[1], 7:8][w1] - 1
		if (length(w2) > 0)
			myClusters[w:dim(myClusters)[1], 7:8][w2] <- a[w]
	}
	
	if (all) { # also renumber columns 1 to 3
		count <- 0L
		myClusters[, 1:2] <- myClusters[, 7:8]
		for (i in 1:dim(myClusters)[1]) {
			if (myClusters[i, 7] > 0 && myClusters[i, 8] > 0) {
				myClusters[i, 1] <- myClusters[myClusters[i, 7], 3]
				myClusters[i, 2] <- myClusters[myClusters[i, 8], 3]
				count <- count + 1L
				myClusters[i, 3] <- count
			} else if (myClusters[i, 7] > 0) {
				myClusters[i, 1] <- myClusters[myClusters[i, 7], 3]
				myClusters[i, 3] <- myClusters[i, 1]
			} else if (myClusters[i, 8] > 0) {
				myClusters[i, 2] <- myClusters[myClusters[i, 8], 3]
				myClusters[i, 3] <- myClusters[i, 2]
			} else {
				count <- count + 1L
				myClusters[i, 3] <- count
			}
		}
	}
	
	return(myClusters)
}

.swapBranches <- function(myClusters, r1, c1, r2, c2) {
	# swap branch [r1, c1] with [r2, c2]
	temp <- myClusters[r1, c1]
	myClusters[r1, c1] <- myClusters[r2, c2]
	myClusters[r2, c2] <- temp
	temp <- myClusters[r1, c1 - 3]
	myClusters[r1, c1 - 3] <- myClusters[r2, c2 - 3]
	myClusters[r2, c2 - 3] <- temp
	
	myClusters <- .reorderClusters(myClusters)
	
	return(myClusters)
}

MODELS <- list(Nucleotide=c("JC69",
		"JC69+G4",
		"K80",
		"K80+G4",
		"F81",
		"F81+G4",
		"HKY85",
		"HKY85+G4",
		"T92",
		"T92+G4",
		"TN93",
		"TN93+G4",
		"SYM",
		"SYM+G4",
		"GTR",
		"GTR+G4"),
	Protein=c("AB",
		"AB+G4",
		"BLOSUM62",
		"BLOSUM62+G4",
		"cpREV",
		"cpREV+G4",
		"cpREV64",
		"cpREV64+G4",
		"Dayhoff",
		"Dayhoff+G4",
		"DCMut-Dayhoff",
		"DCMut-Dayhoff+G4",
		"DCMut-JTT",
		"DCMut-JTT+G4",
		"DEN",
		"DEN+G4",
		"FLAVI",
		"FLAVI+G4",
		"FLU",
		"FLU+G4",
		"gcpREV",
		"gcpREV+G4",
		"HIVb",
		"HIVb+G4",
		"HIVw",
		"HIVw+G4",
		"JTT",
		"JTT+G4",
		"LG",
		"LG+G4",
		"MtArt",
		"MtArt+G4",
		"mtDeu",
		"mtDeu+G4",
		"mtInv",
		"mtInv+G4",
		"mtMam",
		"mtMam+G4",
		"mtMet",
		"mtMet+G4",
		"mtOrt",
		"mtOrt+G4",
		"mtREV",
		"mtREV+G4",
		"mtVer",
		"mtVer+G4",
		"MtZoa",
		"MtZoa+G4",
		"PMB",
		"PMB+G4",
		"Q.bird",
		"Q.bird+G4",
		"Q.insect",
		"Q.insect+G4",
		"Q.LG",
		"Q.LG+G4",
		"Q.mammal",
		"Q.mammal+G4",
		"Q.pfam",
		"Q.pfam+G4",
		"Q.plant",
		"Q.plant+G4",
		"Q.yeast",
		"Q.yeast+G4",
		"rtREV",
		"rtREV+G4",
		"stmtREV",
		"stmtREV+G4",
		"VT",
		"VT+G4",
		"WAG",
		"WAG+G4",
		"WAGstar",
		"WAGstar+G4"))

.splitClusters <- function(x, y) {
	clusterNum <- 0L
	X <- integer(length(x))
	u.y <- unique(y)
	for (i in u.y) {
		w.y <- which(y==i)
		u.x <- unique(x[w.y])
		for (j in u.x) {
			clusterNum <- clusterNum + 1L
			w.x <- which(x[w.y]==j)
			X[w.y[w.x]] <- clusterNum
		}
	}
	return(X)
}

.root <- function(x1, root) {
	# if root is zero then midpoint root the tree
	# otherwise outgroup root based on root index
	# (note: output columns 1:3 are uncorrected)
	
	n <- nrow(x1)
	if (root==0) { # midpoint root
		# find the leaf at minimum height
		r1 <- which(x1[, 7] < 0)
		h1 <- x1[r1, 6] - x1[r1, 4]
		r2 <- which(x1[, 8] < 0)
		h2 <- x1[r2, 6] - x1[r2, 5]
		r <- c(r1, r2) # row number of leaf
		z <- rep(c(7L, 8L), c(length(r1), length(r2)))
		h <- c(h1, h2) # height of leaf
		
		# reorder by sequence number
		o <- order(x1[cbind(r, z)],
			decreasing=TRUE)
		h <- h[o]
		minH <- which.min(h) # index of lowest leaf
	} else { # outgroup root
		w <- which(x1[n, 7:8]==-root)
		if (length(w) > 0) { # already outgroup rooted
			# extend the root node
			x1[n, 6] <- x1[n, 6] + x1[n, 3 + w]
			x1[n, 6 - w] <- x1[n, 6 - w] + x1[n, 3 + w]
			x1[n, 3 + w] <- 0
			return(x1)
		}
		minH <- root # index of root
	}
	
	# find most distant leaf from minH
	longest <- numeric(n) # length of longest path
	merged <- logical(n) # whether merged yet
	index <- numeric(n) # column back to minH
	for (i in seq_len(n)) {
		b1 <- x1[i, 7]
		if (b1 < 0) { # merged with leaf
			if (b1==-minH) {
				merged[i] <- TRUE
				b1 <- NA_real_
			} else {
				l1 <- x1[i, 4]
			}
		} else { # merged with node
			if (merged[b1]) {
				merged[i] <- TRUE
				b1 <- NA_real_
			} else {
				l1 <- longest[b1] + x1[i, 4]
			}
		}
		
		b2 <- x1[i, 8]
		if (b2 < 0) { # merged with leaf
			if (b2==-minH) {
				merged[i] <- TRUE
				b2 <- NA_real_
			} else {
				l2 <- x1[i, 5]
			}
		} else { # merged with node
			if (merged[b2]) {
				merged[i] <- TRUE
				b2 <- NA_real_
			} else {
				l2 <- longest[b2] + x1[i, 5]
			}
		}
		
		if (is.na(b1)) { # b1 contains minH
			longest[i] <- l2
			# leave index[i] at zero
		} else if (is.na(b2)) { # b2 contains minH
			longest[i] <- l1
			index[i] <- 1
		} else if (l1 >= l2) {
			longest[i] <- l1
			# index[i] not needed
		} else { # l2 > l1
			longest[i] <- l2
			# index[i] not needed
		}
	}
	
	if (root==0) { # determine height of the midpoint
		w <- which(merged)
		longest <- longest + x1[, 6] - h[minH]
		m <- w[which.max(longest[w])]
		midH <- longest[m]/2
		if (isTRUE(all.equal(x1[n, 6], midH)))
			return(x1) # already midpoint rooted
		
		# find the edge containing the midpoint
		lowH <- x1[m, 6] - x1[m, 4 + index[m]]
		while (lowH > midH) { # descend the tree
			m <- x1[m, 7 + index[m]]
			lowH <- x1[m, 6] - x1[m, 4 + index[m]]
		}
	} else { # root at tip of outgroup
		w <- which(x1[, 7:8]==-root, arr.ind=TRUE)
		midH <- x1[w[1], 6] - x1[w[1], 3 + w[2]]
		m <- w[1]
	}
	
	# invert and lower nodes above rotation point
	.dropH <- function(i, delta) {
		stack <- integer(n)
		pos <- 1L
		stack[pos] <- i
		while (pos > 0) {
			i <- stack[pos]
			x1[i, 6] <<- x1[i, 6] - delta
			pos <- pos - 1L
			if (x1[i, 7] > 0) {
				pos <- pos + 1L
				stack[pos] <- x1[i, 7]
			}
			if (x1[i, 8] > 0) {
				pos <- pos + 1L
				stack[pos] <- x1[i, 8]
			}
		}
	}
	up <- integer(n) # pointers up tree
	w <- which(x1[, 7:8] > 0, arr.ind=TRUE)
	if (length(w) > 0)
		up[x1[, 7:8][w]] <- w[, "row"]
	remove <- logical(n) # replaced nodes
	x2 <- x1 # new rooted tree
	count <- n # row in x2
	# make new root node
	delta <- x1[m, 6] - midH
	x2[count, 4:10] <- c(x1[m, 4 + index[m]] - delta,
		delta,
		midH,
		x1[m, 7 + index[m]],
		count - 1,
		x1[m, 9 + index[m]],
		-1)
	if (up[m]) {
		while (up[m]) {
			count <- count - 1
			delta <- x1[m, 6] - midH
			remove[m] <- TRUE
			x2[count, 4:10] <- c(x1[m, 5 - index[m]],
				x1[up[m], 4 + index[up[m]]],
				midH - delta,
				x1[m, 8 - index[m]],
				count - 1,
				x1[m, 10 - index[m]],
				-1)
			if (x1[m, 8 - index[m]] > 0)
				.dropH(x1[m, 8 - index[m]], 2*delta)
			m <- up[m]
		}
		delta <- x1[m, 6] - midH
		x2[count, 5] <- sum(x1[m, 4:5])
	}
	remove[m] <- TRUE
	keep <- which(!remove)
	x2[count, 8] <- x1[m, 8 - index[m]]
	if (x2[count, 8] > 0)
		x2[count, 8] <- match(x2[count, 8], keep)
	x2[count, 10] <- x1[m, 10 - index[m]]
	if (x1[m, 8 - index[m]] > 0)
		.dropH(x1[m, 8 - index[m]], 2*delta)
	if (length(keep) > 0) {
		x2[1:(count - 1),] <- x1[keep,]
		w <- which(x2[1:(count - 1), 7:8] > 0)
		if (length(w) > 0)
			x2[1:(count - 1), 7:8][w] <- match(x2[1:(count - 1), 7:8][w], keep)
		w <- which(x2[n:count, 7] %in% keep)
		x2[n:count, 7][w] <- match(x2[n:count, 7][w], keep)
	}
	
#	w <- which(x2[, 7] > 0)
#	if (length(w) > 0)
#		x2[w, 4] <- x2[w, 6] - x2[x2[w, 7], 6]
#	w <- which(x2[, 8] > 0)
#	if (length(w) > 0)
#		x2[w, 5] <- x2[w, 6] - x2[x2[w, 8], 6]
	
	if (root > 0)
		x2[, 6] <- x2[, 6] - min(x2[, 6] - x2[, 4], x2[, 6] - x2[, 5])
	
	return(x2)
}

.applyMidpoints <- function(dend, dim) {
	# initialize a stack of maximum length (dim)
	stack <- vector("list", dim)
	visit <- logical(dim) # node already visited
	parent <- integer(dim) # index of parent node
	index <- integer(dim) # index in parent node
	pos <- 1L # current position in the stack
	stack[[pos]] <- dend
	while (pos > 0L) { # more nodes to visit
		if (visit[pos]) { # ascending tree
			visit[pos] <- FALSE # reset visit
			
			members <- sapply(stack[[pos]],
				function(x) {
					m <- attr(x, "members")
					if (is.null(m)) {
						return(1L)
					} else {
						return(m)
					}
				})
			
			l <- length(stack[[pos]])
			if (is.leaf(stack[[pos]][[1]]) && is.leaf(stack[[pos]][[l]])) {
				attr(stack[[pos]], "midpoint") <- (sum(members) - 1)/2
			} else if (is.leaf(stack[[pos]][[1]])) {
				attr(stack[[pos]], "midpoint") <- (sum(members[-l]) + attr(stack[[pos]][[l]], "midpoint"))/2
			} else if (is.leaf(stack[[pos]][[l]])) {
				attr(stack[[pos]], "midpoint") <- (attr(stack[[pos]][[1]], "midpoint") + sum(members[-l]))/2
			} else {
				attr(stack[[pos]], "midpoint") <- (sum(members[-l]) + attr(stack[[pos]][[1]], "midpoint") + attr(stack[[pos]][[l]], "midpoint"))/2
			}
			
			# replace self in parent
			if (parent[pos] > 0)
				stack[[parent[pos]]][[index[pos]]] <- stack[[pos]]
			pos <- pos - 1L # pop off of stack
		} else { # descending tree
			visit[pos] <- TRUE
			p <- pos
			for (i in seq_along(stack[[p]])) {
				if (!is.leaf(stack[[p]][[i]])) {
					# push subtree onto stack
					pos <- pos + 1L
					stack[[pos]] <- stack[[p]][[i]]
					parent[[pos]] <- p
					index[[pos]] <- i
				}
			}
		}
	}
	return(stack[[1L]])
}

.adjustTreeHeights <- function(myClusters) {
	# given myClusters return adjusted heights
	cumHeight <- numeric(max(myClusters[, 3]))
	for (i in 1:dim(myClusters)[1]) {
		if (myClusters[i, 1] < 0 && myClusters[i, 2] < 0) {
			cumHeight[myClusters[i, 3]] <- max(myClusters[i, 4], myClusters[i, 5])
			myClusters[i, 6] <- cumHeight[myClusters[i, 3]]
		} else if (myClusters[i, 1] > 0 && myClusters[i, 2] > 0) {
			cumHeight[myClusters[i, 3]] <- max(myClusters[i, 4] + cumHeight[myClusters[i, 1]],
				myClusters[i, 5] + cumHeight[myClusters[i, 2]])
			myClusters[i, 6] <- cumHeight[myClusters[i, 3]]
		} else if (myClusters[i, 1] > 0) {
			cumHeight[myClusters[i, 3]] <- cumHeight[myClusters[i, 1]] + myClusters[i, 4]
			if (myClusters[i, 5] > cumHeight[myClusters[i, 3]])
				cumHeight[myClusters[i, 3]] <- myClusters[i, 5]
			myClusters[i, 6] <- cumHeight[myClusters[i, 3]]
		} else {
			cumHeight[myClusters[i, 3]] <- cumHeight[myClusters[i, 2]] + myClusters[i, 5]
			if (myClusters[i, 4] > cumHeight[myClusters[i, 3]])
				cumHeight[myClusters[i, 3]] <- myClusters[i, 4]
			myClusters[i, 6] <- cumHeight[myClusters[i, 3]]
		}
	}
	
	myClusters <- .Call("adjustHeights",
		myClusters,
		PACKAGE="DECIPHER")
	return(myClusters)
}

.getParams <- function(myClusters) {
	j <- which(myClusters[, 4:5] < 0)
	myClusters[, 4:5][j] <- -1*myClusters[, 4:5][j]
	lower <- params <- upper <- as.vector(myClusters[, 4:5])
	if (length(j) > 0) {
		j <- c(j,
			ifelse(j > nrow(myClusters),
				j - nrow(myClusters),
				j + nrow(myClusters)))
		j <- c(j, which(myClusters[, 7:8] %in% j))
		lower[j] <- 0
		upper[j] <- 3*max(params[j]) + 0.01
	}
	list(myClusters, lower, params, upper)
}

.giveParamsDNA <- function(model_params, model, .rates) {
	model_params <- as.numeric(model_params[1:12])
	if (!is.na(model_params[1])) {
		if (!is.na(model_params[2])) { # all bases free
			model_params[4] <- 1 - sum(model_params[1:3])
		} else {
			model_params[1:4] <- c(model_params[1],
				rep((1 - 2*model_params[1])/2, 2),
					model_params[1])
		}
	}
	if (is.na(model_params[7]))
		model_params[7] <- model_params[6]
	w <- which(is.na(model_params))
	if (length(w) > 0)
		model_params[w] <- c(0.25, 0.25, 0.25, 0.25, 0, 1, 1, 1, 1, 1, 0, NA)[w]
	
	model_params[1:5] <- model_params[1:5]/sum(model_params[1:5])
	if (is.na(model_params[12])) {
		model_params <- c(model_params[1:11], 1, 1)
	} else {
		model_params <- c(model_params[1:11],
			.rates(model_params[12],
				as.integer(sub(".+(\\+G(\\d+)).*",
					"\\2",
					model))))
	}
}

.giveParamsAA <- function(model_params, model, .rates) {
	model_params <- as.numeric(model_params[1:213])
	if (is.na(model_params[211]))
		model_params[211:212] <- 0
	model_params[191:211] <- model_params[191:211]/sum(model_params[191:211])
	if (is.na(model_params[213])) {
		model_params <- c(model_params[1:212], 1, 1)
	} else {
		model_params <- c(model_params[1:212],
			.rates(model_params[213],
				as.integer(sub(".+(\\+G(\\d+)).*",
					"\\2",
					model))))
	}
}

.colEdge <- function(dend, dim, labels, r, c) {
	# color edges by cluster
	
	# initialize a stack of maximum length (dim)
	stack <- vector("list", dim)
	visit <- logical(dim) # node already visited
	parent <- integer(dim) # index of parent node
	index <- integer(dim) # index in parent node
	pos <- 1L # current position in the stack
	stack[[pos]] <- dend
	while (pos > 0L) { # more nodes to visit
		if (visit[pos]) { # ascending tree
			visit[pos] <- FALSE # reset visit
			
			for (i in seq_along(stack[[pos]])) {
				if (is.leaf(stack[[pos]][[i]])) {
					a <- attributes(stack[[pos]][[i]])
					num <- c$cluster[which(labels==as.character(a$label))]
					attr(stack[[pos]][[i]], "edgePar") <- list(col=r[num %% 15 + 1])
				}
			}
			
			# replace self in parent
			if (parent[pos] > 0)
				stack[[parent[pos]]][[index[pos]]] <- stack[[pos]]
			pos <- pos - 1L # pop off of stack
		} else { # descending tree
			visit[pos] <- TRUE
			p <- pos
			for (i in seq_along(stack[[p]])) {
				if (!is.leaf(stack[[p]][[i]])) {
					# push subtree onto stack
					pos <- pos + 1L
					stack[[pos]] <- stack[[p]][[i]]
					parent[[pos]] <- p
					index[[pos]] <- i
				}
			}
		}
	}
	return(stack[[1L]])
}

.reorder <- function(dend, dim, c) {
	# initialize a stack of maximum length (dim)
	stack <- vector("list", dim)
	visit <- logical(dim) # node already visited
	parent <- integer(dim) # index of parent node
	index <- integer(dim) # index in parent node
	pos <- 1L # current position in the stack
	stack[[pos]] <- dend
	while (pos > 0L) { # more nodes to visit
		if (visit[pos]) { # ascending tree
			visit[pos] <- FALSE # reset visit
			
			members <- lapply(stack[[pos]], unlist)
			# sort tree by ascending cluster number
			o <- sort.list(sapply(members,
					function(x)
						min(c[x, 1])))
			stack[[pos]][] <- stack[[pos]][o]
			
			# replace self in parent
			if (parent[pos] > 0)
				stack[[parent[pos]]][[index[pos]]] <- stack[[pos]]
			pos <- pos - 1L # pop off of stack
		} else { # descending tree
			visit[pos] <- TRUE
			p <- pos
			for (i in seq_along(stack[[p]])) {
				if (!is.leaf(stack[[p]][[i]])) {
					# push subtree onto stack
					pos <- pos + 1L
					stack[[pos]] <- stack[[p]][[i]]
					parent[[pos]] <- p
					index[[pos]] <- i
				}
			}
		}
	}
	return(stack[[1L]])
}

.mask <- function(myXStringSet) {
	v <- TerminalChar(myXStringSet)
	v <- mapply(function(i, j, l) {
			if (i == j) # possibly all gaps
				j <- 0
			IRanges(c(seq_len(i),
					l - seq_len(j) + 1),
				width=1)
		},
		v[, 1],
		v[, 2],
		width(myXStringSet)[1])
	myXStringSet <- replaceAt(myXStringSet,
		IRangesList(unname(v)),
		"+")
	v <- vmatchPattern("-", myXStringSet)
	keep <- lapply(v,
		function(x)
			start(x)[c(1L, which(diff(start(x)) > 1) + 1L)])
	keep <- unique(unlist(keep))
	v <- lapply(v,
		function(x)
			x[!(start(x) %in% keep)])
	myXStringSet <- replaceAt(myXStringSet,
		IRangesList(unname(v)),
		"+")
	myXStringSet
}

.getClusters <- function(i, l, myClusters) {
	if (l <= 0)
		return(i)
	l <- l - sum(myClusters[i, 7:8] > 0)
	if (myClusters[i, 7] > 0) {
		k1 <- .getClusters(myClusters[i, 7],
			l - 1L,
			myClusters)
	} else {
		k1 <- integer()
	}
	if (myClusters[i, 8] > 0) {
		k2 <- .getClusters(myClusters[i, 8],
			l - 1L,
			myClusters)
	} else {
		k2 <- integer()
	}
	return(c(i, k1, k2))
}

.cophenetic <- function(myClusters) {
	n <- nrow(myClusters)
	one <- two <- vector("list", n)
	d <- matrix(1e-6, n + 1, n + 1) # initialize positive
	for (i in seq_len(n)) {
		j <- myClusters[i, 7]
		if (j > 0) {
			one[[i]] <- c(one[[j]], two[[j]])
		} else {
			one[[i]] <- -j
		}
		j <- one[[i]]
		d[j, -j] <- d[j, -j] + myClusters[i, 4]
		d[-j, j] <- d[-j, j] + myClusters[i, 4]
		
		j <- myClusters[i, 8]
		if (j > 0) {
			two[[i]] <- c(one[[j]], two[[j]])
		} else {
			two[[i]] <- -j
		}
		j <- two[[i]]
		d[j, -j] <- d[j, -j] + myClusters[i, 5]
		d[-j, j] <- d[-j, j] + myClusters[i, 5]
	}
	as.dist(d)
}

.rowSums <- function(d) {
	n <- attr(d, "Size")
	r <- numeric(n)
	for (k in seq_along(d)) {
		j <- floor((2*n + 1 - sqrt((2*n - 1)^2 - 8*(k - 1)))/2)
		i <- j + k - (2*n - j)*(j - 1)/2
		r[i] <- r[i] + d[k]
		r[j] <- r[j] + d[k]
	}
	r
}

.makeS <- function(Q, c, v=0.1, noise=0, resolution=10) {
	# make substitution matrix
	if (c <= 5) { # nucleotides
		s <- seq_len(c)
		type <- 1L
		if (noise > 0) {
			Q[1:12] <- abs(Q[1:12] + rnorm(12, sd=Q[1:12]*noise))
			Q[s] <- Q[s]/sum(Q[s]) # normalize state frequencies
		}
	} else { # amino acids
		s <- seq_len(c) + 190
		type <- 3L
		if (noise > 0) {
			Q[1:212] <- abs(Q[1:212] + rnorm(212, sd=Q[1:212]*noise))
			Q[s] <- Q[s]/sum(Q[s]) # normalize state frequencies
		}
	}
	
	# log((p(x)*p(x->y|v))/(p(x)*p(y)))
	P <- .Call("expM", v, Q, type) # convert rate matrix to probabilities
	S <- log(P[seq_len(c), seq_len(c)]/Q[s]) # substitution matrix
	diag(S) <- 0 # no score for no change
	S <- abs(S) # more positive is worse score
	
	# factor into powers of -2 for floating point reproducibility
	powers <- 2^(seq(0, 1 - resolution))
	for (i in seq_len(nrow(S) - 1L)) {
		for (j in (i + 1L):nrow(S)) {
			res <- logical(resolution)
			num <- S[i, j]
			for (k in seq_len(resolution)) {
				res[k] <- floor(num/powers[k])
				num <- num - res[k]*powers[k]
			}
			S[i, j] <- S[j, i] <- sum(res*2^(1 - seq_along(res)))
		}
	}
	
	S
}

.Sankoff <- function(C, z, S, weights, scoreOnly=TRUE, add=0, processors=1L) {
	if (is.double(C))
		mode(C) <- "integer"
	
	c <- nrow(S)
	l <- ncol(z)
	n <- nrow(C)
	
	result <- .Call("clusterMP",
		C,
		z,
		S,
		c(nrow(S), ncol(z), nrow(C), nrow(z)),
		scoreOnly,
		add,
		weights,
		processors)
	
	result
}

.reorder2 <- function(myClusters) {
	repeat {
		a <- pmax(myClusters[, 1], myClusters[, 2])
		w <- which(a > seq_len(nrow(myClusters)))
		if (length(w) == 0)
			break
		w <- w[1]
		temp <- myClusters[w,]
		myClusters[w:(a[w] - 1),] <- myClusters[(w + 1):a[w],]
		myClusters[a[w],] <- temp
		w1 <- which(myClusters[w:dim(myClusters)[1],] %in% (w + 1):a[w])
		w2 <- which(myClusters[w:dim(myClusters)[1],] %in% w)
		if (length(w1) > 0)
			myClusters[w:dim(myClusters)[1],][w1] <- myClusters[w:dim(myClusters)[1],][w1] - 1
		if (length(w2) > 0)
			myClusters[w:dim(myClusters)[1],][w2] <- a[w]
	}
	myClusters
}

.clusterMP <- function(z, # integer encoded XStringSet
	S, # substitution matrix
	seed, # order of addition
	NNIs=10, # how often to attempt NNIs
	weights,
	processors) {
	l <- nrow(z)
	
	C <- matrix(NA_integer_, l - 1, 2)
	C[1,] <- -seed[1:2]
	for (i in 1L:(l - 2L)) {
		if (i %% NNIs == 0) { # attempt NNIs
			res <- .Sankoff(C[seq_len(i),, drop=FALSE],
				z,
				S,
				weights,
				TRUE,
				-1, # add
				processors)
			
			count <- 1L
			Corg <- C
			for (j in rev(seq_len(i))) {
				if (Corg[j, 1] > 0) {
					count <- count + 1L
					if (res[1] > res[count] &&
						C[j, 1] > 0 &&
						all(C[c(j, C[j, 1]),] == Corg[c(j, C[j, 1]),])) {
						# swap left-left with right
						temp <- C[C[j, 1], 1]
						C[C[j, 1], 1] <- C[j, 2]
						C[j, 2] <- temp
						C <- .reorder2(C)
					}
					
					count <- count + 1L
					if (res[1] > res[count] &&
						C[j, 1] > 0 &&
						all(C[c(j, C[j, 1]),] == Corg[c(j, C[j, 1]),])) {
						# swap left-right with right
						temp <- C[C[j, 1], 2]
						C[C[j, 1], 2] <- C[j, 2]
						C[j, 2] <- temp
						C <- .reorder2(C)
					}
				}
				
				if (Corg[j, 2] > 0) {
					count <- count + 1L
					if (res[1] > res[count] &&
						C[j, 2] > 0 &&
						all(C[c(j, C[j, 2]),] == Corg[c(j, C[j, 2]),])) {
						# swap right-left with left
						temp <- C[C[j, 2], 1]
						C[C[j, 2], 1] <- C[j, 1]
						C[j, 1] <- temp
						C <- .reorder2(C)
					}
					
					count <- count + 1L
					if (res[1] > res[count] &&
						C[j, 2] > 0 &&
						all(C[c(j, C[j, 2]),] == Corg[c(j, C[j, 2]),])) {
						# swap right-right with left
						temp <- C[C[j, 2], 2]
						C[C[j, 2], 2] <- C[j, 1]
						C[j, 1] <- temp
						C <- .reorder2(C)
					}
				}
			}
		}
		
		res <- .Sankoff(C[seq_len(i),, drop=FALSE],
			z,
			S,
			weights,
			TRUE,
			seed[i + 2],
			processors)
		
		w <- which.min(res[-1])
		if (w > i) {
			w <- c(w - i, 2L)
		} else {
			w <- c(w, 1L)
		}
		
		# add sequence to the tree
		if (w[1] > i) { # add to root
			C[w[1],] <- c(i, -seed[i + 2])
		} else {
			C[w[1]:i,][C[w[1]:i,] %in% w[1]:i] <- C[w[1]:i,][C[w[1]:i,] %in% w[1]:i] + 1L
			C[(w[1] + 1):(i + 1),] <- C[w[1]:i,]
			C[w[1], 3 - w[2]] <- -seed[i + 2]
			C[w[1] + 1, w[2]] <- w[1]
		}
	}
	
	c(list(C), .Sankoff(C, z, S, weights, FALSE, processors=processors))
}

.guessParams <- function(a, indels) {
	if (indels) {
		b <- a[c(1:4, 16),, drop=FALSE]
	} else {
		b <- a[1:4,, drop=FALSE]
	}
	fr <- rowSums(b)
	fr <- ifelse(fr == 0, 1, fr)
	fr <- fr/sum(fr)
	if (!indels)
		fr <- c(fr, 0)
	fr <- unname(fr)
	
	k0 <- sum(b[3,]*b[4,])
	if (k0 == 0)
		k0 <- 1
	k1 <- sum(b[1,]*b[3,])/k0
	if (k1 == 0)
		k1 <- 1/k0
	k2 <- sum(b[2,]*b[4,])/k0
	if (k2 == 0)
		k2 <- 1/k0
	k3 <- sum(b[1,]*b[2,])/k0
	if (k3 == 0)
		k3 <- 1/k0
	k4 <- sum(b[1,]*b[4,])/k0
	if (k4 == 0)
		k4 <- 1/k0
	k5 <- sum(b[2,]*b[3,])/k0
	if (k5 == 0)
		k5 <- 1/k0
	
	c(fr, k1, k2, k3, k4, k5, 1, 1)
}

.extractClades <- function(myClusters) {
	n <- nrow(myClusters)
	first <- as.integer(myClusters[, 7L])
	second <- as.integer(myClusters[, 8L])
	one <- two <- vector("list", n)
	for (i in seq_len(n)) {
		j <- first[i]
		if (j > 0) {
			m1 <- min(one[[j]])
			m2 <- min(two[[j]])
			if (m1 < m2) {
				one[[i]] <- c(one[[j]], two[[j]])
			} else {
				one[[i]] <- c(two[[j]], one[[j]])
			}
		} else {
			one[[i]] <- -j
		}
		
		j <- second[i]
		if (j > 0) {
			m1 <- min(one[[j]])
			m2 <- min(two[[j]])
			if (m1 < m2) {
				two[[i]] <- c(one[[j]], two[[j]])
			} else {
				two[[i]] <- c(two[[j]], one[[j]])
			}
		} else {
			two[[i]] <- -j
		}
	}
	c(one, two)
}

.localBranches <- function(myClusters,
	myXStringSet,
	model_params,
	weights_ML,
	processors,
	h=0.1, # starting neighborhood around center point
	epsilon=1e-5, # convergence precision
	tol=0.001, # convergence tolerance on likelihood
	absTol=0.1, # only converge when within absTol of LnL or better
	minIterations=5,
	maxIterations=20,
	lower=1e-6) {
	# optimize local branch lengths for all possible NNIs
	n <- nrow(myClusters)
	w <- which(myClusters[, 7:8] > 0, arr.ind=TRUE)
	w <- unname(w)
	
	Up <- integer(n)
	Up[myClusters[, 7:8][w]] <- w[, 1L]
	
	if (any(myClusters[nrow(myClusters), 7:8] < 0))
		w <- w[w[, 1L] != nrow(myClusters),, drop=FALSE] # trifurcation at root
	w <- rbind(cbind(w, 1L), cbind(w, -1L)) # row, col, flip
	
	ind <- matrix(NA_integer_, nrow=5, ncol=nrow(w))
	# fill rows: center, opposite, down-left, down-right, up
	for (i in seq_len(nrow(w))) {
		if (w[i, 2L] == 1L) { # center branch is left
			ind[1L, i] <- w[i, 1L]
			if (w[i, 1L] == nrow(myClusters)) { # root
				ind[2L, i] <- as.integer(myClusters[w[i, 1L], 8]) # opposite-left
			} else {
				ind[2L, i] <- ind[1L, i] + n
			}
		} else { # center branch is right
			ind[1L, i] <- w[i, 1L] + n
			if (w[i, 1L] == nrow(myClusters)) { # root
				ind[2L, i] <- as.integer(myClusters[w[i, 1L], 7]) # opposite-left
			} else {
				ind[2L, i] <- w[i, 1L]
			}
		}
		if (w[i, 3L] > 0) { # swap down-left with right
			if (w[i, 2L] == 1L) {
				ind[3L, i] <- as.integer(myClusters[w[i, 1L], 7])
			} else {
				ind[3L, i] <- as.integer(myClusters[w[i, 1L], 8])
			}
			ind[4L, i] <- ind[3L, i] + n
		} else { # swap down-right with right
			if (w[i, 2L] == 1L) {
				ind[4L, i] <- as.integer(myClusters[w[i, 1L], 7])
			} else {
				ind[4L, i] <- as.integer(myClusters[w[i, 1L], 8])
			}
			ind[3L, i] <- ind[4L, i] + n
		}
		if (w[i, 1L] == nrow(myClusters)) { # root
			ind[5L, i] <- ind[2L, i] + n # opposite-right
		} else {
			ind[5L, i] <- Up[w[i, 1L]]
			if (myClusters[ind[5L, i], 7] != w[i, 1L])
				ind[5L, i] <- ind[5L, i] + n
		}
	}
	
	# initialize branch lengths at mean expected value per quartet
	x <- matrix(myClusters[, 4:5][as.vector(ind)],
		nrow=5,
		ncol=nrow(w))
	coefs1 <- matrix(c(-0.01, 0.22, 0.64, 0.2, 0.61,
			0.07, -0.09, 1.02, -0.13, 0.1,
			0.02, 1.07, -0.12, 0.17, -0.17,
			0.05, 0.14, -0.17, 1.03, -0.11,
			-0.01, -0.13, 0.17, -0.07, 1.06),
		nrow=5,
		ncol=5)
	coefs2 <- matrix(c(-0.04, 0.22, 0.27, 0.61, 0.58,
			0.05, -0.08, -0.14, 1.04, 0.12,
			0.04, 1.05, 0.17, -0.14, -0.19,
			0.05, 0.13, 1.03, -0.15, -0.1,
			0.02, -0.14, -0.09, 0.16, 1.05),
		nrow=5,
		ncol=5)
	W <- w[, 3L] > 0
	x[, W] <- coefs1 %*% x[, W, drop=FALSE]
	W <- w[, 3L] < 0
	x[, W] <- coefs2 %*% x[, W, drop=FALSE]
	# impose reasonable stating conditions
	x[x < 0.001] <- 0.001
	x[x > 1] <- 1
	x <- X <- log(x) # enforce positive lengths during optimization
	
	log_lower <- log(lower)
	W <- seq_len(nrow(w))
	res <- rep(Inf, nrow(w))
	count <- 0L
	while (length(W) > 0) {
		count <- count + 1L
		H <- h/count # adaptively shrink h
		# apply a mask at around point (x +/- h)
		mask <- c(0, 0, 0, 0, 0,
			H, 0, 0, 0, 0,
			0, H, 0, 0, 0,
			0, 0, H, 0, 0,
			0, 0, 0, H, 0,
			0, 0, 0, 0, H,
			-H, 0, 0, 0, 0,
			0, -H, 0, 0, 0,
			0, 0, -H, 0, 0,
			0, 0, 0, -H, 0,
			0, 0, 0, 0, -H)
		reps <- rep(W, each=11)
		lengths <- x[, reps] + mask
		branches <- rep(ind[1L, W]*w[W, 3L], each=11)
		
		LnL <- .Call("clusterML",
			myClusters,
			myXStringSet,
			model_params,
			branches,
			exp(lengths), # prevents negative lengths
			0,
			ifelse(is(myXStringSet, "AAStringSet"), 3L, 1L),
			weights_ML,
			processors,
			PACKAGE="DECIPHER")
		
		f <- LnL[-1L][c(TRUE, rep(FALSE, 10))]
		change <- res[W] - f
		
		# only record likeihood improvements
		keep <- f <= res[W] # minimum observed
		res[W[keep]] <- f[keep]
		X[, W[keep]] <- x[, W[keep]]
		
		f <- rep(f, each=5)
		f_left <- LnL[-1L][c(rep(FALSE, 6), rep(TRUE, 5))]
		f_right <- LnL[-1L][c(FALSE, rep(TRUE, 5), rep(FALSE, 5))]
		
		LnL <- LnL[1L]
		
		# calculate the gradient
		df <- (f_right - f_left)/(2*H)
		ddf <- (f_right - 2*f + f_left)/H^2
		
		delta <- ifelse(ddf == 0, H, abs(df/ddf))
		delta[delta > 1] <- 1
		delta <- ifelse(df <= 0, delta, -delta)
		x[, W] <- x[, W] + delta
		x[, W] <- ifelse(x[, W] < log_lower, log_lower, x[, W])
		
		# `change` is positive when converging
		if (count < minIterations) {
			conv <- matrix(df > -epsilon & df < epsilon,
				ncol=length(W))
			W <- W[colSums(conv) != 5 | # early convergence
				abs(change) > tol]
		} else if (count < maxIterations) {
			# project best foreseeable improvement in likelihood
			projected <- df*delta + (ddf*delta^2)/2
			projected <- f + (maxIterations - count)*projected
			conv <- matrix(projected > LnL + absTol, # extrapolation
				ncol=length(W))
			W <- W[colSums(conv) != 5 & # unlikely to surpass LnL
				abs(change) > tol] # normal convergence
		} else { # must continue converging
			# project likelihood in next iteration
			projected <- f + df*delta + (ddf*delta^2)/2
			# only interested in cases where better than LnL
			conv <- matrix(projected > LnL,
				ncol=length(W))
			W <- W[colSums(conv) != 5 & # improvement projected
				change > tol] # prevent oscillations
		}
	}
	
	list(w, res, exp(X), LnL)
}

.rNNIs <- function(myClusters, fracRandomNNIs, prob) {
	if (fracRandomNNIs == 0)
		return(myClusters)
	
	w <- which(myClusters[, 7:8, drop=FALSE] > 0, arr.ind=TRUE)
	s <- sample(seq_along(prob),
		ceiling(length(prob)*fracRandomNNIs),
		prob=prob)
	s <- s[order(w[s, 1L])]
	flip <- sample(0:1, length(s), replace=TRUE)
	for (j in seq_along(s)) {
		i <- w[s[j], 1L]
		side <- w[s[j], 2L]
		if (i == nrow(myClusters)) { # root
			if (myClusters[i, 7] > 0 && myClusters[i, 8] > 0) {
				if (flip[j] > 0) {
					myClusters <- .swapBranches(myClusters,
						myClusters[i, 7], 8,
						myClusters[i, 8], 7)
				} else {
					myClusters <- .swapBranches(myClusters,
						myClusters[i, 7], 8,
						myClusters[i, 8], 8)
				}
			}
		} else if (side == 1L) {
			if (myClusters[i, 7] > 0) {
				if (flip[j] > 0) {
					myClusters <- .swapBranches(myClusters,
						myClusters[i, 7], 7,
						i, 8)
				} else {
					myClusters <- .swapBranches(myClusters,
						myClusters[i, 7], 8,
						i, 8)
				}
			}
		} else if (myClusters[i, 8] > 0) {
			if (flip[j] > 0) {
				myClusters <- .swapBranches(myClusters,
					myClusters[i, 8], 7,
					i, 7)
			} else {
				myClusters <- .swapBranches(myClusters,
					myClusters[i, 8], 8,
					i, 7)
			}
		}
	}
	
	myClusters
}

.redundancy <- function(y, i) {
	# downweight redundancy in y
	c <- 1 - cor(t(y))^2 # 1 - R^2
	count <- 1L
	w <- numeric(nrow(c))
	w[i] <- 1L
	while(count < nrow(c)) {
		count <- count + 1L
		c[, i] <- NA_real_
		j <- which.min(c[i,]) # highest correlation
		w[j] <- c[i, j]
		i <- j
	}
	w
}

.estimate <- function(x, # dependent variable
	y, # independent variables
	base, # base of weight function
	tol=0.9, # minimum cumulative stdev of principal components
	interval=0, # exclusion interval for normal quantile [0, 1)
	#eps=0.1, # standard error of random projection
	mult, # bounds of multiplier on noise (>= 0)
	N=1) {
	#j <- nrow(y)
	#k <- ncol(y)
	#orderPCA <- j*k^2 + k^3
	#n <- ceiling(9*log(j)/(eps^2 - 2*eps^3/3) + 1)
	#orderRP <- (j + 2)*k*n + n*j^2
	
	# perform dimensionality-reduction in log space
	y <- log(y) # error ~ lognormal
	
	#if (orderPCA < orderRP) {
		PCA <- prcomp(y)
		
		# select random principal components
		w <- which(PCA$sdev > 0)
		w <- sample(w, prob=PCA$sdev[w])
		s <- cumsum(PCA$sdev[w])
		s <- s/s[length(s)]
		n <- which.max(s >= tol - 1e-8)
		if (n < 3)
			n <- 3
		w <- w[seq_len(n)]
		X <- PCA$x[, w]
		R <- t(PCA$rotation[, w])
		C <- PCA$center
	#} else {
	#	R <- matrix(rnorm(k*n, sd=1/sqrt(n)), k, n) # sqrt(rowSums(R^2)) ~= 1
	#	#R <- R/sqrt(rowSums(R^2))
	#	
	#	X <- y %*% R
	#	R <- t(R)
	#	C <- rep(0, k)
	#}
	
	# determine weights
	w <- base^(-1000*(x/min(x) - 1)) # downweight worse scores
	w <- w*.redundancy(X, which.max(x)) # downweight redundant scores
	w <- w/sum(w)
	means <- colSums(w*X)
	sdevs <- sqrt(colSums(w*t(t(X) - means)^2)/(1 - sum(w^2)))
	
	Y <- matrix(qnorm(runif(n*N, interval),
			means,
			runif(n*N, mult[1L], mult[2L])*sdevs,
			lower.tail=sample(c(TRUE, FALSE), n, replace=TRUE)),
		ncol=n,
		byrow=TRUE)
	Y <- Y %*% R
	Y <- t(t(Y) + C)
	
	Y <- exp(Y)
	Y[Y < 1e-8] <- 1e-8
	Y[Y > 20] <- 20
	Y
}

.support <- function(myClusters, Trees) {
	t1 <- .extractClades(myClusters)
	w <- which(myClusters[, 7:8, drop=FALSE] > 0)
	t1 <- t1[w]
	t1 <- lapply(t1, sort)
	dim <- nrow(myClusters) + 1L
	W <- which(lengths(Trees) > 0)
	support <- numeric(length(t1))
	for (j in seq_along(W)) {
		t2 <- .extractClades(Trees[[W[j]]])
		t2 <- lapply(t2, sort)
		t2 <- c(t2,
			lapply(t2,
				function(x)
					seq_len(dim)[-x]))
		support <- support + tail(duplicated(c(t2, t1)), length(t1))
	}
	support/length(W)
}

TreeLine <- function(myXStringSet=NULL,
	myDistMatrix=NULL,
	method="ML",
	type="dendrogram",
	model=MODELS,
	cutoff=-Inf,
	showPlot=FALSE,
	collapse=-1,
	reconstruct=FALSE,
	root=0,
	informationCriterion="AICc",
	maxGenerations=20,
	maxTime=Inf,
	quadrature=FALSE,
	costMatrix=1 - diag(length(x)),
	processors=1,
	verbose=TRUE) {
	
	# initialize variables
	time.1 <- Sys.time()
	
	# error checking
	if (length(method) != 1)
		stop("Only one method can be specified.")
	METHODS <- c("NJ","UPGMA", "ML", "complete", "single", "WPGMA", "MP")
	method <- pmatch(method, METHODS)
	if (is.na(method))
		stop("Invalid method.")
	if (method==-1)
		stop("Ambiguous method.")
	if (!is.logical(reconstruct) &&
		!is.numeric(reconstruct))
		stop("reconstruct must be a logical or numeric.")
	if (is.numeric(reconstruct)) {
		if (reconstruct <= 0)
			stop("reconstruct must be be greater than zero.")
		if (reconstruct > 1)
			stop("reconstruct can be at most one.")
		if (method == 7)
			stop("reconstruct must be a logical when method is 'MP'.")
	}
	TYPES <- c("clusters", "dendrogram", "both")
	type <- pmatch(type[1], TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	if (!is.numeric(cutoff))
		stop("cutoff must be a numeric.")
	if (is.integer(cutoff))
		cutoff <- as.numeric(cutoff)
	if (!is.logical(showPlot))
		stop("showPlot must be a logical.")
	if (length(cutoff) > 1 && type > 1)
		warning("Multiple cutoffs can only be supplied when type is 'clusters'.")
	ASC <- TRUE
	if (is.unsorted(cutoff)) {
		if (is.unsorted(rev(cutoff))) {
			stop("cutoff must be sorted.")
		} else {
			ASC <- FALSE
		}
	}
	if (!is.numeric(collapse))
		stop("collapse must be a numeric.")
	if (length(informationCriterion) != 1L)
		stop("Only one informationCriterion can be specified.")
	ICs <- c("AICc", "BIC")
	informationCriterion <- pmatch(informationCriterion, ICs)
	if (is.na(informationCriterion))
		stop("Invalid informationCriterion.")
	if (informationCriterion == -1)
		stop("Ambiguous informationCriterion.")
	informationCriterion <- ICs[informationCriterion]
	if (!is.numeric(maxGenerations))
		stop("maxGenerations must be a numeric.")
	if (maxGenerations <= 0)
		stop("maxGenerations must be greater than zero.")
	if (floor(maxGenerations) != maxGenerations)
		stop("maxGenerations must be a whole number.")
	maxGenerations <- as.integer(maxGenerations)
	if (length(maxTime) != 1)
		stop("maxTime must be a single numeric.")
	if (!is.numeric(maxTime))
		stop("maxTime must be a numeric.")
	if (!is.logical(quadrature))
		stop("quadrature must be a logical.")
	if (method == 7) {
		sexpr <- substitute(costMatrix)
		if (!is.numeric(sexpr)) {
			if (is.name(sexpr)) { # function name
				sexpr <- call(as.character(sexpr),
					as.name("x")) # pass in 'x'
			} else if (!(is.call(sexpr) && # call
				"x" %in% all.vars(sexpr))) { # containing 'x'
				stop("costMatrix must be a call containing 'x'.")
			}
		}
	}
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (!is.null(processors) && !is.numeric(processors))
		stop("processors must be a numeric.")
	if (!is.null(processors) && floor(processors) != processors)
		stop("processors must be a whole number.")
	if (!is.null(processors) && processors < 1)
		stop("processors must be at least 1.")
	if (is.null(processors)) {
		processors <- detectCores()
	} else {
		processors <- as.integer(processors)
	}
	
	if (!((method == 3 || method == 7) &&
		is.null(myDistMatrix))) {
		if (is(myDistMatrix, "matrix")) {
			dim <- dim(myDistMatrix)
			if (dim[2] != dim[1])
				stop("myDistMatrix is not square.")
			dim <- dim[1]
		} else if (is(myDistMatrix, "dist")) {
			dim <- attr(myDistMatrix, "Size")
		} else {
			stop(paste("myDistMatrix must be a matrix for method '", METHODS[method], "'.", sep=""))
		}
		if (dim < 2)
			stop("myDistMatrix is too small.")
		if (typeof(myDistMatrix)=="integer")
			myDistMatrix[] <- as.numeric(myDistMatrix)
	} else {
		dim <- length(myXStringSet)
	}
	if (type > 1) {
		if (!is.numeric(root))
			stop("root must be a numeric.")
		if (length(root) != 1)
			stop("root must be a single numeric.")
		if (floor(root) != root)
			stop("root must be an integer.")
		if (root < 0)
			stop("root must be at least 0.")
		if (root > dim)
			stop(paste("root cannot be greater than ", dim, ".", sep=""))
	}
	if (method == 3 ||
		method == 7 ||
		(reconstruct && type > 1)) {
		if (is(myXStringSet, "DNAStringSet")) {
			typeX <- 1L
		} else if (is(myXStringSet, "RNAStringSet")) {
			typeX <- 2L
		} else if (is(myXStringSet, "AAStringSet")) {
			typeX <- 3L
		} else {
			stop("myXStringSet must be an AAStringSet, DNAStringSet, or RNAStringSet.")
		}
		if (length(myXStringSet) != dim)
			stop("myDistMatrix must have as many rows as the number of sequences.")
		if (length(unique(width(myXStringSet))) != 1)
			stop("All sequences in myXStringSet must be the same width (aligned).")
	}
	
	orgXStringSet <- myXStringSet # keep unmasked for reconstruction
	absTol <- 1 # (initial) absolute convergence tolerance
	relTol <- 0.001 # (initial) relative convergence tolerance
	
	if (method == 3L || # ML needs a model
		(type > 1 && # return a dendrogram
		reconstruct && # reconstruction uses a model
		method != 7L)) { # MP reconstruction does not use a model
		if (is.list(model)) {
			if (typeX == 3L) { # amino acids
				model <- model$Protein
				if (is.null(model))
					stop("model must contain a list element named 'Protein'.")
				if (!is.character(model))
					stop("The list element 'Protein' in model must be a character vector.")
			} else {
				model <- model$Nucleotide
				if (is.null(model))
					stop("model must contain a list element named 'Nucleotide'.")
				if (!is.character(model))
					stop("The list element 'Nucleotide' in model must be a character vector.")
			}
		} else if (!is.character(model)) {
			stop("model must be a list or character vector.")
		}
		if (length(model) == 0L)
			stop("No model(s) specified.")
		if (any(is.na(model)))
			stop("model cannot be NA.")
		model <- unique(model)
		submodels <- model
		indels <- grepl("\\+Indels?", submodels, ignore.case=TRUE)
		if (any(indels)) {
			if (!all(indels))
				stop("Models with indels cannot be compared to indel-free models.")
			submodels <- gsub("\\+Indels?", "", submodels, ignore.case=TRUE)
			indels <- TRUE
			myXStringSet <- .mask(myXStringSet)
		} else {
			indels <- FALSE
		}
		rates <- sub("([^+]*)(\\+G(\\d+))?", "\\3", submodels, ignore.case=TRUE)
		submodels <- sub("([^+]*)(\\+G(\\d+))?", "\\1", submodels, ignore.case=TRUE)
		empirical <- grepl("\\+F", submodels, ignore.case=TRUE)
		if (any(empirical)) {
			submodels <- gsub("\\+F", "", submodels, ignore.case=TRUE)
			if (any(empirical & (submodels %in% c("JC69", "K80", "SYM"))))
				stop("model cannot be 'JC69', 'K80', or 'SYM' with '+F'.")
		}
		
		if (typeX == 3L) { # amino acids
			MODELS <- MODELS$Protein
			ProtModels <- matrix(c(0.1784266, 0.0929129, 0.782913, 1.241095, 0.05795374, 7.185182, 0.008929181, 0.1821885, 1.374268e-06, 0.02340019, 0.1992269, 1.923901, 0.08705989, 0.1843856, 1.046446e-08, 0.9521821, 0.06273863, 0.5038373, 7.426619, 7.519215e-11, 3.691671, 1.851951, 1.0894, 0.4868901, 2.1124, 0.05891123, 0.0551634, 1.38937, 5.241316, 10.4955, 14.05444, 11.26995, 3.963388, 8.908434, 7.29808, 9.139518, 0.1140412, 0.3245175, 1.762721, 0.03916999, 0.0006594967, 6.712736e-06, 0.0001029959, 0.03560482, 4.706586, 0.06969101, 0.3932002, 0.02769442, 0.03020502, 0.006079219, 0.6802781, 0.001283121, 0.02157936, 5.879103, 1.601123, 0.07388355, 7.54924, 6.190318, 0.06622772, 3.722878e-16, 3.030805, 3.608816, 0.055044, 1.455741, 0.5059793, 0.02158451, 0.06299271, 0.3362326, 0.03972173, 0.03357577, 0.007213178, 0.001233336, 0.07659566, 0.02187264, 2.298295, 10.96748, 5.647985, 1.238634, 0.1130146, 0.08208677, 0.1955446, 0.1031734, 0.1993818, 0.00149661, 0.05288625, 0.1984772, 5.642309, 2.714705, 3.390618, 0.004649035, 3.94794, 1.800713, 0.3498713, 0.007342554, 0.1509482, 0.004878395, 0.7426909, 0.02889815, 0.07915056, 10.49496, 0.05016568, 1.149931, 0.009948994, 0.07417279, 0.3556198, 0.9988358, 1.926435, 7.348346, 0.5822988, 0.2639482, 0.0005906405, 0.06776709, 0.9984215, 5.439116, 0.6007607, 0.1580539, 0.08688405, 0.01861354, 0.9813064, 1.284651, 2.912317, 1.135258, 2.147175, 0.1516881, 3.225214e-06, 0.1202094, 0.06016624, 0.07862767, 3.443285, 3.087152, 0.5702792, 1.039298, 1.415612, 0.03674486, 0.9057112, 3.058575, 0.07939549, 0.5724286, 0.0007310937, 0.01423897, 0.4440833, 4.332983e-05, 0.02252612, 0.1386853, 7.01389, 0.06318748, 0.3378544, 0.008024263, 0.1011149, 0.2199856, 0.005516074, 0.1385142, 0.01412361, 0.1433528, 0.1711315, 2.622763, 0.9078338, 0.7741612, 0.02737091, 0.1240642, 0.2295842, 20.55414, 0.2903165, 0.152132, 0.07109973, 0.002246759, 7.074464, 0.1992133, 0.8104751, 0.09984255, 0.6121284, 3.774477, 0.1366145, 0.04931206, 0.4076074, 0.02243512, 0.009047737, 0.5795409, 0.42282, 6.890244, 7.926675, 3.59531, 0.0349344, 4.39672, 1.643946, 0.2217442, 0.07477041, 0.2166054, 0.09663569, 0.5010635, 0.06541704, 0.04708366, 0.03168984, 0.04688141, 0.02150693, 0.04240711, 0.02842211, 0.1005278, 0.009812606, 0.03424424, 0.06222565, 0.04844488, 0.0176037, 0.03478555, 0.03962469, 0.1280566, 0.08199314, 0.03393045, 0.07586119, 0.04948141,
					0.735790389698, 0.485391055466, 1.297446705134, 0.543161820899, 0.500964408555, 3.180100048216, 1.45999531047, 0.227826574209, 0.397358949897, 0.240836614802, 1.199705704602, 3.020833610064, 1.839216146992, 1.190945703396, 0.32980150463, 1.1709490428, 1.36057419042, 1.24048850864, 3.761625208368, 0.140748891814, 5.528919177928, 1.95588357496, 0.418763308518, 1.355872344485, 0.798473248968, 0.418203192284, 0.609846305383, 0.423579992176, 0.716241444998, 1.456141166336, 2.414501434208, 0.778142664022, 0.354058109831, 2.43534113114, 1.626891056982, 0.539859124954, 0.605899003687, 0.232036445142, 0.283017326278, 0.418555732462, 0.774894022794, 0.236202451204, 0.186848046932, 0.189296292376, 0.252718447885, 0.800016530518, 0.622711669692, 0.211888159615, 0.218131577594, 0.831842640142, 0.580737093181, 0.372625175087, 0.217721159236, 0.348072209797, 3.890963773304, 1.295201266783, 5.411115141489, 1.593137043457, 1.032447924952, 0.285078800906, 3.945277674515, 2.802427151679, 0.752042440303, 1.022507035889, 0.406193586642, 0.445570274261, 1.253758266664, 0.983692987457, 0.648441278787, 0.222621897958, 0.76768882348, 2.494896077113, 0.55541539747, 0.459436173579, 0.984311525359, 3.364797763104, 6.030559379572, 1.073061184332, 0.492964679748, 0.371644693209, 0.354861249223, 0.281730694207, 0.441337471187, 0.14435695975, 0.291409084165, 0.368166464453, 0.714533703928, 1.517359325954, 2.064839703237, 0.266924750511, 1.77385516883, 1.173275900924, 0.448133661718, 0.494887043702, 0.730628272998, 0.356008498769, 0.858570575674, 0.926563934846, 0.504086599527, 0.527007339151, 0.388355409206, 0.374555687471, 1.047383450722, 0.454123625103, 0.233597909629, 4.325092687057, 1.12278310421, 2.904101656456, 1.582754142065, 1.197188415094, 1.934870924596, 1.769893238937, 1.509326253224, 1.11702976291, 0.35754441246, 0.352969184527, 1.752165917819, 0.918723415746, 0.540027644824, 1.169129577716, 1.729178019485, 0.914665954563, 1.898173634533, 0.934187509431, 1.119831358516, 1.277480294596, 1.071097236007, 0.641436011405, 0.585407090225, 1.17909119726, 0.915259857694, 1.303875200799, 1.488548053722, 0.488206118793, 1.005451683149, 5.15155629227, 0.465839367725, 0.426382310122, 0.191482046247, 0.145345046279, 0.527664418872, 0.758653808642, 0.407635648938, 0.508358924638, 0.30124860078, 0.34198578754, 0.6914746346, 0.332243040634, 0.888101098152, 2.074324893497, 0.252214830027, 0.387925622098, 0.513128126891, 0.718206697586, 0.720517441216, 0.538222519037, 0.261422208965, 0.470237733696, 0.95898974285, 0.596719300346, 0.308055737035, 4.218953969389, 0.674617093228, 0.811245856323, 0.7179934869, 0.951682162246, 6.747260430801, 0.369405319355, 0.796751520761, 0.801010243199, 4.054419006558, 2.187774522005, 0.438388343772, 0.312858797993, 0.258129289418, 1.116352478606, 0.530785790125, 0.524253846338, 0.25334079019, 0.20155597175, 8.311839405458, 2.231405688913, 0.498138475304, 2.575850755315, 0.838119610178, 0.496908410676, 0.561925457442, 2.253074051176, 0.266508731426, 1, 0.074, 0.052, 0.045, 0.054, 0.025, 0.034, 0.054, 0.074, 0.026, 0.068, 0.099, 0.058, 0.025, 0.047, 0.039, 0.057, 0.051, 0.013, 0.032, 0.073,
					105, 227, 357, 175, 43, 4435, 669, 823, 538, 10, 157, 1745, 768, 400, 10, 499, 152, 1055, 3691, 10, 3122, 665, 243, 653, 431, 303, 133, 379, 66, 715, 1405, 331, 441, 1269, 162, 19, 145, 136, 168, 10, 280, 92, 148, 40, 29, 197, 203, 113, 10, 396, 286, 82, 20, 66, 1745, 236, 4482, 2430, 412, 48, 3313, 2629, 263, 305, 345, 218, 185, 125, 61, 47, 159, 202, 113, 21, 10, 1772, 1351, 193, 68, 53, 97, 22, 726, 10, 145, 25, 127, 454, 1268, 72, 327, 490, 87, 173, 170, 285, 323, 185, 28, 152, 117, 219, 302, 100, 43, 2440, 385, 2085, 590, 2331, 396, 568, 691, 303, 216, 516, 868, 93, 487, 1202, 1340, 314, 1393, 266, 576, 241, 369, 92, 32, 1040, 156, 918, 645, 148, 260, 2151, 14, 230, 40, 18, 435, 53, 63, 82, 69, 42, 159, 10, 86, 468, 49, 73, 29, 56, 323, 754, 281, 1466, 391, 142, 10, 1971, 89, 189, 247, 215, 2370, 97, 522, 71, 346, 968, 92, 83, 75, 592, 54, 200, 91, 25, 4797, 865, 249, 475, 317, 122, 167, 760, 10, 119, 0.0755, 0.0621, 0.041, 0.0371, 0.0091, 0.0382, 0.0495, 0.0838, 0.0246, 0.0806, 0.1011, 0.0504, 0.022, 0.0506, 0.0431, 0.0622, 0.0543, 0.0181, 0.0307, 0.066,
					6.5, 4.5, 10.6, 84.3, 9.5, 643.2, 19.5, 353.7, 10.9, 10.7, 6.1, 486.3, 18, 11.6, 0.1, 74.5, 21.5, 13, 437.4, 0.1, 342.6, 118.1, 183.9, 17.4, 150.3, 86.8, 7.1, 161.9, 2.8, 346.6, 345.3, 202.4, 111.8, 450.1, 6.2, 2.2, 1.5, 50.6, 25.6, 5.6, 3.4, 3.6, 4.3, 2.5, 8.4, 3.9, 36.9, 2.4, 5.9, 20.3, 26.1, 5.1, 3.4, 17.3, 205, 4.2, 712.1, 639.2, 10.1, 0.1, 500.5, 426.6, 29.3, 9.2, 37.9, 10.8, 13.4, 53.5, 9.9, 3.8, 10.5, 9.5, 9.6, 3.8, 3.6, 534.9, 142.8, 83.6, 4.3, 5, 8.7, 7.5, 238, 2.4, 7.7, 3.1, 11, 61, 542.3, 9.4, 3.8, 91.2, 69, 3.5, 13.4, 6.5, 145.6, 8.1, 2.6, 133.9, 2.1, 155.8, 21.2, 10.5, 12.6, 251.1, 82.9, 271.4, 34.8, 471.9, 10.7, 16.4, 136.7, 19.2, 36.2, 160.3, 23.9, 6.2, 249.4, 348.6, 467.5, 82.5, 215.5, 8, 7.4, 5.4, 11.6, 6.3, 3.8, 266.2, 10.7, 140.2, 295.2, 3.6, 181.2, 144.8, 3.4, 171.8, 6.1, 3.5, 518.6, 17, 9.1, 49, 5.7, 3.3, 98.8, 2.3, 11.1, 34.1, 1.1, 56.3, 1.5, 2.2, 4.3, 69.9, 202.9, 579.1, 9.4, 9.1, 2.1, 889.2, 10.8, 9.6, 20.1, 3.4, 255.9, 5.6, 264.3, 3.3, 21.7, 363.2, 8.4, 1.6, 10.3, 37.8, 5.1, 21.6, 76, 1.1, 595, 155.8, 9.2, 191.9, 102.2, 7.7, 10.1, 36.8, 5, 7.2, 0.061007, 0.060799, 0.043028, 0.038515, 0.011297, 0.035406, 0.050764, 0.073749, 0.024609, 0.085629, 0.10693, 0.046704, 0.023382, 0.056136, 0.043289, 0.073994, 0.052078, 0.018023, 0.036043, 0.05862,
					27, 98, 32, 120, 0, 905, 36, 23, 0, 0, 89, 246, 103, 134, 0, 198, 1, 148, 1153, 0, 716, 240, 9, 139, 125, 11, 28, 81, 23, 240, 535, 86, 28, 606, 43, 10, 65, 64, 77, 24, 44, 18, 61, 0, 7, 41, 15, 34, 0, 0, 73, 11, 7, 44, 257, 26, 464, 318, 71, 0, 153, 83, 27, 26, 46, 18, 72, 90, 1, 0, 0, 114, 30, 17, 0, 336, 527, 243, 18, 14, 14, 0, 0, 0, 0, 15, 48, 196, 157, 0, 92, 250, 103, 42, 13, 19, 153, 51, 34, 94, 12, 32, 33, 17, 11, 409, 154, 495, 95, 161, 56, 79, 234, 35, 24, 17, 96, 62, 46, 245, 371, 26, 229, 66, 16, 53, 34, 30, 22, 192, 33, 136, 104, 13, 78, 550, 0, 201, 23, 0, 0, 0, 0, 0, 27, 0, 46, 0, 0, 76, 0, 75, 0, 24, 8, 95, 0, 96, 0, 22, 0, 127, 37, 28, 13, 0, 698, 0, 34, 42, 61, 208, 24, 15, 18, 49, 35, 37, 54, 44, 889, 175, 10, 258, 12, 48, 30, 157, 0, 28, 0.087127, 0.040904, 0.040432, 0.046872, 0.033474, 0.038255, 0.04953, 0.088612, 0.033618, 0.036886, 0.085357, 0.080482, 0.014753, 0.039772, 0.05068, 0.069577, 0.058542, 0.010494, 0.029916, 0.064718,
					0.267828, 0.984474, 0.327059, 1.199805, 0, 8.931515, 0.360016, 0.232374, 0, 0, 0.887753, 2.439939, 1.028509, 1.348551, 0, 1.961167, 0, 1.493409, 11.388659, 0, 7.086022, 2.386111, 0.087791, 1.385352, 1.240981, 0.107278, 0.281581, 0.811907, 0.228116, 2.383148, 5.290024, 0.868241, 0.282729, 6.011613, 0.439469, 0.106802, 0.653416, 0.632629, 0.768024, 0.239248, 0.438074, 0.180393, 0.609526, 0, 0.076981, 0.406431, 0.154924, 0.341113, 0, 0, 0.730772, 0.11288, 0.071514, 0.443504, 2.556685, 0.258635, 4.610124, 3.148371, 0.716913, 0, 1.519078, 0.830078, 0.267683, 0.270475, 0.460857, 0.180629, 0.71784, 0.896321, 0, 0, 0, 1.127499, 0.304803, 0.170372, 0, 3.332732, 5.230115, 2.411739, 0.183641, 0.136906, 0.138503, 0, 0, 0, 0, 0.153478, 0.475927, 1.951951, 1.56516, 0, 0.92186, 2.48592, 1.028313, 0.419244, 0.13394, 0.18755, 1.526188, 0.507003, 0.347153, 0.933709, 0.119152, 0.316258, 0.335419, 0.170205, 0.110506, 4.05187, 1.53159, 4.885892, 0.956097, 1.598356, 0.561828, 0.793999, 2.322243, 0.353643, 0.247955, 0.171432, 0.954557, 0.619951, 0.459901, 2.427202, 3.680365, 0.265745, 2.271697, 0.66093, 0.162366, 0.525651, 0.340156, 0.306662, 0.226333, 1.900739, 0.33109, 1.350599, 1.031534, 0.136655, 0.782857, 5.436674, 0, 2.001375, 0.224968, 0, 0, 0, 0, 0, 0.270564, 0, 0.461776, 0, 0, 0.762354, 0, 0.740819, 0, 0.244139, 0.078012, 0.94694, 0, 0.953164, 0, 0.214717, 0, 1.2654, 0.374834, 0.286572, 0.132142, 0, 6.952629, 0, 0.336289, 0.417839, 0.60807, 2.059564, 0.240368, 0.158067, 0.178316, 0.484678, 0.346983, 0.36725, 0.538165, 0.438715, 8.810038, 1.745156, 0.10385, 2.565955, 0.123606, 0.485026, 0.303836, 1.561997, 0, 0.279379, 0.087127, 0.040904, 0.040432, 0.046872, 0.033474, 0.038255, 0.04953, 0.088612, 0.033619, 0.036886, 0.085357, 0.080481, 0.014753, 0.039772, 0.05068, 0.069577, 0.058542, 0.010494, 0.029916, 0.064718,
					0.531678, 0.557967, 0.451095, 0.827445, 0.154899, 5.54953, 0.574478, 1.019843, 0.313311, 0.105625, 0.556725, 3.021995, 0.768834, 0.521646, 0.091304, 1.066681, 0.318483, 0.578115, 7.766557, 0.053907, 3.417706, 1.740159, 1.359652, 0.773313, 1.272434, 0.546389, 0.231294, 1.115632, 0.21997, 3.210671, 4.025778, 1.032342, 0.724998, 5.68408, 0.243768, 0.201696, 0.361684, 0.239195, 0.491003, 0.115968, 0.150559, 0.07827, 0.111773, 0.053769, 0.181788, 0.310007, 0.372261, 0.137289, 0.061486, 0.164593, 0.709004, 0.097485, 0.069492, 0.540571, 2.335139, 0.369437, 6.529255, 2.529517, 0.282466, 0.049009, 2.966732, 1.731684, 0.26984, 0.525096, 0.202562, 0.146481, 0.469395, 0.431045, 0.33072, 0.190001, 0.409202, 0.456901, 0.175084, 0.130379, 0.32966, 4.831666, 3.856906, 0.624581, 0.138293, 0.065314, 0.073481, 0.032522, 0.678335, 0.045683, 0.043829, 0.050212, 0.453428, 0.77709, 2.500294, 0.024521, 0.436181, 1.959599, 0.710489, 0.121804, 0.127164, 0.123653, 1.608126, 0.191994, 0.208081, 1.141961, 0.09858, 1.060504, 0.216345, 0.164215, 0.148483, 3.887095, 1.001551, 5.057964, 0.589268, 2.155331, 0.548807, 0.312449, 1.874296, 0.743458, 0.405119, 0.592511, 0.474478, 0.285564, 0.943971, 2.788406, 4.582565, 0.650282, 2.351311, 0.425159, 0.469823, 0.523825, 0.331584, 0.316862, 0.477355, 2.553806, 0.272514, 0.965641, 2.114728, 0.138904, 1.176961, 4.777647, 0.084329, 1.257961, 0.0277, 0.057466, 1.104181, 0.172206, 0.114381, 0.54418, 0.128193, 0.13451, 0.530324, 0.089134, 0.201334, 0.537922, 0.069965, 0.310927, 0.080556, 0.139492, 0.235601, 0.700693, 0.453952, 2.114852, 0.254745, 0.063452, 0.0525, 5.8484, 0.303445, 0.241094, 0.087904, 0.18987, 5.484236, 0.11385, 0.628608, 0.201094, 0.747889, 2.924161, 0.171995, 0.164525, 0.315261, 0.621323, 0.179771, 0.465271, 0.47014, 0.121827, 9.533943, 1.761439, 0.124066, 3.038533, 0.593478, 0.211561, 0.408532, 1.14398, 0.239697, 0.165473, 0.076862, 0.051057, 0.042546, 0.051269, 0.020279, 0.041061, 0.06182, 0.074714, 0.022983, 0.052569, 0.091111, 0.059498, 0.023414, 0.04053, 0.050532, 0.068225, 0.058518, 0.014336, 0.032303, 0.066374,
					0.41136921, 0.272108441, 0.690697816, 0.537979882, 0.128916229, 6.419184592, 1.927407726, 0.637641072, 0.427332366, 0.097266653, 1.000822747, 3.125742015, 1.740085998, 0.577858431, 0.068442448, 1.064958047, 0.364765559, 0.532651549, 6.209372754, 0.011334452, 4.340145193, 1.726700741, 0.580834919, 1.078463067, 0.938958085, 0.489548568, 0.230619324, 0.724478494, 0.341369311, 2.524906667, 4.818725162, 1.119743587, 0.595915643, 6.48313597, 0.401796545, 0.241270083, 0.160754676, 0.165240985, 0.363112134, 0.021046133, 0.274642348, 0.09567887, 0.050706819, 0.021111304, 0.146772416, 0.385451233, 0.314453376, 0.071487784, 0.02242492, 0.517303194, 0.779312312, 0.07924615, 0.05334276, 0.468341277, 3.884679638, 0.495045475, 8.155672927, 2.465965958, 0.291646536, 0.022531075, 3.349091, 2.182616219, 0.221198251, 0.599064533, 0.197661824, 0.138352381, 0.812936578, 0.43490116, 0.308117286, 0.034608487, 0.574523883, 1.345172045, 0.145134696, 0.088317968, 0.362877569, 4.63178066, 4.973633974, 0.607340896, 0.271938118, 0.05859432, 0.097921003, 0.032902914, 1.201541857, 0.050161785, 0.023936456, 0.088993297, 0.73647602, 1.301154678, 3.997048474, 0.03800583, 1.407808655, 1.279061375, 0.356075939, 0.154293877, 0.400956836, 0.075042846, 0.879248211, 0.376861606, 0.158645548, 0.748197796, 0.084385318, 0.514925623, 0.349497923, 0.083301449, 0.100796273, 4.636757164, 0.86902097, 5.073408819, 1.205097628, 2.403450536, 1.196393794, 0.557335893, 1.497680467, 0.93120842, 0.112057128, 0.478892595, 0.658694216, 0.263297738, 0.779492979, 1.949090428, 3.104455681, 0.506483243, 1.81902479, 0.36642353, 0.720288888, 0.888663974, 0.467561373, 0.090869516, 0.510002979, 2.209774633, 0.249127345, 0.939599909, 1.858374216, 0.154629354, 0.72029647, 5.295672294, 0.127904694, 0.503935361, 0.030958389, 0.031968231, 1.183140601, 0.176037651, 0.070145131, 0.222343146, 0.374648992, 0.078073244, 0.436140848, 0.045734005, 0.411607444, 1.578424111, 0.078951166, 0.222802108, 0.098596321, 0.241089433, 0.320273276, 1.003153977, 0.296230395, 1.487182263, 0.323840895, 0.146012022, 0.062432558, 11.70025569, 0.254715265, 0.340819545, 0.135388887, 0.418123202, 9.692174642, 0.098514909, 0.542362405, 0.229175918, 2.071201519, 3.76162852, 0.164533836, 0.096783732, 0.11352745, 1.472220719, 0.226873568, 0.346774522, 0.207470551, 0.105889494, 12.05200144, 1.735555404, 0.167239757, 2.334685771, 0.810584804, 0.267492174, 0.114366882, 1.637649103, 0.160563059, 0.263012065, 0.0698408, 0.0527365, 0.0386741, 0.0434241, 0.0184966, 0.0342301, 0.0625993, 0.079547, 0.0206979, 0.0599927, 0.0926412, 0.0627311, 0.0329818, 0.0337404, 0.0441391, 0.0628714, 0.0706965, 0.0265006, 0.0240127, 0.069446,
					0.077462, 0.078037, 2e-05, 0.550515, 0.089476, 8.801355, 0.114675, 0.572845, 2e-05, 2e-05, 0.09049, 3.856678, 0.093133, 0.183601, 2e-05, 0.560685, 0.020614, 2e-05, 7.603314, 2e-05, 1.058066, 0.478097, 1.281339, 0.147801, 1.359685, 0.221287, 2e-05, 2.331269, 0.168719, 2.704348, 3.326034, 0.543235, 0.240631, 8.958548, 0.025386, 0.032838, 2e-05, 0.15125, 0.538307, 2e-05, 2e-05, 2e-05, 2e-05, 2e-05, 0.037032, 0.002576, 0.114188, 2e-05, 2e-05, 0.085459, 0.789751, 0.023842, 0.037784, 0.375616, 1.571343, 2e-05, 19.413475, 2.764154, 2e-05, 2e-05, 1.545589, 2.99296, 0.027884, 2e-05, 0.253247, 0.03856, 2e-05, 0.151495, 2e-05, 2e-05, 2e-05, 2e-05, 0.031428, 0.009423, 2e-05, 8.538171, 2.895172, 0.340189, 2e-05, 2e-05, 2e-05, 2e-05, 0.773082, 2e-05, 2e-05, 2e-05, 0.092713, 1.368786, 6.296781, 2e-05, 2e-05, 1.256282, 0.292376, 2e-05, 2e-05, 2e-05, 1.571014, 2e-05, 0.029602, 0.723179, 2e-05, 1.182805, 0.097941, 2e-05, 0.035125, 2.104417, 0.947496, 9.494063, 0.127754, 1.516267, 0.140959, 2e-05, 1.735582, 0.079928, 0.211717, 1.00888, 0.131204, 0.073511, 1.821906, 4.417267, 6.315368, 0.199572, 1.32982, 2e-05, 2e-05, 0.057328, 2e-05, 0.044927, 0.144097, 6.767791, 0.033401, 0.567475, 2.228245, 2e-05, 1.543179, 3.290965, 2e-05, 0.791312, 2e-05, 2e-05, 2.305902, 2e-05, 2e-05, 0.323982, 0.206837, 2e-05, 0.176249, 2e-05, 0.173206, 0.148283, 0.074193, 0.154086, 2e-05, 0.043583, 2e-05, 1.247952, 0.771837, 4.085694, 2e-05, 0.029684, 2e-05, 16.468674, 2e-05, 2e-05, 0.071741, 0.228965, 5.907642, 2e-05, 0.431256, 2e-05, 2e-05, 8.030561, 2e-05, 2e-05, 0.439809, 0.149716, 2e-05, 0.346045, 0.583259, 0.067473, 16.913225, 1.074192, 0.006629, 3.855848, 1.176807, 0.03675, 0.052613, 0.134211, 2e-05, 0.234116, 0.0775, 0.053813, 0.03395, 0.034973, 0.014056, 0.030139, 0.054825, 0.086284, 0.01821, 0.063272, 0.103857, 0.059646, 0.040389, 0.03363, 0.036649, 0.060915, 0.076327, 0.030152, 0.020069, 0.071343,
					0.138658764751059, 0.0533665787145181, 0.161000889039552, 0.584852305649886, 0.00677184253227681, 7.73739287051356, 0.0264470951166826, 0.16720700818221, 1.30249856764315e-05, 0.014132062548787, 0.353753981649393, 3.29271694159791, 0.530642655337477, 0.145469388422239, 0.00254733397966779, 1.4842345032161, 0.124897616909194, 0.0616521921873234, 5.37051127867923, 3.91106992668137e-11, 1.19562912226203, 1.13231312248046, 1.19062446519178, 0.322524647863997, 1.93483278448943, 0.116941459124876, 0.108051341246072, 1.59309882471598, 0.214757862168721, 1.87956993845887, 1.38709603234116, 0.887570549414031, 0.0218446166959521, 5.33031341222104, 0.256491863423002, 0.0587745274250666, 0.149926734229061, 0.246117171830255, 0.21857197541607, 0.0140859174993809, 0.00111215807314139, 0.0288399502994541, 0.0142107118685268, 1.62662283098296e-05, 0.243190142026506, 0.0231169515264061, 0.296045557460629, 0.000835873174542931, 0.00573068208525287, 0.00561362724916376, 1.02036695531654, 0.016499535540562, 0.00651622937676521, 0.321611693603646, 3.51207228207807, 0.474333610192982, 15.3000966197798, 2.6468479652886, 0.290042980143818, 3.83228119049152e-06, 2.559587177122, 3.88148880863814, 0.264148929349066, 0.347302791211758, 0.227707997165566, 0.129223639195248, 0.0587454231508643, 0.890162345593224, 0.00525168778853117, 0.0417629637305017, 0.111457310321926, 0.190259181297527, 0.313974351356074, 0.00150046692269255, 0.00127350890508147, 9.01795420287895, 6.74693648486614, 1.33129161941264, 0.0804909094320368, 0.0160550314767596, 0.000836445615590923, 1.0600102849456e-06, 0.10405366623526, 0.0326806570137471, 0.00100350082518749, 0.00123664495412902, 0.119028506158521, 1.46335727834648, 2.98680003596399, 0.319895904499071, 0.279910508981581, 0.659311477863896, 0.154027179890711, 0.0364417719063219, 0.188539456415654, 1.59312060172652e-13, 0.712769599068934, 0.319558828428154, 0.0386317614553493, 0.924466914225534, 0.0805433268150369, 0.634308520867322, 0.195750631825315, 0.0568693216513547, 0.0071324304661639, 3.01134451903854, 0.950138410087378, 3.88131053061457, 0.338372183381345, 0.336263344504404, 0.487822498528951, 0.307140298031341, 1.58564657669139, 0.580704249811294, 0.290381075260226, 0.570766693213698, 0.283807671568883, 0.00702658828739369, 0.996685669575839, 2.08738534433198, 5.4182981753166, 0.183076905018197, 2.14033231636063, 0.135481232622983, 0.011975265782196, 0.60234096342392, 0.2801248951174, 0.0188080299490973, 0.368713573381758, 2.90405228596936, 0.0449263566753846, 1.52696419998775, 2.03151132062208, 0.000134906239484254, 0.54225109402693, 2.2068599339404, 0.195966354027106, 1.36942940801512, 0.000536284040016542, 1.4893873721753e-05, 0.0941066800969967, 0.0440205200833047, 0.155245492137294, 0.196486447133033, 0.0223729191088972, 0.0321321499585514, 0.431277662888057, 4.97641445484395e-05, 0.0704600385245663, 0.814753093809928, 0.000431020702277328, 0.0998357527014247, 0.207066205546908, 0.0182892882245349, 0.0998554972524385, 0.373101926513925, 0.525398542949365, 0.601692431136271, 0.0722059354079545, 0.104092870343653, 0.0748149970972622, 6.44895444648517, 0.273934263183281, 0.340058468374384, 0.0124162215506117, 0.874272174533394, 5.39392424532822, 0.000182294881489116, 0.392552239890831, 0.124898020409882, 0.42775543040588, 3.53200526987468, 0.103964386383736, 0.0102575172450253, 0.297123975243582, 0.0549045639492389, 0.406697814049488, 0.285047948309311, 0.337229618868315, 0.0986313546653266, 14.3940521944257, 0.890598579382591, 0.0731279296372675, 4.90484223478739, 0.592587985458668, 0.0589719751511691, 0.0882564232979724, 0.654109108255219, 0.256900461407996, 0.167581646770807, 0.0470718, 0.0509102, 0.0742143, 0.0478596, 0.0250216, 0.0333036, 0.0545874, 0.0763734, 0.0199642, 0.0671336, 0.0714981, 0.0567845, 0.0181507, 0.0304961, 0.0506561, 0.0884091, 0.0743386, 0.0185237, 0.0314741, 0.0632292,
					12.2, 59.1, 313, 87.7, 40.6, 5083.8, 699.4, 1867.9, 133.7, 125.3, 199.5, 2321.3, 827.8, 387.6, 195.1, 288.3, 51.9, 230.9, 3638.7, 102.4, 2965.3, 677.9, 391.2, 446.5, 464.8, 433.6, 40.8, 393.8, 32, 1228.7, 2156.8, 527.7, 574.8, 1889.6, 144.6, 35.6, 39.5, 105.4, 146, 16.8, 96, 48.5, 23.4, 11.8, 38.1, 33.8, 124.9, 11.6, 8.7, 446.5, 235.4, 26.8, 4.3, 46.1, 1657.9, 78.3, 4666.4, 1930.8, 78.3, 89.2, 2571.2, 1870.9, 193.4, 96.9, 125.7, 37.6, 197.1, 112.1, 41.4, 33.1, 288.5, 341.8, 33.3, 11.8, 54.5, 2493.4, 1389.6, 220.1, 59.3, 26.1, 13.4, 15.9, 1540.1, 10.7, 30.7, 10.1, 87.2, 354.4, 1585.8, 9, 183.4, 505.8, 89.4, 18.7, 43.6, 124.9, 223.1, 39.4, 11.8, 185.1, 63.1, 347.7, 49.8, 19.2, 82.1, 2443.9, 269.6, 2505.4, 142.1, 3029.8, 246.9, 156.7, 709, 277, 75.6, 584.9, 209.3, 32.9, 790.5, 1154.5, 1646.2, 173.4, 1373.3, 66.1, 516.8, 128.5, 230, 65.5, 77.1, 1227.2, 62.9, 648.5, 1143.5, 29.8, 322.8, 2042.5, 30, 257.8, 19.9, 35, 562.6, 24.6, 57.2, 92.6, 58.5, 21.1, 218.2, 15.9, 70.5, 558.3, 18.5, 102.6, 22, 27.9, 137, 316.2, 272.1, 1571.2, 179.1, 79.2, 20.6, 4117.5, 35.8, 82, 84.2, 41.1, 2360.2, 15.8, 501.9, 64.6, 296.3, 1143.8, 62, 27.6, 18.5, 722.2, 50.2, 103.8, 72.6, 22.8, 6006.5, 661.7, 66.2, 465.7, 180.6, 75.5, 61.2, 583.7, 10.8, 37.3, 0.07951, 0.056001, 0.040459, 0.03322, 0.009051, 0.037505, 0.049675, 0.080233, 0.02188, 0.080496, 0.107512, 0.049324, 0.020776, 0.047731, 0.039916, 0.07382, 0.053615, 0.016705, 0.03079, 0.071781,
					0.307507, 0.005, 0.295543, 1.45504, 0.005, 17.6612, 0.123758, 0.351721, 0.0860642, 0.005, 0.0551128, 3.4215, 0.672052, 0.005, 0.005, 1.48135, 0.0749218, 0.0792633, 10.5872, 0.005, 2.5602, 2.13536, 3.65345, 0.323401, 2.83806, 0.897871, 0.0619137, 3.92775, 0.0847613, 9.04044, 7.64585, 1.9169, 0.240073, 7.05545, 0.11974, 0.005, 0.005, 0.677289, 0.680565, 0.0176792, 0.005, 0.005, 0.00609079, 0.005, 0.103111, 0.215256, 0.701427, 0.005, 0.00876048, 0.129777, 1.49456, 0.005, 0.005, 1.74171, 5.95879, 0.005, 20.45, 7.90443, 0.005, 0.005, 6.54737, 4.61482, 0.521705, 0.005, 0.322319, 0.0814995, 0.0186643, 2.51394, 0.005, 0.005, 0.005, 0.303676, 0.175789, 0.005, 0.005, 11.2065, 5.31961, 1.28246, 0.0141269, 0.005, 0.005, 0.005, 9.29815, 0.005, 0.005, 0.291561, 0.145558, 3.39836, 8.52484, 0.0342658, 0.188025, 2.12217, 1.28355, 0.00739578, 0.0342658, 0.005, 4.47211, 0.0120226, 0.005, 2.45318, 0.0410593, 2.07757, 0.0313862, 0.005, 0.005, 2.46633, 3.4791, 13.1447, 0.52823, 4.69314, 0.116311, 0.005, 4.38041, 0.382747, 1.21803, 0.927656, 0.504111, 0.005, 0.956472, 5.37762, 15.9183, 2.86868, 6.88667, 0.274724, 0.739969, 0.243589, 0.289774, 0.369615, 0.711594, 8.61217, 0.0437673, 4.67142, 4.94026, 0.0141269, 2.01417, 8.93107, 0.005, 0.991338, 0.005, 0.005, 2.63277, 0.026656, 0.005, 1.21674, 0.0695179, 0.005, 0.748843, 0.005, 0.089078, 0.829343, 0.0444506, 0.0248728, 0.005, 0.005, 0.00991826, 1.76417, 0.674653, 7.57932, 0.113033, 0.0792633, 0.005, 18.6943, 0.148168, 0.111986, 0.005, 0.005, 15.34, 0.0304381, 0.648024, 0.105652, 1.28022, 7.61428, 0.0812454, 0.026656, 1.04793, 0.420027, 0.0209153, 1.02847, 0.953155, 0.005, 17.7389, 1.41036, 0.265829, 6.8532, 0.723274, 0.005, 0.0749218, 0.709226, 0.005, 0.0410593, 0.060490222, 0.066039665, 0.044127815, 0.042109048, 0.020075899, 0.053606488, 0.071567447, 0.072308239, 0.022293943, 0.069730629, 0.098851122, 0.056968211, 0.019768318, 0.028809447, 0.046025282, 0.05060433, 0.053636813, 0.033011601, 0.028350243, 0.061625237,
					0.0744808, 0.617509, 0.16024, 4.43521, 0.0674539, 29.4087, 0.167653, 2.86364, 0.0604932, 0.005, 0.005, 10.6746, 0.342068, 0.005, 0.005, 5.56325, 0.0251632, 0.201526, 12.1233, 0.005, 3.20656, 1.8685, 13.4379, 0.0604932, 10.3969, 0.0489798, 0.0604932, 14.7801, 0.005, 6.84405, 8.59876, 2.31779, 0.005, 18.5465, 0.005, 0.005, 0.005, 1.34069, 0.987028, 0.145124, 0.005, 0.0342252, 0.0390512, 0.005, 0.005, 0.16024, 0.586757, 0.005, 0.005, 0.005, 2.89048, 0.129839, 0.0489798, 1.76382, 9.10246, 0.592784, 39.8897, 10.6655, 0.894313, 0.005, 13.0705, 23.9626, 0.279425, 0.22406, 0.817481, 0.005, 0.005, 3.28652, 0.201526, 0.005, 0.005, 0.005, 0.005, 0.0489798, 0.005, 17.3064, 11.3839, 4.09564, 0.597923, 0.005, 0.005, 0.005, 0.362959, 0.005, 0.005, 0.005, 0.005, 1.48288, 7.48781, 0.005, 0.005, 1.00981, 0.404723, 0.344848, 0.005, 0.005, 3.04502, 0.005, 0.005, 13.9444, 0.005, 9.83095, 0.111928, 0.005, 0.0342252, 8.5942, 8.35024, 14.5699, 0.427881, 1.12195, 0.16024, 0.005, 6.27966, 0.725157, 0.740091, 6.14396, 0.005, 0.392575, 4.27939, 14.249, 24.1422, 0.928203, 4.54206, 0.630395, 0.005, 0.203091, 0.458743, 0.0489798, 0.95956, 9.36345, 0.005, 4.04802, 7.41313, 0.114512, 4.33701, 6.34079, 0.005, 5.96564, 0.005, 0.005, 5.49894, 0.0443298, 0.005, 2.8258, 0.005, 0.005, 1.37031, 0.005, 0.005, 0.005, 0.005, 1.10156, 0.005, 0.005, 0.005, 5.06475, 2.28154, 8.34835, 0.005, 0.005, 0.005, 47.4889, 0.114512, 0.005, 0.005, 0.579198, 4.12728, 0.005, 0.933142, 0.490608, 0.005, 24.8094, 0.279425, 0.0744808, 2.91786, 0.005, 0.005, 2.19952, 2.79622, 0.827479, 24.8231, 2.95344, 0.128065, 14.7683, 2.28, 0.005, 0.862637, 0.005, 0.005, 1.35482, 0.0377494, 0.057321, 0.0891129, 0.0342034, 0.0240105, 0.0437824, 0.0618606, 0.0838496, 0.0156076, 0.0983641, 0.0577867, 0.0641682, 0.0158419, 0.0422741, 0.0458601, 0.0550846, 0.0813774, 0.019597, 0.0205847, 0.0515639,
					58, 54, 45, 81, 16, 528, 56, 113, 34, 10, 57, 310, 86, 49, 9, 105, 29, 58, 767, 5, 323, 179, 137, 81, 130, 59, 26, 119, 27, 328, 391, 112, 69, 597, 26, 23, 36, 22, 47, 11, 17, 9, 12, 6, 16, 30, 38, 12, 7, 23, 72, 9, 6, 56, 229, 35, 646, 263, 26, 7, 292, 181, 27, 45, 21, 14, 54, 44, 30, 15, 31, 43, 18, 14, 33, 479, 388, 65, 15, 5, 10, 4, 78, 4, 5, 5, 40, 89, 248, 4, 43, 194, 74, 15, 15, 14, 164, 18, 24, 115, 10, 102, 21, 16, 17, 378, 101, 503, 59, 223, 53, 30, 201, 73, 40, 59, 47, 29, 92, 285, 475, 64, 232, 38, 42, 51, 32, 33, 46, 245, 25, 103, 226, 12, 118, 477, 9, 126, 8, 4, 115, 18, 10, 55, 8, 9, 52, 10, 24, 53, 6, 35, 12, 11, 20, 70, 46, 209, 24, 7, 8, 573, 32, 24, 8, 18, 536, 10, 63, 21, 71, 298, 17, 16, 31, 62, 20, 45, 47, 11, 961, 180, 14, 323, 62, 23, 38, 112, 25, 16, 0.076748, 0.051691, 0.042645, 0.051544, 0.019803, 0.040752, 0.06183, 0.073152, 0.022944, 0.053761, 0.091904, 0.058676, 0.023826, 0.040126, 0.050901, 0.068765, 0.058565, 0.014261, 0.032102, 0.066005,
					0.425093, 0.276818, 0.751878, 0.395144, 0.123954, 5.076149, 2.489084, 0.534551, 0.528768, 0.062556, 0.969894, 2.807908, 1.695752, 0.523386, 0.084808, 1.038545, 0.36397, 0.541712, 5.24387, 0.003499, 4.128591, 2.06604, 0.390192, 1.437645, 0.844926, 0.569265, 0.267959, 0.348847, 0.358858, 2.426601, 4.509238, 0.927114, 0.640543, 4.813505, 0.423881, 0.311484, 0.14983, 0.126991, 0.191503, 0.01069, 0.320627, 0.072854, 0.044265, 0.008705, 0.108882, 0.395337, 0.301848, 0.068427, 0.015076, 0.594007, 0.582457, 0.069673, 0.044261, 0.366317, 4.145067, 0.536518, 6.326067, 2.145078, 0.282959, 0.013266, 3.234294, 1.807177, 0.296636, 0.697264, 0.159069, 0.1375, 1.124035, 0.484133, 0.371004, 0.025548, 0.89368, 1.672569, 0.173735, 0.139538, 0.442472, 4.273607, 6.312358, 0.656604, 0.253701, 0.052722, 0.089525, 0.017416, 1.105251, 0.035855, 0.018811, 0.089586, 0.682139, 1.112727, 2.592692, 0.023918, 1.798853, 1.177651, 0.332533, 0.161787, 0.394456, 0.075382, 0.624294, 0.419409, 0.196961, 0.508851, 0.078281, 0.24906, 0.390322, 0.099849, 0.094464, 4.727182, 0.858151, 4.008358, 1.240275, 2.784478, 1.223828, 0.611973, 1.73999, 0.990012, 0.064105, 0.182287, 0.748683, 0.34696, 0.361819, 1.338132, 2.139501, 0.578987, 2.000679, 0.42586, 1.14348, 1.080136, 0.604545, 0.129836, 0.584262, 1.033739, 0.302936, 1.136863, 2.020366, 0.165001, 0.571468, 6.472279, 0.180717, 0.593607, 0.045376, 0.02989, 0.670128, 0.236199, 0.077852, 0.268491, 0.597054, 0.11166, 0.619632, 0.049906, 0.696175, 2.457121, 0.095131, 0.248862, 0.140825, 0.218959, 0.31444, 0.612025, 0.135107, 1.165532, 0.257336, 0.120037, 0.054679, 5.306834, 0.232523, 0.299648, 0.131932, 0.481306, 7.803902, 0.089613, 0.400547, 0.245841, 3.151815, 2.54787, 0.170887, 0.083688, 0.037967, 1.959291, 0.210332, 0.245034, 0.076701, 0.119013, 10.649107, 1.702745, 0.185202, 1.898718, 0.654683, 0.296501, 0.098369, 2.188158, 0.18951, 0.249313, 0.079066, 0.055941, 0.041977, 0.053052, 0.012937, 0.040767, 0.071586, 0.057337, 0.022355, 0.062157, 0.099081, 0.0646, 0.022951, 0.042302, 0.04404, 0.061197, 0.053287, 0.012066, 0.034155, 0.069147,
					0.2, 0.2, 0.2, 1, 4, 500, 254, 36, 98, 11, 0.2, 154, 262, 0.2, 0.2, 0.2, 0.2, 183, 862, 0.2, 262, 200, 0.2, 121, 12, 81, 3, 44, 0.2, 41, 180, 0.2, 12, 314, 15, 0.2, 26, 2, 21, 7, 63, 11, 7, 3, 0.2, 4, 2, 13, 1, 79, 16, 2, 1, 6, 515, 0.2, 209, 467, 2, 0.2, 349, 106, 0.2, 0.2, 3, 4, 121, 5, 79, 0.2, 312, 67, 0.2, 56, 0.2, 515, 885, 106, 13, 5, 20, 0.2, 184, 0.2, 0.2, 1, 14, 118, 263, 11, 322, 49, 0.2, 17, 0.2, 0.2, 39, 8, 0.2, 1, 0.2, 12, 17, 5, 15, 673, 3, 398, 44, 664, 52, 31, 226, 11, 7, 8, 144, 112, 36, 87, 244, 0.2, 166, 0.2, 183, 44, 43, 0.2, 19, 204, 48, 70, 289, 14, 47, 660, 0.2, 0.2, 8, 0.2, 22, 7, 11, 2, 0.2, 0.2, 21, 16, 71, 54, 0.2, 2, 0.2, 1, 4, 251, 0.2, 72, 87, 8, 9, 191, 12, 20, 117, 71, 792, 18, 30, 46, 38, 340, 0.2, 23, 0.2, 350, 0.2, 14, 3, 0.2, 1855, 85, 26, 281, 52, 32, 61, 544, 0.2, 2, 0.054116, 0.018227, 0.039903, 0.02016, 0.009709, 0.018781, 0.024289, 0.068183, 0.024518, 0.092638, 0.148658, 0.021718, 0.061453, 0.088668, 0.041826, 0.09103, 0.049194, 0.029786, 0.039443, 0.057701,
					0.579999479, 0.539999739, 0.449999033, 0.810000261, 0.16, 5.28000938, 0.560000261, 1.12999942, 0.34, 0.1, 0.570000521, 3.099998065, 0.86, 0.490000776, 0.09, 1.05, 0.290000193, 0.58, 7.669990688, 0.05, 3.230000982, 1.790001042, 1.37000058, 0.810001172, 1.3, 0.590001515, 0.259999509, 1.19, 0.270000521, 3.279990714, 3.910001172, 1.120000388, 0.690001515, 5.969989203, 0.26, 0.230000547, 0.360000261, 0.219999613, 0.469998828, 0.11, 0.17, 0.09, 0.12, 0.06, 0.16, 0.3, 0.380000387, 0.12, 0.07, 0.23, 0.719999018, 0.09, 0.06, 0.559998257, 2.290005766, 0.35, 6.460002708, 2.629991793, 0.259999224, 0.07, 2.920003926, 1.809995148, 0.269999453, 0.45, 0.209999814, 0.140000435, 0.539999739, 0.439999226, 0.3, 0.15, 0.31, 0.430000982, 0.18, 0.140000273, 0.33, 4.789996466, 3.880005223, 0.65, 0.15, 0.05, 0.1, 0.04, 0.77999798, 0.04, 0.05, 0.05, 0.4, 0.890000186, 2.48000087, 0.04, 0.430000839, 1.939998436, 0.739999226, 0.15, 0.15, 0.14, 1.640000491, 0.18, 0.240000273, 1.15, 0.1, 1.020000218, 0.210000682, 0.16, 0.17, 3.779994267, 1.009999807, 5.029991793, 0.590000776, 2.230000505, 0.530000982, 0.3, 2.010006562, 0.729999128, 0.4, 0.590000435, 0.469999659, 0.29, 0.919999502, 2.850002947, 4.75, 0.639999226, 2.32, 0.379999612, 0.42, 0.509999509, 0.32, 0.330000547, 0.459998257, 2.44999163, 0.25, 1.030000341, 2.260001679, 0.12, 1.180000393, 4.769999273, 0.09, 1.260000774, 0.08, 0.04, 1.150002525, 0.18, 0.1, 0.55, 0.08, 0.09, 0.520000218, 0.1, 0.24, 0.530000498, 0.06, 0.349999273, 0.12, 0.11, 0.2, 0.7, 0.459999224, 2.090001515, 0.24, 0.07, 0.08, 5.72999477, 0.319999628, 0.240000435, 0.08, 0.18, 5.359991028, 0.1, 0.629999273, 0.209999146, 0.709999299, 2.979999479, 0.17, 0.16, 0.309999224, 0.62000202, 0.2, 0.45, 0.469999453, 0.11, 9.609996094, 1.799997824, 0.14, 3.230000839, 0.619999502, 0.229999411, 0.38, 1.12, 0.25, 0.16, 0.076748, 0.051691, 0.042645, 0.051544, 0.019803, 0.040752, 0.06183, 0.073152, 0.022944, 0.053761, 0.091904, 0.058676, 0.023826, 0.040126, 0.050901, 0.068765, 0.058565, 0.014261, 0.032102, 0.066005,
					0.06530644, 0.023031335, 0.146561047, 0.096092436, 0.06502478, 3.415327194, 1.28072497, 0.679182126, 0.275519494, 0.058426262, 0.086341586, 2.473707967, 1.291208568, 0.197009971, 0.168046158, 0.146368301, 0.125827127, 1.180176072, 6.92394539, 0.042897894, 2.547773154, 1.928381054, 0.213949898, 0.560014605, 0.578004602, 0.864617688, 0.114539362, 0.580214588, 0.039132614, 1.26906793, 1.709002394, 0.363026845, 0.322083466, 4.583452693, 0.213642835, 0.060608885, 0.205907738, 0.029217332, 0.306194977, 0.015208038, 0.284235413, 0.023875024, 0.021322965, 0.024247325, 0.025707591, 0.152674891, 0.049494471, 0.074390198, 0.010586961, 0.383509846, 0.210077561, 0.033400761, 0.048874207, 0.081325959, 2.170827545, 0.012900409, 1.670774092, 2.67407197, 0.093130848, 0.004637869, 1.957280491, 1.165876969, 0.155194943, 0.242401141, 0.103232226, 0.082086239, 0.591066469, 0.023380636, 0.572039434, 0.046079767, 0.474164931, 0.396395563, 0.191400752, 0.203954151, 0.069031445, 2.819570428, 3.505272096, 0.699837777, 0.110781139, 0.017509637, 0.19398434, 0.016810861, 0.881389471, 0.0516323, 0.024540037, 0.102740696, 0.150689339, 0.873314375, 1.72595797, 0.098459452, 0.933701119, 0.620795088, 0.130350711, 0.225691914, 0.105680626, 0.069572145, 0.575590493, 0.121314656, 0.027354382, 0.302117256, 0.067838585, 0.15052721, 0.180712097, 0.0906773, 0.10276184, 3.395875107, 0.291424058, 2.086218508, 0.496388403, 3.082524883, 0.531596087, 0.578799869, 2.00759745, 0.358354415, 0.123017672, 0.221745642, 0.886469072, 0.637433394, 0.396402546, 1.032873046, 3.279987047, 0.115087339, 1.172036188, 0.142148796, 0.761882175, 0.537251034, 0.235249866, 0.043502996, 0.262083807, 1.917454187, 0.301511758, 0.513154153, 2.162937509, 0.119606516, 0.732132049, 3.844613102, 0.024276525, 0.312473934, 0.072153678, 0.07241172, 0.731894552, 0.118514715, 0.103406554, 0.232603054, 0.080190757, 0.063931898, 0.29314594, 0.106683898, 0.25043733, 0.481853846, 0.053766077, 0.200878795, 0.04812456, 0.021305188, 0.176243737, 1.124968964, 0.189540114, 1.321924864, 0.468774829, 0.134931891, 0.069261223, 2.292162589, 0.173437273, 0.159571011, 0.355513641, 0.393151804, 3.250281971, 0.08116755, 0.341886365, 0.174187774, 0.639221889, 2.724636643, 0.072179463, 0.066248854, 0.088554073, 1.567135308, 0.066173756, 0.23163582, 0.387735993, 0.00874705, 8.543938496, 1.064655893, 0.058182013, 1.497831855, 0.528477555, 0.16919167, 0.381877235, 1.908390779, 0.120430237, 0.080536907, 0.0321172, 0.0110775, 0.0616225, 0.0162975, 0.0134929, 0.0148012, 0.0222659, 0.0477397, 0.0117794, 0.0948758, 0.149561, 0.0439519, 0.0770013, 0.101961, 0.0265993, 0.105144, 0.0430699, 0.0207374, 0.0466226, 0.059281,
					32, 2, 4, 11, 0, 864, 0, 186, 0, 0, 0, 246, 8, 49, 0, 0, 0, 0, 569, 0, 274, 78, 18, 47, 79, 0, 0, 22, 8, 232, 458, 11, 305, 550, 22, 0, 75, 0, 19, 0, 41, 0, 0, 0, 0, 21, 6, 0, 0, 27, 20, 0, 0, 26, 232, 0, 50, 408, 0, 0, 242, 215, 0, 0, 6, 4, 76, 0, 21, 0, 0, 22, 0, 0, 0, 378, 609, 59, 0, 0, 6, 5, 7, 0, 0, 0, 0, 57, 246, 0, 11, 53, 9, 33, 2, 0, 51, 0, 0, 53, 5, 43, 18, 0, 17, 342, 3, 446, 16, 347, 30, 21, 112, 20, 0, 74, 65, 47, 90, 202, 681, 0, 110, 0, 114, 0, 4, 0, 1, 360, 34, 50, 691, 8, 78, 614, 5, 16, 6, 0, 65, 0, 0, 0, 0, 0, 12, 0, 13, 0, 7, 17, 0, 0, 0, 156, 0, 530, 54, 0, 1, 1525, 16, 25, 67, 0, 682, 8, 107, 0, 14, 398, 0, 0, 10, 0, 33, 20, 5, 0, 2220, 100, 0, 832, 6, 0, 0, 237, 0, 0, 0.0692, 0.0184, 0.04, 0.0186, 0.0065, 0.0238, 0.0236, 0.0557, 0.0277, 0.0905, 0.1675, 0.0221, 0.0561, 0.0611, 0.0536, 0.0725, 0.087, 0.0293, 0.034, 0.0428,
					0.058078195, 0.03289392, 0.141364275, 0.119156855, 0.049700412, 4.658420071, 0.633255848, 0.739813857, 0.29328103, 0.077419374, 0.052454947, 2.673108089, 0.832791533, 0.131355702, 0.152595208, 0.179163888, 0.080835481, 0.812241124, 6.033788982, 0.050609064, 2.236617623, 1.46586228, 0.219967124, 0.543750757, 0.630753299, 0.91412559, 0.072395536, 0.768853295, 0.03019213, 1.522256865, 1.738679644, 0.479791112, 0.6038339, 4.518450891, 0.105414735, 0.025252656, 0.367600449, 0.012428576, 0.244934765, 0.010668856, 0.235804245, 0.008875686, 0.014004526, 0.013781055, 0.017140139, 0.109872766, 0.058180015, 0.046318944, 0.005529144, 0.299518997, 0.254452467, 0.019157619, 0.027264554, 0.111638937, 1.897974368, 0.020509508, 1.057185633, 2.53039843, 0.049007456, 0.015753762, 1.827218186, 1.379217783, 0.134187175, 0.135153663, 0.064936611, 0.06132452, 0.653363993, 0.013494034, 0.399827723, 0.026109947, 0.492340144, 0.237094366, 0.128410054, 0.145331466, 0.032834314, 2.918353208, 3.425553709, 0.65931076, 0.062762255, 0.008043958, 0.138759291, 0.012560743, 0.925810864, 0.026306325, 0.017716308, 0.068139281, 0.090353067, 0.750900541, 1.811101233, 0.097125534, 0.748424997, 0.408077053, 0.155008566, 0.080313958, 0.044609563, 0.029408411, 0.849512435, 0.048786299, 0.005914206, 0.519954375, 0.024850021, 0.270260781, 0.121234921, 0.032714699, 0.054271889, 2.771686015, 0.197379185, 2.634378514, 0.360804781, 3.283014871, 0.384800284, 0.363104466, 1.746570145, 0.297586084, 0.096272864, 0.311525131, 0.695088128, 0.458734096, 0.499349901, 1.231180819, 6.73088516, 0.056079813, 0.961285093, 0.102136221, 0.338668196, 0.274195947, 0.134802671, 0.02455829, 0.221010609, 2.453458143, 0.253366704, 0.393851704, 3.035215726, 0.053947743, 0.73460491, 3.114742907, 0.013623416, 0.370819892, 0.049019408, 0.040920528, 1.018410485, 0.12314062, 0.086028795, 0.233963371, 0.037480927, 0.028656797, 0.253243013, 0.073508963, 0.167575318, 0.330781928, 0.029433866, 0.169212029, 0.014378616, 0.014501407, 0.127519332, 1.020785491, 0.160289958, 1.967371255, 0.319105788, 0.093214721, 0.046746341, 3.907918551, 0.135319461, 0.123555332, 0.281699174, 0.316599031, 3.209083303, 0.054012183, 0.374184286, 0.091031787, 0.481044316, 2.815163085, 0.041063684, 0.051741627, 0.084589029, 1.394528044, 0.027669233, 0.227827051, 0.417148954, 0.003511008, 10.953425842, 0.958273743, 0.055461435, 2.562484895, 0.466243442, 0.054078533, 0.267109465, 1.514059674, 0.093136256, 0.06996454, 0.0437932, 0.0129578, 0.0570013, 0.016899, 0.0113305, 0.0180181, 0.0225385, 0.0470501, 0.0171837, 0.0897794, 0.155226, 0.0399135, 0.0674443, 0.088448, 0.0375282, 0.0937522, 0.063579, 0.0226713, 0.0415682, 0.0533174,
					0.039528087, 0.025400085, 0.22223502, 0.142029882, 0.088828203, 6.008326657, 0.638545114, 1.628754606, 0.567060562, 0.24304131, 0.092724227, 3.010946438, 1.443003411, 0.295205229, 0.252172496, 0.202461465, 0.027923182, 1.64404753, 10.551917849, 0.24541372, 2.495973822, 1.052275544, 0.211135339, 0.691309787, 1.133088063, 1.252025822, 0.082042453, 1.243211952, 0.066631625, 1.429437093, 2.547192176, 0.460290243, 0.235759491, 5.730432288, 0.118159164, 0.078732538, 0.081733538, 0.04402363, 0.443600397, 0.046466009, 0.27602785, 0.053149759, 0.063589602, 0.037294183, 0.101971413, 0.056354236, 0.067365714, 0.070560457, 0.024143465, 0.284975871, 0.226359102, 0.048862861, 0.02729188, 0.140577603, 1.588591037, 0.003057567, 1.427283663, 3.947911674, 0.131824102, 0.001108702, 3.979585233, 3.329630344, 0.130289331, 0.252933131, 0.081965793, 0.069714632, 0.410547061, 0.019979828, 0.427881496, 0.047234589, 0.267657471, 0.245833715, 0.291712104, 0.172669968, 0.144769632, 2.954441199, 3.979827417, 0.800287959, 0.058434718, 5.3751e-05, 0.123799323, 0.016596667, 1.498574002, 0.061292528, 0.06582573, 0.091192406, 0.166957163, 1.0003315, 2.16083819, 0.032956833, 0.63123039, 0.521438001, 0.419988041, 0.271362149, 0.064778927, 0.003598032, 0.947363601, 0.110752819, 0.010476132, 0.733232602, 0.072930098, 0.289575606, 0.427737046, 0.071284209, 0.081313817, 3.263437232, 0.24643457, 2.50950419, 0.494306172, 4.039429963, 0.452085737, 0.565283789, 2.472249816, 0.199735044, 0.182221764, 0.45360495, 3.040594995, 0.749806894, 0.618513589, 1.512976341, 4.0422603, 0.086100759, 1.769881797, 0.070259694, 0.338507365, 0.153003303, 0.238367009, 0.032149901, 0.237733268, 2.138952103, 0.268803775, 0.711144381, 3.346621035, 0.093374439, 1.114690474, 3.905596222, 0.027788876, 1.208615324, 0.095835633, 0.189869703, 2.089960072, 0.011311881, 0.184724205, 0.247942431, 0.018917835, 0.072472002, 0.419246587, 0.322553479, 0.207645878, 0.516028013, 0.06038615, 0.340305601, 0.015955802, 0.015138802, 0.224420502, 1.440857722, 0.406456046, 3.725176326, 0.751567432, 0.268495741, 0.086015251, 4.152621629, 0.171278693, 0.205722504, 0.39321817, 0.30392798, 3.901760795, 0.21839741, 0.519913887, 0.125832809, 0.777361058, 2.311962101, 0.064960081, 0.085124839, 0.13038429, 1.688525926, 0.005788447, 0.32454558, 0.536618725, 0.002501554, 9.408676331, 0.84405672, 0.051816346, 2.493990344, 0.661521677, 0.082600308, 0.384783028, 1.507320622, 0.135636288, 0.160946765, 0.0371933, 0.00998407, 0.0577903, 0.0166619, 0.0110449, 0.0163803, 0.0181385, 0.0416754, 0.0139646, 0.111614, 0.158532, 0.0331297, 0.0877754, 0.0883578, 0.0265907, 0.0979643, 0.0558081, 0.0177844, 0.0456523, 0.0539583,
					23.18, 26.95, 13.24, 17.67, 1.9, 794.38, 59.93, 103.33, 58.94, 1.9, 1.9, 220.99, 173.56, 55.28, 75.24, 9.77, 1.9, 63.05, 583.55, 1.9, 313.56, 120.71, 23.03, 53.3, 56.77, 30.71, 6.75, 28.28, 13.9, 165.23, 496.13, 113.99, 141.49, 582.4, 49.12, 1.9, 96.49, 1.9, 27.1, 4.34, 62.73, 8.34, 3.31, 5.98, 12.26, 25.46, 15.58, 15.16, 1.9, 25.65, 39.7, 1.9, 2.41, 11.49, 329.09, 8.36, 141.4, 608.7, 2.31, 1.9, 465.58, 313.86, 22.73, 127.67, 19.57, 14.88, 141.88, 1.9, 65.41, 1.9, 6.18, 47.37, 1.9, 1.9, 11.97, 517.98, 537.53, 91.37, 6.37, 4.69, 15.2, 4.98, 70.8, 19.11, 2.67, 1.9, 48.16, 84.67, 216.06, 6.44, 90.82, 54.31, 23.64, 73.31, 13.43, 31.26, 137.29, 12.83, 1.9, 60.97, 20.63, 40.1, 50.1, 18.84, 17.31, 387.86, 6.04, 494.39, 69.02, 277.05, 54.11, 54.71, 125.93, 77.46, 47.7, 73.61, 105.79, 111.16, 64.29, 169.9, 480.72, 2.08, 238.46, 28.01, 179.97, 94.93, 14.82, 11.17, 44.78, 368.43, 126.4, 136.33, 528.17, 33.85, 128.22, 597.21, 1.9, 21.95, 10.68, 19.86, 33.6, 1.9, 1.9, 10.92, 7.08, 1.9, 32.44, 24, 21.71, 7.84, 4.21, 38.58, 9.99, 6.48, 1.9, 191.36, 21.21, 254.77, 38.82, 13.12, 3.21, 670.14, 25.01, 44.15, 51.17, 39.96, 465.58, 16.21, 64.92, 38.73, 26.25, 195.06, 7.64, 1.9, 1.9, 1.9, 19, 21.14, 2.53, 1.9, 1222.94, 91.67, 1.9, 387.54, 6.35, 8.23, 1.9, 204.54, 5.37, 1.9, 0.072, 0.019, 0.039, 0.019, 0.006, 0.025, 0.024, 0.056, 0.028, 0.088, 0.169, 0.023, 0.054, 0.061, 0.054, 0.072, 0.086, 0.029, 0.033, 0.043,
					0.061426075, 0.031831916, 0.149057853, 0.150297046, 0.063019188, 8.577252478, 0.258187878, 1.170363549, 0.306981436, 0.124361236, 0.027354705, 3.133756862, 0.416576378, 0.084904707, 0.140528069, 0.187890775, 0.088305718, 0.339883214, 6.634089838, 0.030812588, 2.06066417, 1.113145063, 0.288931657, 0.573764727, 0.961896424, 0.910743527, 0.041354439, 1.151834508, 0.027049022, 2.1220079, 2.472501269, 0.805468382, 1.447715191, 4.324671886, 0.081983093, 0.01280823, 0.419106653, 0.001696908, 0.12218513, 0.00267315, 0.074268458, 0.004845909, 0.001457333, 0.002980311, 0.021069514, 0.084964207, 0.078245354, 0.00518453, 0.002808236, 0.182360515, 0.253147623, 0.009125664, 0.00689027, 0.147175841, 1.41388806, 0.020710956, 0.431137268, 2.487731335, 0.02045269, 0.018885125, 2.031041265, 2.113012391, 0.110065675, 0.182068253, 0.006852808, 0.021134276, 0.77609134, 0.002391939, 0.038343725, 0.006156273, 0.081020095, 0.124183971, 0.059893588, 0.034003939, 0.025326733, 2.73580491, 3.155941253, 0.43579328, 0.069418424, 0.00232623, 0.011832201, 0.005370941, 1.071104374, 0.010980483, 0.001569405, 0.015417345, 0.183163356, 0.516552726, 2.390574145, 0.010488183, 0.147101242, 0.285706397, 0.221106403, 0.039069406, 0.032172065, 0.023520451, 0.915068649, 0.023665106, 0.001588822, 0.673949441, 0.01898617, 0.357953716, 0.152366262, 0.03130126, 0.068701823, 2.539091702, 0.127840652, 3.777847616, 0.376137787, 3.378191969, 0.256671584, 0.104799233, 1.238301314, 0.450133322, 0.042952526, 0.435065499, 0.193583353, 0.148305649, 0.869391551, 1.922238895, 5.457787758, 0.028428747, 0.866028161, 0.088772981, 0.162800805, 0.118976512, 0.097890386, 0.007678823, 0.151157734, 2.411521891, 0.181367201, 0.408372307, 3.695874029, 0.064760754, 0.576413595, 2.908579763, 0.008350305, 0.544282004, 0.007528026, 0.024653916, 1.764289698, 0.122739514, 0.062084664, 0.253335423, 0.033143637, 0.001928621, 0.188124443, 0.031311764, 0.059361129, 0.095489234, 0.019550695, 0.141663396, 0.003372908, 0.013142627, 0.145341785, 0.855444034, 0.174581453, 5.129621829, 0.24379326, 0.039546237, 0.020949593, 9.019615949, 0.047323608, 0.106687168, 0.071183277, 0.079180592, 4.071651775, 0.072344705, 0.53309475, 0.086650158, 0.29464987, 3.176975964, 0.034013497, 0.013012894, 0.123720116, 0.377893834, 0.009298536, 0.215600533, 0.444794286, 0.002695769, 13.419120278, 0.867600218, 0.02631879, 4.378921288, 0.344319078, 0.025284157, 0.058827343, 1.341486679, 0.045879219, 0.027607483, 0.0706288, 0.0138991, 0.0455021, 0.014849, 0.00674192, 0.026439, 0.0214826, 0.0440199, 0.0241895, 0.0908219, 0.172674, 0.0273258, 0.0563431, 0.0495703, 0.0542482, 0.0746629, 0.109035, 0.0254891, 0.0264554, 0.0456223,
					3.3, 1.7, 33.6, 16.1, 3.2, 617, 272.5, 61.1, 94.6, 9.5, 7.3, 231, 190.3, 19.3, 49.1, 17.1, 6.4, 174, 883.6, 3.4, 349.4, 289.3, 7.2, 99.3, 26, 82.4, 8.9, 43.1, 2.3, 61.7, 228.9, 55.6, 37.5, 421.8, 14.9, 7.4, 33.2, 0.2, 24.3, 1.5, 48.8, 0.2, 7.3, 3.4, 1.6, 15.6, 4.1, 7.9, 0.5, 59.7, 23, 1, 3.5, 6.6, 425.2, 0.2, 292.3, 413.4, 0.2, 0.2, 334, 163.2, 10.1, 23.9, 8.4, 6.7, 136.5, 3.8, 73.7, 0.2, 264.8, 83.9, 0.2, 52.2, 7.1, 449.7, 636.3, 83, 26.5, 0.2, 12.9, 2, 167.8, 9.5, 0.2, 5.8, 13.1, 90.3, 234.2, 16.3, 215.6, 61.8, 7.5, 22.6, 0.2, 8.1, 52.2, 20.6, 1.3, 15.6, 2.6, 11.4, 24.3, 5.4, 10.5, 644.9, 11.8, 420.2, 51.4, 656.3, 96.4, 38.4, 257.1, 23.1, 7.2, 15.2, 144.9, 95.3, 32.2, 79.7, 378.1, 3.2, 184.6, 2.3, 199, 39.4, 34.5, 5.2, 19.4, 222.3, 50, 75.5, 305.1, 19.3, 56.9, 666.3, 3.1, 16.9, 6.4, 0.2, 36.1, 6.1, 3.5, 12.3, 4.5, 9.7, 27.2, 6.6, 48.7, 58.2, 1.3, 10.3, 3.6, 2.1, 13.8, 141.6, 13.9, 76.7, 52.3, 10, 4.3, 266.5, 13.1, 5.7, 45, 41.4, 590.5, 4.2, 29.7, 29, 79.8, 321.9, 5.1, 7.1, 3.7, 243.8, 9, 16.3, 23.7, 0.3, 1710.6, 126.1, 11.1, 279.6, 59.6, 17.9, 49.5, 396.4, 13.7, 15.6, 0.06888, 0.021037, 0.03039, 0.020696, 0.009966, 0.018623, 0.024989, 0.071968, 0.026814, 0.085072, 0.156717, 0.019276, 0.050652, 0.081712, 0.044803, 0.080535, 0.056386, 0.027998, 0.037404, 0.066083,
					0.674995699, 0.589645178, 1.189067034, 0.462499504, 0.605460903, 3.573373315, 1.065445546, 0.31444833, 0.589852457, 0.246951424, 1.111766964, 2.967840934, 2.299755865, 1.686058219, 0.245163782, 1.046334652, 1.201770702, 1.277836748, 4.399995525, 0.091071867, 4.15967899, 1.587964372, 0.523770553, 1.374854049, 0.734992057, 0.31706632, 0.596789898, 0.463812837, 0.580830874, 1.457127446, 2.283037894, 0.839348444, 0.411543728, 1.812173605, 0.877842609, 0.476331437, 0.464590585, 0.35964586, 0.426069419, 0.266775558, 0.417547309, 0.315256838, 0.30421529, 0.180198883, 0.285186418, 0.804404505, 0.520701585, 0.41009447, 0.269124919, 0.450795211, 0.625792937, 0.32078471, 0.259854426, 0.363981358, 4.162454693, 0.831998835, 4.956476453, 2.037575629, 1.114178954, 0.274163536, 3.521346591, 2.415974716, 0.581001076, 0.985885486, 0.374784947, 0.498011337, 1.546725076, 0.81346254, 0.737846301, 0.341932741, 0.618614612, 2.067388546, 0.531773639, 0.465349326, 0.380925433, 3.65807012, 5.002338375, 0.661095832, 0.546169219, 0.303437244, 0.425193716, 0.219005213, 0.669206193, 0.406042546, 0.224154698, 0.35402891, 0.576231691, 1.495264661, 2.392638293, 0.269496317, 2.306919847, 1.241586045, 0.65577338, 0.711495595, 0.775624818, 0.198679914, 0.850116543, 0.794584081, 0.588254139, 0.456058589, 0.366232942, 0.430073179, 1.036079005, 0.337502282, 0.481144863, 3.452308792, 0.910144334, 2.572577221, 1.440896785, 0.99870098, 1.348272505, 1.205509425, 1.402122097, 0.799966711, 0.530641901, 0.402471997, 1.234648153, 0.945453716, 0.613230817, 1.217683028, 1.751412803, 0.89517149, 1.823161023, 0.994227284, 0.847312432, 1.320626678, 0.949599791, 0.542185658, 0.83039281, 1.114132523, 0.779827336, 1.290709079, 1.551488041, 0.718895136, 0.780913179, 4.448982584, 0.35011051, 0.618778365, 0.422407388, 0.362495245, 0.445669347, 0.72038474, 0.261258229, 0.37874827, 0.72436751, 0.516260502, 0.794797115, 0.43340962, 0.768395107, 3.29519344, 0.499869138, 0.496334956, 0.38372361, 0.573154753, 0.628599063, 0.720013799, 0.436220437, 0.55626163, 0.728970584, 0.50720003, 0.284727562, 2.210952064, 0.570562395, 0.811019594, 0.664884513, 0.93253606, 5.894735673, 0.433748126, 0.593795813, 0.523549536, 2.996248013, 2.063050067, 0.388680158, 0.474418852, 0.275658381, 0.998911631, 0.634408285, 0.527640634, 0.314700907, 0.305792277, 8.002789424, 2.113077156, 0.526184203, 1.737356217, 0.983844803, 0.551333603, 0.507506011, 1.89965079, 0.429570747, 0.716795463, 0.075559, 0.053764, 0.037684, 0.044693, 0.028483, 0.033893, 0.053472, 0.078036, 0.03004, 0.059869, 0.095793, 0.051997, 0.021908, 0.045001, 0.042029, 0.068212, 0.056413, 0.015725, 0.035974, 0.071456,
					0.086772353, 0.041489234, 0.145522693, 0.684872641, 0.038182359, 4.588509426, 0.214514127, 2.310001798, 0.094392928, 0.085866464, 0.10407131, 4.326376802, 0.119230449, 0.051178971, 0.037945162, 0.768984102, 0.089892499, 0.080461407, 5.114292481, 0.028175939, 1.988992665, 2.215762795, 1.93785533, 0.295994975, 2.04182192, 1.073717504, 0.105321639, 1.420133178, 0.080617856, 6.927918567, 3.653639043, 1.071049041, 2.585629399, 6.123585849, 0.058085432, 0.108840477, 0.046084289, 0.203642191, 0.46507458, 0.025094905, 0.142117407, 0.031239309, 0.016280777, 0.029795942, 0.047732049, 0.116389136, 0.336043365, 0.021519291, 0.026444908, 0.269064319, 0.527627637, 0.03852806, 0.04732831, 0.607702969, 1.93576598, 0.09363087, 5.781174762, 2.258646191, 0.026415889, 0.011551984, 1.637753572, 1.761938656, 0.087005765, 0.004280896, 0.148647886, 0.044628305, 0.235654912, 0.732190858, 0.04898364, 0.016856105, 0.129334019, 0.168730118, 0.077835842, 0.072799197, 0.03542684, 7.398376755, 3.334268838, 0.708348805, 0.092004559, 0.061734224, 0.047373731, 0.043468817, 1.575630749, 0.045946765, 0.043087788, 0.04957537, 0.141447452, 0.962852023, 2.950112395, 0.040069839, 0.134340383, 2.022526978, 0.407828168, 0.021912937, 0.057096272, 0.124231767, 1.446258879, 0.060414502, 0.063511937, 1.485608296, 0.042135216, 1.450911972, 0.054326874, 0.048270788, 0.072731066, 2.352357809, 0.961653423, 5.657997764, 0.161901882, 2.848589688, 0.086837961, 0.05617398, 2.55862178, 0.248112786, 0.287491221, 0.653587912, 0.105167714, 0.055804307, 0.936346554, 2.818055661, 9.18921051, 0.711833578, 1.856292472, 0.041327738, 0.117750493, 0.075637987, 0.075603922, 0.08545894, 0.064922501, 3.67091655, 0.086159175, 0.827391264, 6.908142924, 0.086710531, 1.546051817, 3.320463486, 0.067881151, 2.240359335, 0.02907552, 0.024018688, 2.572562133, 0.667514386, 0.057387721, 0.561105174, 0.03249847, 0.045679288, 0.68715629, 0.024744157, 0.188378332, 0.454841057, 0.05477351, 0.269153898, 0.048484013, 0.06242815, 0.069988164, 0.632160724, 0.452258389, 5.266555093, 0.049244809, 0.032418984, 0.035199737, 6.764718412, 0.083433936, 0.094733583, 0.036104364, 0.04382335, 3.894213263, 0.056354401, 0.461507798, 0.048169854, 0.308629484, 5.126996382, 0.063299974, 0.037743293, 0.240637305, 0.199254052, 0.052218599, 0.307081397, 0.848692284, 0.058399592, 13.397205976, 2.035426441, 0.054281453, 7.643787012, 0.780383422, 0.08101024, 0.08262015, 0.391832532, 0.082144494, 0.059610809, 0.066363, 0.054021, 0.037784, 0.047511, 0.022651, 0.048841, 0.071571, 0.058368, 0.025403, 0.045108, 0.100181, 0.061361, 0.021069, 0.03823, 0.053861, 0.089298, 0.053536, 0.012313, 0.027173, 0.065359,
					0.245103884, 0.396680459, 0.60259659, 0.388851211, 0.169418318, 4.002207051, 1.687379028, 0.888477494, 0.707887462, 0.137794758, 0.813368876, 2.004043618, 1.229351544, 0.369955633, 0.2306078, 0.797312266, 0.224720163, 0.413998478, 7.046040747, 0.061906941, 2.610915158, 2.83518824, 0.639959824, 1.963743484, 1.231717882, 1.043646688, 0.393592403, 0.553182381, 0.362258381, 2.615490291, 4.57326098, 0.797011018, 0.896230099, 5.6751528, 0.330930746, 0.37812445, 0.228513969, 0.141228769, 0.208490043, 0.05317395, 0.46042895, 0.116150694, 0.064283924, 0.068021312, 0.140874297, 0.2869539, 0.240779967, 0.09866768, 0.058795467, 0.493466139, 0.643789822, 0.079309054, 0.089730309, 0.50099201, 3.38156008, 0.325994754, 5.655671758, 1.541274942, 0.211396072, 0.082097991, 2.054032646, 0.761883191, 0.341199781, 0.396451596, 0.14488845, 0.103298167, 0.854941063, 0.469236454, 0.344916944, 0.082042638, 0.564344515, 1.049351234, 0.169182004, 0.214023208, 0.337275045, 4.60765751, 5.779760276, 0.501772274, 0.202912169, 0.083503104, 0.104497835, 0.062211382, 1.265119124, 0.061813267, 0.03737175, 0.130615907, 0.693502085, 0.880367224, 2.142804923, 0.038837539, 1.068089269, 2.098306482, 0.448091631, 0.295630473, 0.393309252, 0.214758763, 1.22694582, 0.411528885, 0.517165713, 0.865174971, 0.16490386, 0.371130819, 0.335358981, 0.264037782, 0.131445137, 4.447741802, 0.609758408, 3.798086697, 0.731000508, 3.467889377, 0.737657764, 0.379154866, 2.320367575, 0.697254375, 0.107585891, 0.228053049, 0.417871055, 0.387566384, 0.273240209, 1.842353422, 3.173624632, 0.474607818, 2.267368199, 0.310254294, 1.12563823, 0.924249331, 0.410818877, 0.401697351, 0.49076449, 1.310627756, 0.287148699, 0.827466845, 2.36890891, 0.164073566, 1.031968282, 5.706447436, 0.113618297, 0.421211392, 0.052343656, 0.053998922, 1.063048958, 0.118187043, 0.058461221, 0.275756639, 0.303820579, 0.088741408, 0.404916338, 0.067974475, 0.484518123, 1.412153565, 0.13174532, 0.257205952, 0.121214423, 0.159463305, 0.209585392, 0.485673084, 0.153134973, 1.990179836, 0.174718362, 0.083149696, 0.137582926, 6.378949466, 0.188267059, 0.280539346, 0.067140987, 0.264484864, 10.165067155, 0.155778076, 0.294006077, 0.191009158, 1.616501717, 2.466363084, 0.145305726, 0.149535968, 0.109375281, 1.295115095, 0.258250207, 0.287452076, 0.247675666, 0.17741256, 10.861831869, 1.518647299, 0.149072919, 2.157390423, 0.496751123, 0.367029379, 0.180581516, 1.639361417, 0.103038434, 0.189789265, 0.063003, 0.049585, 0.04755, 0.048622, 0.015291, 0.044058, 0.072012, 0.03781, 0.022358, 0.066563, 0.107325, 0.080621, 0.023976, 0.041578, 0.028532, 0.081767, 0.055167, 0.009698, 0.032219, 0.072265,
					0.42405754, 0.271250376, 0.764825201, 0.401980425, 0.140794893, 5.062142812, 2.32745344, 0.519625401, 0.507639701, 0.071616388, 1.072849735, 3.123615724, 1.879768813, 0.664578644, 0.104402693, 1.069798507, 0.402468665, 0.59533039, 5.60103324, 0.000110712, 4.776313857, 1.793212924, 0.339928462, 1.257892965, 0.741295542, 0.488516386, 0.272681884, 0.332128789, 0.392890955, 2.558126695, 4.639398318, 0.987385467, 0.677157762, 5.468173089, 0.49825795, 0.296507145, 0.13801491, 0.132223556, 0.192515195, 0.013803595, 0.328676783, 0.094404276, 0.045491635, 0.009271038, 0.121398119, 0.367358512, 0.278306063, 0.065286972, 0.014185552, 0.536563257, 0.592378803, 0.071724886, 0.041522603, 0.366698667, 3.930811583, 0.563978915, 6.845986165, 2.226004549, 0.353825055, 0.002829803, 3.797608224, 2.027714606, 0.282308308, 0.80198353, 0.169487062, 0.141910872, 1.128143826, 0.489462975, 0.408708404, 0.02380644, 0.904437688, 1.839460723, 0.19667133, 0.125388186, 0.497496251, 4.359960777, 6.201341249, 0.706300529, 0.236258684, 0.055016203, 0.092864269, 0.016998748, 0.957089079, 0.037738461, 0.018574182, 0.075225605, 0.703446369, 1.021152682, 2.22610419, 0.019952549, 1.726515479, 1.158823278, 0.344518589, 0.176377788, 0.403179438, 0.080797749, 0.708201457, 0.456455246, 0.180072253, 0.538060086, 0.083244488, 0.235273346, 0.427986716, 0.111914514, 0.092858734, 4.667976453, 0.835065436, 3.955582326, 1.24115439, 2.628965883, 1.346684133, 0.67039214, 1.526697708, 1.06796622, 0.054379608, 0.183116838, 0.828615345, 0.354736216, 0.335899646, 1.341323449, 2.159899413, 0.617964388, 2.022264246, 0.477448489, 1.129405944, 1.248400502, 0.659780378, 0.124572237, 0.64913983, 1.017998492, 0.303426657, 1.225974643, 2.088967227, 0.169151604, 0.60596983, 6.645417766, 0.15186279, 0.47705214, 0.039141184, 0.025562335, 0.529678067, 0.221132079, 0.066497056, 0.187213189, 0.555296865, 0.099536892, 0.479069315, 0.048841785, 0.606743708, 1.854514258, 0.076842003, 0.202206939, 0.123868108, 0.192903795, 0.274767702, 0.526873728, 0.124040819, 1.020373268, 0.276139426, 0.11234297, 0.044355728, 5.074607897, 0.213437583, 0.276570096, 0.143396629, 0.480201901, 6.641252774, 0.081147599, 0.358925565, 0.235545698, 2.37326614, 2.511548974, 0.18066567, 0.088670492, 0.042435978, 1.853035143, 0.252220059, 0.270634816, 0.078641075, 0.148865811, 11.081231885, 1.642242475, 0.209350089, 2.027034834, 0.628362476, 0.310510022, 0.091930966, 2.286257438, 0.164177306, 0.267739693, 0.080009, 0.052947, 0.041171, 0.050146, 0.015018, 0.035929, 0.061392, 0.064793, 0.021709, 0.063895, 0.106292, 0.057047, 0.02344, 0.047712, 0.039604, 0.06298, 0.052863, 0.014987, 0.037434, 0.070634,
					0.164520503, 0.13378666, 0.301753825, 0.667759687, 0.102671327, 4.669240542, 0.375561674, 1.87855413, 0.218506732, 0.132842206, 0.172392722, 3.619009961, 0.250359878, 0.124852096, 0.095864103, 0.781197529, 0.20835749, 0.202345038, 5.379197337, 0.056958677, 1.746841825, 2.0763385, 1.797395963, 0.6794945, 1.991708503, 1.179185026, 0.168866582, 1.45353378, 0.161623558, 4.710380602, 3.940278094, 0.956160014, 1.79641177, 5.600293701, 0.164776154, 0.227985715, 0.203510704, 0.217922574, 0.520603634, 0.071028795, 0.233803448, 0.064010956, 0.054619397, 0.086763961, 0.094641478, 0.17892575, 0.403058144, 0.048270492, 0.042419644, 0.297464704, 0.583097899, 0.070620812, 0.086464596, 0.627289885, 2.383373587, 0.150019883, 5.862144662, 2.011938106, 0.10985077, 0.049580542, 1.705745478, 1.499396112, 0.184898463, 0.108290367, 0.192452773, 0.074328424, 0.426954296, 0.737495653, 0.124928924, 0.058185276, 0.192839562, 0.273324822, 0.142904022, 0.172462927, 0.119658259, 7.422424301, 4.121859707, 0.760101625, 0.156046906, 0.098468307, 0.074017139, 0.075986849, 1.649041128, 0.061484905, 0.05991523, 0.098183513, 0.275397853, 1.183422487, 2.957377967, 0.060573566, 0.338235196, 1.899027463, 0.525891518, 0.068009985, 0.123677157, 0.200110404, 1.590125986, 0.137838674, 0.136748511, 1.445166153, 0.107833373, 1.381061149, 0.129943936, 0.137065096, 0.146012906, 2.359299936, 0.843914331, 5.427387027, 0.263292147, 2.507254736, 0.189228647, 0.139146074, 2.427076788, 0.376188219, 0.330866332, 0.498124825, 0.193144225, 0.150312624, 0.827642694, 2.786459548, 7.524532664, 0.539886133, 1.955763604, 0.13270577, 0.233037817, 0.165733855, 0.155868082, 0.198123473, 0.145627442, 3.528654954, 0.153102275, 0.889109054, 4.911310946, 0.119728464, 1.52998822, 3.075670543, 0.101039395, 2.131494369, 0.058369398, 0.043905868, 2.144674342, 0.5214196, 0.088844939, 0.579178348, 0.134177423, 0.08868466, 0.683790243, 0.064089109, 0.27624561, 0.603795093, 0.100007515, 0.264443242, 0.07963088, 0.091694198, 0.139818982, 0.655149301, 0.420781565, 5.68528532, 0.126113134, 0.070592995, 0.083304052, 7.28708651, 0.172110578, 0.17628843, 0.080506454, 0.165573576, 5.877878503, 0.112513257, 0.458878972, 0.102291161, 0.645002709, 4.300445029, 0.107177686, 0.082093349, 0.228482448, 0.278600785, 0.094238794, 0.325642029, 0.868890049, 0.081513165, 12.950251994, 1.879650916, 0.09739139, 7.15849766, 0.733119107, 0.145626117, 0.152049202, 0.573328824, 0.113423271, 0.125003102, 0.067997, 0.055503, 0.036288, 0.046867, 0.021435, 0.050281, 0.068935, 0.055323, 0.02641, 0.041953, 0.101191, 0.060037, 0.019662, 0.036237, 0.055146, 0.096864, 0.057136, 0.011785, 0.02473, 0.066223,
					0.531344742, 0.266631781, 0.610524242, 0.479415354, 0.145193836, 4.395589145, 2.490407258, 0.797726764, 0.617331366, 0.086320436, 1.058818226, 2.850794598, 1.685541958, 0.623180282, 0.163023963, 1.178483844, 0.358512949, 0.572153867, 4.775857514, 0.004224284, 4.045465925, 1.897932882, 0.427923043, 1.417171473, 0.993358642, 0.723368327, 0.349584077, 0.412573692, 0.453464468, 2.765967911, 4.395995003, 0.944775779, 1.220289602, 4.992256584, 0.444851397, 0.432227777, 0.176454495, 0.103023046, 0.192557924, 0.012280201, 0.599859067, 0.090487083, 0.045755066, 0.025135568, 0.108027826, 0.41943365, 0.307712278, 0.070917051, 0.019106538, 0.827369996, 0.609556427, 0.066812844, 0.070095729, 0.420152907, 4.31603981, 0.501174376, 5.070603955, 2.126974783, 0.311739636, 0.042113153, 3.211891588, 1.628511729, 0.323774881, 0.779002069, 0.188076678, 0.128912693, 1.17507728, 0.565654138, 0.405987508, 0.044371788, 1.330027314, 1.704580053, 0.21768989, 0.196370346, 0.488138895, 5.05239799, 6.674964742, 0.790881216, 0.266730243, 0.050284344, 0.098902029, 0.02328159, 1.570975979, 0.044860498, 0.016021778, 0.116848629, 0.75484032, 1.333037626, 2.539936322, 0.022355802, 1.944915128, 1.371519672, 0.396818483, 0.23086586, 0.538351193, 0.117103191, 0.764761023, 0.532587532, 0.323334201, 0.635697033, 0.101194285, 0.285369684, 0.446883087, 0.146079999, 0.106999973, 3.745215188, 0.711941869, 3.559567653, 1.18899172, 3.084090802, 1.178336151, 0.587105765, 1.645856748, 1.047872072, 0.078854093, 0.182164122, 0.756113129, 0.449949835, 0.319772739, 1.499944856, 2.134156546, 0.664457704, 1.923647217, 0.498902533, 1.382576296, 1.172710577, 0.733559201, 0.232722239, 0.64897747, 1.067438502, 0.303042511, 1.169845557, 2.176841262, 0.185666747, 0.667935113, 5.538833003, 0.198118657, 0.503943335, 0.048383536, 0.030638664, 0.964172901, 0.185024339, 0.04402957, 0.1999174, 0.592439554, 0.143667666, 0.614180565, 0.040725071, 0.765641182, 2.174974075, 0.133590865, 0.217347672, 0.125958817, 0.208747809, 0.263834003, 0.570488409, 0.134714779, 1.863419357, 0.264729772, 0.10044741, 0.074554161, 5.545635324, 0.271724216, 0.338670344, 0.138599247, 0.65118087, 7.474120415, 0.108442089, 0.374198514, 0.267599595, 2.60441128, 2.688525915, 0.201107356, 0.119873351, 0.052396485, 2.371412865, 0.282178057, 0.297627071, 0.134209258, 0.22973234, 11.786948184, 2.030164484, 0.222793132, 2.397008325, 0.758096789, 0.362295352, 0.127446562, 2.284500453, 0.201953893, 0.33768812, 0.085788, 0.057731, 0.042028, 0.056462, 0.010447, 0.039548, 0.067799, 0.064861, 0.02104, 0.055398, 0.100413, 0.059401, 0.019898, 0.042789, 0.039579, 0.069262, 0.055498, 0.01443, 0.033233, 0.064396,
					0.061995451, 0.071787018, 0.324146307, 0.48272325, 0.017012888, 5.640569094, 0.523802485, 2.773840824, 0.412259505, 0.072474815, 0.26669947, 3.319795598, 0.641219533, 0.10139171, 0.105145707, 0.848857861, 0.02737889, 0.203417806, 6.983296655, 0.004325163, 3.676487445, 2.099517664, 0.923332845, 0.989386407, 1.214523767, 1.225782488, 0.189953993, 0.893983216, 0.086540313, 4.386107333, 5.633764157, 1.133417157, 3.025390235, 8.135368256, 0.117555479, 0.12776016, 0.047711951, 0.148339563, 0.470979475, 0.007582115, 0.267133294, 0.022845535, 0.010786746, 0.011450563, 0.074259162, 0.132461376, 0.242868565, 0.008688476, 0.007966889, 0.387974066, 0.713657288, 0.029290769, 0.013401682, 0.666044528, 3.178490964, 0.135231094, 6.796751125, 3.253593795, 0.018324171, 0.000109001, 2.561068546, 1.166641775, 0.150166421, 0.040648681, 0.143807298, 0.02226565, 0.487603228, 0.622545469, 0.164436842, 0.009922368, 0.133016192, 0.383963916, 0.121953673, 0.072782198, 0.060192574, 6.137883228, 5.109951861, 0.997263914, 0.101998846, 0.018497482, 0.03348949, 0.011945428, 2.787245776, 0.000109001, 0.003428084, 0.016479873, 0.304517393, 1.260751123, 3.06753877, 0.004231422, 0.367280212, 1.679096158, 0.423800539, 0.020342871, 0.080963818, 0.083580934, 1.394642593, 0.104298769, 0.064165663, 1.15090747, 0.036382379, 0.55976194, 0.059246444, 0.049078824, 0.049489758, 3.724953683, 0.724887668, 4.517199656, 0.281847349, 5.612663729, 0.304709235, 0.112097795, 1.591990129, 0.447466853, 0.164991657, 0.54964555, 0.209386705, 0.171531722, 0.574024731, 2.359864552, 4.753201153, 0.655493224, 2.669910479, 0.105335369, 0.389742063, 0.300884387, 0.195263436, 0.123875362, 0.116237656, 2.552275429, 0.089626134, 1.172764365, 4.055778486, 0.062695238, 1.066677979, 5.735630021, 0.016624844, 0.808384672, 0.005560145, 0.000347713, 2.878143952, 0.17646184, 0.013106289, 0.196264065, 0.06688306, 0.000109001, 0.582218341, 0.002231252, 0.229612944, 0.597613653, 0.009791567, 0.216845648, 0.000109001, 0.033775073, 0.054821001, 0.7534187, 0.255631501, 4.340476213, 0.036965535, 0.020444242, 0.012699715, 7.284279143, 0.053844351, 0.084699285, 0.016433002, 0.066804579, 8.185198097, 0.032534641, 0.405822992, 0.051593479, 0.524555683, 3.226243245, 0.041421498, 0.061468976, 0.1470043, 0.560195764, 0.088095759, 0.303691165, 0.297902118, 0.066305354, 10.911840366, 1.918251057, 0.050145945, 2.585318015, 0.627525728, 0.11460809, 0.058704709, 0.972725592, 0.026833885, 0.133140453, 0.074923, 0.0505, 0.038734, 0.053195, 0.0113, 0.037499, 0.068513, 0.059627, 0.021204, 0.058991, 0.102504, 0.067306, 0.022371, 0.043798, 0.037039, 0.084451, 0.04785, 0.012322, 0.030777, 0.077097,
					0.289760345, 0.342709634, 0.718300668, 0.367886518, 0.0725626, 4.19952265, 3.691718604, 0.710404342, 0.490823952, 0.052871678, 1.043523521, 2.008049414, 1.553982959, 0.449223552, 0.070963821, 1.093168815, 0.139582694, 0.508273036, 5.624270628, 0.003118838, 3.707570018, 2.882498207, 0.540565704, 1.98785624, 1.230079486, 1.182218128, 0.398459557, 0.542851891, 0.41737659, 2.637569081, 4.921959964, 0.793976563, 0.971187101, 5.448912766, 0.394563504, 0.492571299, 0.182603367, 0.132608825, 0.173586122, 0.017141968, 0.595831116, 0.083860488, 0.036837325, 0.027627162, 0.116377119, 0.334353188, 0.205566131, 0.06882953, 0.021633405, 0.79620177, 0.465713591, 0.062744857, 0.043967523, 0.378253079, 3.366451402, 0.48377469, 6.326126715, 1.772882594, 0.19303756, 0.025984034, 2.595024049, 0.980236499, 0.351568331, 0.667015573, 0.10815372, 0.090350953, 0.980143366, 0.618539982, 0.354791389, 0.052464777, 0.984440777, 1.12234819, 0.146156313, 0.190256883, 0.392801743, 4.449605122, 5.986085778, 0.578286681, 0.224857913, 0.047417869, 0.073768891, 0.02029185, 1.625802715, 0.031817027, 0.015954525, 0.094986525, 0.713520963, 0.952803363, 2.116139937, 0.012660509, 1.303449237, 1.902789461, 0.327908181, 0.288401051, 0.418241394, 0.119283083, 0.994116603, 0.453784839, 0.361756388, 0.717288404, 0.097506652, 0.200252021, 0.383751236, 0.146206205, 0.091990755, 4.877234074, 0.595474107, 3.112134717, 0.924744908, 3.584581583, 1.23959233, 0.512666905, 1.86305726, 0.955917766, 0.058056067, 0.241269232, 0.551359124, 0.486660696, 0.253426381, 1.572179322, 2.789256813, 0.577913042, 2.17234998, 0.379983795, 1.746675706, 1.0041827, 0.576198956, 0.295694231, 0.636427012, 1.138771705, 0.291850288, 0.980917255, 1.905694316, 0.158541858, 0.735801862, 5.167267165, 0.082544434, 0.297203203, 0.023885222, 0.009802221, 1.098149199, 0.061054054, 0.02127307, 0.13204892, 0.320030703, 0.103698871, 0.41152475, 0.030961092, 0.442192031, 1.474552928, 0.057108109, 0.136000409, 0.072755518, 0.17650539, 0.179523334, 0.380460546, 0.103068007, 1.987376163, 0.160373469, 0.057915259, 0.060539607, 6.419218432, 0.171959625, 0.282755433, 0.046896769, 0.361827346, 7.660749919, 0.08548255, 0.275364684, 0.182111094, 1.812622556, 2.936397672, 0.163865947, 0.135944973, 0.08035138, 2.828941339, 0.271020707, 0.294499028, 0.133659895, 0.255089465, 11.291813609, 1.513322757, 0.167243117, 2.085093473, 0.620211936, 0.348307577, 0.145275987, 2.18899524, 0.093529663, 0.240182683, 0.059954, 0.042032, 0.052518, 0.054641, 0.008189, 0.040467, 0.070691, 0.039935, 0.018393, 0.069555, 0.109563, 0.081967, 0.018694, 0.046979, 0.031382, 0.091102, 0.055887, 0.010241, 0.033496, 0.064313,
					34, 51, 35, 10, 30, 384, 439, 92, 128, 1, 32, 221, 236, 78, 70, 81, 10, 79, 542, 1, 372, 135, 41, 94, 61, 48, 18, 70, 30, 90, 320, 91, 124, 387, 34, 68, 1, 24, 35, 1, 104, 33, 1, 1, 34, 45, 18, 15, 5, 110, 54, 21, 3, 51, 385, 38, 593, 123, 20, 16, 309, 141, 30, 76, 34, 23, 235, 57, 1, 1, 156, 158, 1, 37, 116, 375, 581, 134, 1, 7, 49, 1, 70, 1, 1, 7, 141, 64, 179, 14, 247, 97, 24, 33, 55, 1, 68, 52, 17, 44, 10, 22, 43, 1, 11, 460, 102, 294, 136, 75, 225, 95, 152, 183, 4, 24, 77, 1, 20, 134, 258, 64, 148, 55, 117, 146, 82, 7, 49, 72, 25, 110, 131, 69, 62, 671, 5, 13, 16, 1, 55, 10, 17, 23, 48, 39, 47, 6, 111, 182, 9, 14, 1, 55, 47, 28, 1, 131, 45, 1, 21, 307, 26, 64, 1, 74, 1017, 14, 31, 34, 176, 197, 29, 21, 6, 295, 36, 35, 3, 1, 1048, 112, 19, 236, 92, 25, 39, 196, 26, 59, 0.0646, 0.0453, 0.0376, 0.0422, 0.0114, 0.0606, 0.0607, 0.0639, 0.0273, 0.0679, 0.1018, 0.0751, 0.015, 0.0287, 0.0681, 0.0488, 0.0622, 0.0251, 0.0318, 0.0619,
					0.1159435373, 0.2458816714, 0.1355713516, 0.9578712472, 0.0775041665, 8.4408676914, 0.2327281954, 9.137947033, 0.1137687264, 0.0582110367, 0.3309250853, 5.2854173238, 0.1727184754, 0.8191776581, 0.0009722083, 0.6946680829, 0.0966719296, 0.2990806606, 7.3729791633, 0.0005604799, 3.5773486727, 2.8076062202, 3.0815651393, 0.5575702616, 2.2627839242, 1.1721237455, 0.0482085663, 3.3184632572, 0.2275494971, 2.8251848421, 9.522860803, 2.3191131858, 0.0483235836, 4.413871527, 0.0343694246, 0.094838346, 0.0627691644, 0.5712158076, 0.2238609194, 0.0205779319, 0.1527276944, 0.0206129952, 0.0328079744, 0.1239000315, 0.0802374651, 0.030581884, 0.1930408758, 0.054096725, 0.0018843293, 0.2406073246, 0.329945462, 0.0373753435, 0.000591894, 0.119290461, 1.3184058362, 0.2231434272, 6.0541970908, 4.3977466558, 0.1347413792, 0.0001480536, 5.2864094506, 6.8883522181, 0.5345755286, 0.3991624551, 0.2107928508, 0.1055933141, 0.1874527991, 0.2427875732, 0.0433577842, 2.2173e-06, 0.0927357503, 0.01092383, 0.0663619185, 0.0128777966, 0.0722334577, 4.3016010974, 1.1493262595, 0.4773694701, 0.0458112245, 0.031003075, 0.023349397, 8.0023e-06, 0.8419347601, 0.0027817812, 0.0361207581, 0.0490593583, 0.019708953, 0.3634155844, 2.1032860162, 0.0861057517, 0.1735660361, 1.5133910481, 0.7858555362, 0.3000131148, 0.3337627573, 0.0036260499, 1.5386413234, 0.5196922389, 0.0221252552, 1.0171151697, 0.0534088166, 6.037787908, 0.4350064365, 0.1634497017, 0.3545179411, 2.3008246523, 0.7625702322, 1.9431704326, 0.6961369276, 2.3726544756, 0.1837198343, 0.9087013201, 2.5477016916, 0.3081949928, 0.1713464632, 2.7297706102, 0.3416923226, 0.0730798705, 4.0107845583, 8.4630191575, 4.3546170435, 1.0655012755, 1.6534489471, 0.0985354973, 0.1940108923, 0.3415280861, 0.2794040892, 0.1657005971, 0.2704552047, 2.3418182855, 0.0426297282, 1.2152488582, 4.6553742047, 0.0068797851, 1.1613183519, 2.2213527952, 0.0565037747, 6.7852754661, 1.0442e-06, 2.842e-07, 0.9529353202, 0.0009844045, 0.0002705734, 0.5068170211, 9.32799e-05, 0.0050518699, 0.3163744815, 2.328e-06, 0.1010587493, 0.2890102379, 0.0041564377, 0.0495269526, 0.0002026765, 0.0358664532, 0.0714121777, 0.3036789915, 1.3220740967, 1.7972997876, 0.0066458178, 0.3052655031, 0.0174305437, 21.9842817264, 0.1070890246, 0.0770894218, 0.1929529483, 0.0561599188, 1.6748429971, 0.0021338646, 1.8890678523, 0.283432044, 0.3134203648, 3.2116908598, 0.0108028571, 0.0860833645, 0.0426724431, 0.3652373073, 0.0287789552, 0.1484349765, 0.5158740953, 0.005979137, 3.3648305163, 0.8763855707, 0.0776875418, 0.9145670668, 0.3963331926, 0.1080226203, 0.0640951379, 0.2278998021, 0.0388755869, 0.1836950254, 0.0461811, 0.053408, 0.0361971, 0.0233326, 0.023417, 0.0390397, 0.0341284001, 0.0389164, 0.016464, 0.0891534, 0.1617310001, 0.0551341, 0.0233262, 0.0911252, 0.0344713001, 0.0771077, 0.0418603001, 0.0200784, 0.0305429, 0.0643851996,
					1.24126910678762, 1.2184237953499, 1.57207707533269, 1.37593685094412, 0.755065443900121, 7.85842191536894, 2.47312230875449, 1.44142625674284, 0.978467912277413, 0.227248844812147, 2.21551678051375, 5.51208197052487, 3.01432016709248, 1.6562495638176, 0.458746912674614, 2.33799112074951, 1.35424048606131, 2.00934347783981, 9.68834518756851, 0.451916794319267, 6.81246018399377, 3.33865551464577, 1.3121700301622, 2.41176328988618, 1.91420790259902, 1.10346056844725, 0.87761105947655, 1.3860121390169, 0.961584192691084, 4.92386682839453, 6.19743849778841, 2.14596406101338, 1.51967567593807, 7.99432285649465, 1.63600796885224, 0.856124897304504, 0.890820306192551, 0.432300548792552, 0.917929117533152, 0.216166037272559, 0.912666803253931, 0.488273343287992, 0.403549792963333, 0.288807503303749, 0.578793711540794, 1.07784974087641, 0.838670114915827, 0.409831127081601, 0.357420746899852, 1.40813159984137, 1.3318097154194, 0.561071724229475, 0.357866239574553, 1.07650079495621, 6.00191102584264, 1.49320558163725, 10.017330817366, 4.40345475789626, 1.4521790561664, 0.337109178564748, 6.05190852431188, 4.32900865295828, 0.89455636623452, 1.80851360960392, 0.624429752512714, 0.564232288255632, 1.90064559617176, 1.248863868961, 0.937880370616514, 0.40752399260009, 1.22130548008116, 1.91061908276291, 0.74719362180685, 0.595481279174004, 1.38082917100197, 6.75978997720454, 8.03277929474211, 1.71296709769163, 0.688343902687262, 0.422494519727629, 0.504494427332431, 0.167512972455925, 1.6953951980808, 0.357343252249955, 0.231719438769158, 0.369372264098046, 1.36297655010811, 2.28642869493161, 4.36115480635558, 0.391055990383483, 2.32013735462963, 2.73556200899535, 1.30918377824208, 0.710372053197474, 1.07146059795775, 0.432622707864552, 2.30191777283007, 1.51328074162521, 0.774493361813496, 1.83705558520706, 0.481140238791115, 1.00843205198373, 1.39189355935829, 0.495319380867629, 0.374682110796213, 6.42089618591429, 1.92029942623162, 6.12345123968018, 2.21619445967418, 3.63668154087443, 2.31937036432372, 1.82735355877736, 3.06377761937176, 1.96998951873875, 0.604749150750474, 0.895375466926981, 1.97766301409123, 1.06574823180769, 1.10791447006064, 3.54659148436289, 5.28925141697764, 1.33634017405606, 3.88525061059222, 1.50668398729448, 1.75570652058377, 2.15765101034714, 1.58399817085847, 0.714748967626738, 1.61366545732856, 2.63447783844427, 1.01920043725065, 2.55137813126603, 3.36284883604624, 0.688272590887225, 1.94853766731376, 8.84799840612482, 0.548857847810693, 1.51701421539628, 0.180852575260598, 0.249658418815177, 1.62751798912531, 0.895908268154618, 0.41983911481111, 0.934975359559877, 0.63019546843603, 0.560464827406078, 1.51831144346793, 0.585192087949017, 1.4680478689711, 3.34484372397723, 0.432605800143879, 0.679112659593982, 0.451420309937647, 0.541176991665778, 0.891261440456541, 1.08949265815113, 0.744762089178451, 2.1579775140421, 0.918359680141276, 0.581811133178276, 0.337446764972448, 7.7587442309146, 0.862679604415627, 1.24522432245413, 0.783544753371045, 1.08991657709568, 10.3848523331335, 0.481910901964747, 0.954722930595868, 0.856431418469122, 4.53772357904054, 4.65018946918032, 0.780701785580677, 0.458606198171997, 0.459453524166091, 2.26274569962909, 0.636693250139687, 0.894057287554733, 0.619332103417392, 0.533322094403035, 14.8729334615191, 3.54580932766672, 0.780108033599127, 4.05845771567534, 1.70397305226754, 0.598549891298567, 0.930523211302821, 3.42422184508655, 0.565896924903265, 1, 0.0770764620135024, 0.0500819370772208, 0.0462377395993731, 0.0537929860758246, 0.0144533387583345, 0.0408923608974345, 0.0633579339160905, 0.0655672355884439, 0.0218802687005936, 0.0591969699027449, 0.0976461276528445, 0.059207941082273, 0.0220695876653368, 0.041350852183426, 0.0476871596856874, 0.0707295165111524, 0.0567759161524817, 0.0127019797647213, 0.0323746050281867, 0.0669190817443274,
					0.551571, 0.509848, 0.635346, 0.738998, 0.147304, 5.42942, 1.02704, 0.528191, 0.265256, 0.0302949, 0.908598, 3.0355, 1.54364, 0.616783, 0.0988179, 1.58285, 0.439157, 0.947198, 6.17416, 0.021352, 5.46947, 1.41672, 0.584665, 1.12556, 0.865584, 0.306674, 0.330052, 0.567717, 0.316954, 2.13715, 3.95629, 0.930676, 0.248972, 4.29411, 0.570025, 0.24941, 0.193335, 0.186979, 0.554236, 0.039437, 0.170135, 0.113917, 0.127395, 0.0304501, 0.13819, 0.397915, 0.497671, 0.131528, 0.0848047, 0.384287, 0.869489, 0.154263, 0.0613037, 0.499462, 3.17097, 0.906265, 5.35142, 3.01201, 0.479855, 0.0740339, 3.8949, 2.58443, 0.373558, 0.890432, 0.323832, 0.257555, 0.893496, 0.683162, 0.198221, 0.103754, 0.390482, 1.54526, 0.315124, 0.1741, 0.404141, 4.25746, 4.85402, 0.934276, 0.210494, 0.102711, 0.0961621, 0.0467304, 0.39802, 0.0999208, 0.0811339, 0.049931, 0.679371, 1.05947, 2.11517, 0.088836, 1.19063, 1.43855, 0.679489, 0.195081, 0.423984, 0.109404, 0.933372, 0.682355, 0.24357, 0.696198, 0.0999288, 0.415844, 0.556896, 0.171329, 0.161444, 3.37079, 1.22419, 3.97423, 1.07176, 1.40766, 1.02887, 0.704939, 1.34182, 0.740169, 0.31944, 0.344739, 0.96713, 0.493905, 0.545931, 1.61328, 2.12111, 0.554413, 2.03006, 0.374866, 0.512984, 0.857928, 0.822765, 0.225833, 0.473307, 1.45816, 0.326622, 1.38698, 1.51612, 0.171903, 0.795384, 4.37802, 0.113133, 1.16392, 0.0719167, 0.129767, 0.71707, 0.215737, 0.156557, 0.336983, 0.262569, 0.212483, 0.665309, 0.137505, 0.515706, 1.52964, 0.139405, 0.523742, 0.110864, 0.240735, 0.381533, 1.086, 0.325711, 0.543833, 0.22771, 0.196303, 0.103604, 3.87344, 0.42017, 0.398618, 0.133264, 0.428437, 6.45428, 0.216046, 0.786993, 0.291148, 2.48539, 2.00601, 0.251849, 0.196246, 0.152335, 1.00214, 0.301281, 0.588731, 0.187247, 0.118358, 7.8213, 1.80034, 0.305434, 2.05845, 0.649892, 0.314887, 0.232739, 1.38823, 0.365369, 0.31473, 0.0866279, 0.043972, 0.0390894, 0.0570451, 0.0193078, 0.0367281, 0.0580589, 0.0832518, 0.0244313, 0.048466, 0.086209, 0.0620286, 0.0195027, 0.0384319, 0.0457631, 0.0695179, 0.0610127, 0.0143859, 0.0352742, 0.0708956,
					0.589718, 0.514347, 0.67416, 0.731152, 0.159054, 5.30821, 1.21324, 0.568449, 0.233527, 0.0379056, 1.03344, 3.02808, 1.62299, 0.657364, 0.0999068, 1.55788, 0.443685, 1.00122, 6.04299, 0.0284956, 5.6037, 1.41993, 0.629768, 1.12717, 0.88357, 0.312544, 0.346823, 0.588609, 0.317684, 2.31211, 3.9337, 0.958529, 0.341479, 4.87366, 0.599188, 0.279542, 0.214596, 0.187262, 0.527321, 0.0390513, 0.198958, 0.125999, 0.124553, 0.0310522, 0.162975, 0.400822, 0.51821, 0.144354, 0.0869637, 0.451124, 0.873266, 0.154936, 0.067443, 0.508952, 3.1554, 0.881639, 5.74119, 2.88102, 0.480308, 0.0719929, 4.19125, 2.45392, 0.381514, 0.854485, 0.320597, 0.255092, 0.887458, 0.660816, 0.198404, 0.0992829, 0.428648, 1.64018, 0.294481, 0.184545, 0.40117, 3.94646, 4.81956, 0.877057, 0.213179, 0.122792, 0.0848492, 0.0458258, 0.485001, 0.109241, 0.0873936, 0.0552962, 0.631713, 1.06458, 2.10414, 0.0832422, 1.14516, 1.51861, 0.711498, 0.204905, 0.444152, 0.109081, 0.913179, 0.720567, 0.254626, 0.722123, 0.111722, 0.422851, 0.588203, 0.179858, 0.165205, 3.52499, 1.35611, 3.90127, 1.09965, 1.35221, 0.87908, 0.822025, 1.33618, 0.876688, 0.321774, 0.351913, 1.05314, 0.554077, 0.563999, 1.54694, 2.24161, 0.594177, 2.06787, 0.395176, 0.522957, 0.829315, 0.889765, 0.236489, 0.54992, 1.48876, 0.351564, 1.45173, 1.56873, 0.188237, 0.802531, 4.02507, 0.135395, 1.24086, 0.0746093, 0.142159, 0.728065, 0.208163, 0.176397, 0.366467, 0.261223, 0.259584, 0.706082, 0.159261, 0.565299, 1.58681, 0.135024, 0.528249, 0.118584, 0.270321, 0.386714, 1.05269, 0.326191, 0.481954, 0.210494, 0.209621, 0.108982, 4.31772, 0.44009, 0.427718, 0.155623, 0.437069, 6.49269, 0.212945, 0.742154, 0.286443, 2.42261, 1.92496, 0.282892, 0.193323, 0.155419, 1.10899, 0.32893, 0.588443, 0.190095, 0.119749, 7.48376, 1.82105, 0.300343, 2.03324, 0.653015, 0.325745, 0.23769, 1.4088, 0.396884, 0.353358, 0.0866279, 0.043972, 0.0390894, 0.0570451, 0.0193078, 0.0367281, 0.0580589, 0.0832518, 0.0244313, 0.048466, 0.086209, 0.0620286, 0.0195027, 0.0384319, 0.0457631, 0.0695179, 0.0610127, 0.0143859, 0.0352742, 0.0708956),
				nrow=210,
				dimnames=list(c("A/R", "A/N", "R/N", "A/D", "R/D", "N/D", "A/C", "R/C", "N/C", "D/C", "A/Q", "R/Q", "N/Q", "D/Q", "C/Q", "A/E", "R/E", "N/E", "D/E", "C/E", "Q/E", "A/G", "R/G", "N/G", "D/G", "C/G", "Q/G", "E/G", "A/H", "R/H", "N/H", "D/H", "C/H", "Q/H", "E/H", "G/H", "A/I", "R/I", "N/I", "D/I", "C/I", "Q/I", "E/I", "G/I", "H/I", "A/L", "R/L", "N/L", "D/L", "C/L", "Q/L", "E/L", "G/L", "H/L", "I/L", "A/K", "R/K", "N/K", "D/K", "C/K", "Q/K", "E/K", "G/K", "H/K", "I/K", "L/K", "A/M", "R/M", "N/M", "D/M", "C/M", "Q/M", "E/M", "G/M", "H/M", "I/M", "L/M", "K/M", "A/F", "R/F", "N/F", "D/F", "C/F", "Q/F", "E/F", "G/F", "H/F", "I/F", "L/F", "K/F", "M/F", "A/P", "R/P", "N/P", "D/P", "C/P", "Q/P", "E/P", "G/P", "H/P", "I/P", "L/P", "K/P", "M/P", "F/P", "A/S", "R/S", "N/S", "D/S", "C/S", "Q/S", "E/S", "G/S", "H/S", "I/S", "L/S", "K/S", "M/S", "F/S", "P/S", "A/T", "R/T", "N/T", "D/T", "C/T", "Q/T", "E/T", "G/T", "H/T", "I/T", "L/T", "K/T", "M/T", "F/T", "P/T", "S/T", "A/W", "R/W", "N/W", "D/W", "C/W", "Q/W", "E/W", "G/W", "H/W", "I/W", "L/W", "K/W", "M/W", "F/W", "P/W", "S/W", "T/W", "A/Y", "R/Y", "N/Y", "D/Y", "C/Y", "Q/Y", "E/Y", "G/Y", "H/Y", "I/Y", "L/Y", "K/Y", "M/Y", "F/Y", "P/Y", "S/Y", "T/Y", "W/Y", "A/V", "R/V", "N/V", "D/V", "C/V", "Q/V", "E/V", "G/V", "H/V", "I/V", "L/V", "K/V", "M/V", "F/V", "P/V", "S/V", "T/V", "W/V", "Y/V", "FreqA", "FreqR", "FreqN", "FreqD", "FreqC", "FreqQ", "FreqE", "FreqG", "FreqH", "FreqI", "FreqL", "FreqK", "FreqM", "FreqF", "FreqP", "FreqS", "FreqT", "FreqW", "FreqY", "FreqV"),
					c("AB", "BLOSUM62", "cpREV", "cpREV64", "Dayhoff", "DCMut-Dayhoff", "DCMut-JTT", "DEN", "FLAVI", "FLU", "gcpREV", "HIVb", "HIVw", "JTT", "LG", "MtArt", "mtDeu", "mtInv", "mtMam", "mtMet", "mtOrt", "mtREV", "mtVer", "MtZoa", "PMB", "Q.bird", "Q.insect", "Q.LG", "Q.mammal", "Q.pfam", "Q.plant", "Q.yeast", "rtREV", "stmtREV", "VT", "WAG", "WAGstar")))
		} else {
			MODELS <- MODELS$Nucleotide
		}
		g <- !(submodels %in% MODELS)
		if (any(g))
			stop(paste("Unrecognized model",
				ifelse(sum(g) > 1, "s", ""),
				": ",
				paste(model[g], collapse=", "),
				sep=""))
		rates <- as.numeric(rates)
		if (any(floor(rates) != rates, na.rm=TRUE))
			stop('The number rates in the discrete Gamma distribution (i.e., "+G4") should be an integer value.')
		rates <- as.integer(rates)
		if (any(rates > 10, na.rm=TRUE))
			stop('Up to 10 rates are allowed in the discrete Gamma distribution (i.e., "+G10").')
		if (any(rates < 2, na.rm=TRUE))
			stop('A minimum of two rates are required for the discrete Gamma distribution (i.e., "+G2").')
	}
	
	if (is.null(myDistMatrix)) {
		maxDist <- 3 # maximum distance allowed in distance matrix
		myDistMatrix <- DistanceMatrix(myXStringSet,
			correction="F81",
			penalizeGapLetterMatches=method != 7 && indels,
			processors=processors,
			verbose=FALSE)
		myDistMatrix[myDistMatrix > maxDist] <- maxDist
		# impute missing distances with ultrametric method
		w <- which(is.na(myDistMatrix), arr.ind=TRUE)
		for (k in seq_len(nrow(w))) {
			myDistMatrix[w[k, 1], w[k, 2]] <- min(maxDist,
				pmax(myDistMatrix[w[k, 1],],
					myDistMatrix[w[k, 2],]),
				na.rm=TRUE)
		}
	}
	
	if (verbose && method != 3 && method != 7) {
		pBar <- txtProgressBar(min=0, max=100, initial=0, style=ifelse(interactive(), 3, 1))
		time.1 <- Sys.time()
	} else {
		pBar <- NULL
	}
	
	myClusters <- .Call("cluster",
		myDistMatrix,
		ifelse(method == 3 || method == 7,
			-Inf,
			cutoff[1]),
		ifelse(method == 3 || method == 7,
			1L, # NJ
			method),
		ifelse(is(myDistMatrix, "matrix"),
			dim,
			-dim),
		verbose && method != 3 && method != 7,
		pBar,
		processors,
		PACKAGE="DECIPHER")
	
	if (verbose && method != 3 && method != 7) {
		setTxtProgressBar(pBar, 100)
		close(pBar)
	}
	
	if (method == 3 ||
		method == 7 ||
		(reconstruct && type > 1)) {
		myClusters <- .reorderClusters(myClusters,
			all=method != 3)
		Z <- as.matrix(myXStringSet)
		if (method != 7) {
			if (typeX == 3L) {
				.optimizeModel <- .optimizeModelAA
				.giveParams <- .giveParamsAA
				m <- matrix(NA_real_,
					nrow=length(model),
					ncol=216,
					dimnames=list(model,
						c(rownames(ProtModels),
							"FreqI", "Indels",
							"alpha",
							"-LnL", "AICc", "BIC")))
			} else {
				.optimizeModel <- .optimizeModelDNA
				.giveParams <- .giveParamsDNA
				m <- matrix(NA_real_,
					nrow=length(model),
					ncol=15,
					dimnames=list(model,
						c("FreqA", "FreqC", "FreqG", "FreqT", "FreqI",
							"A/G", "C/T", "A/C", "A/T", "C/G", "Indels",
							"alpha",
							"-LnL", "AICc", "BIC")))
			}
			weights_ML <- tabulate(selfmatch(as.data.frame(Z)), ncol(Z))
			if (quadrature) {
				.rates <- .rates1 # use Laguerre quadrature
			} else {
				.rates <- .rates2 # use equal bins
			}
			
			if (typeX == 3L) {
				a <- alphabetFrequency(myXStringSet, collapse=TRUE)
				if (indels) {
					a <- a[c(1:20, 28)]
					N <- sum(a)/dim # average number of sites
					a <- a/sum(a)
				} else {
					a <- a[1:20]
					N <- sum(a)/dim # average number of sites
					a <- c(a/sum(a), 0)
				}
			} else {
				a <- consensusMatrix(myXStringSet)
				defaults <- .guessParams(a, indels)
				if (indels) {
					N <- sum(a[c(1:14, 16),, drop=FALSE])/dim # average number of sites
				} else {
					N <- sum(a[1:14,, drop=FALSE])/dim # average number of sites
				}
			}
			if (verbose) {
				cat("Optimizing model parameters:")
				flush.console()
			}
			if (method == 3) {
				.minimize <- function(x, branches=integer(), lengths=numeric()) {
					# given branch lengths return -LnL
					myClusters[, 4:5] <- x
					.Call("clusterML",
						myClusters,
						myXStringSet,
						model_params,
						branches,
						lengths,
						0,
						typeX,
						weights_ML,
						processors,
						PACKAGE="DECIPHER")
				}
			}
			
			resTime <- matrix(NA_real_, nrow=length(model), ncol=processors)
			for (i in seq_along(model)) {
				start <- Inf
				it <- 1L
				if (verbose && interactive()) {
					cat("\n", model[i], sep="")
					flush.console()
				}
				if (typeX == 3L) {
					if (indels) {
						if (grepl("\\+F", model[i], ignore.case=TRUE)) {
							defaults <- c(ProtModels[1:190, submodels[i]], a, 1, 1)
						} else {
							defaults <- c(ProtModels[, submodels[i]], a[21], 1, 1)
						}
					} else {
						if (grepl("\\+F", model[i], ignore.case=TRUE)) {
							defaults <- c(ProtModels[1:190, submodels[i]], a[1:20], 0, 0, 1)
						} else{
							defaults <- c(ProtModels[, submodels[i]], 0, 0, 1)
						}
					}
				}
				repeat {
					if (verbose && method == 3 && interactive()) {
						cat("\r", model[i],
							" (iteration ", it, ")",
							sep="")
						flush.console()
					}
					
					if (it == 1L) {
						optProcessors <- processors
					} else if (it == 2L) {
						optProcessors <- 1L
					} else {
						w <- which(is.na(resTime[i,]))
						if (length(w) == 0L) {
							optProcessors <- which.min(resTime[i,])
						} else if (length(w) == 1L) {
							optProcessors <- w
						} else {
							o <- order(resTime[i,])
							optProcessors <- w[which.min(abs(mean(o[1:2]) - w))]
						}
					}
					
					gc(FALSE) # necessary for repeatable timing
					time.3 <- Sys.time()
					temp <- .optimizeModel(myClusters,
						model[i],
						myXStringSet,
						N,
						method == 3, # scaleTree
						.rates,
						if (it == 1L) {
							defaults
						} else {
							m[model[i],]
						},
						weights_ML,
						factr=min(1e10,
							ifelse(start < Inf,
								min(absTol/start, relTol)/.Machine$double.eps,
								relTol/.Machine$double.eps)/10),
						processors=optProcessors)
					time.4 <- Sys.time()
					if (is.na(resTime[i, optProcessors])) {
						resTime[i, optProcessors] <- difftime(time.4, time.3, units='secs')/temp[2L]
					} else {
						resTime[i, optProcessors] <- (resTime[i, optProcessors] + difftime(time.4, time.3, units='secs')/temp[2L])/2
					}
					m[model[i],] <- temp[-1:-2]
					
					if (method == 3) {
						myClusters[, 4:5] <- temp[1L]*myClusters[, 4:5] # scale tree
						
						if (start - m[model[i], "-LnL"] < absTol)
							break
						
						start <- m[model[i], "-LnL"]
						model_params <- .giveParams(m[model[i],],
							rownames(m)[i],
							.rates)
						
						params <- as.vector(myClusters[, 4:5])
						params <- .globalBranches(.minimize, params)
						myClusters[, 4:5] <- params
					} else {
						break
					}
					
					it <- it + 1L
				}
				
				if (verbose) {
					cat("\r", model[i],
						paste(rep(" ",
								max(nchar(rownames(m))) - nchar(model[i]) + 1),
							collapse=""),
						"-ln(L) = ",
						round(m[model[i], "-LnL"], 0),
						", AICc = ",
						round(m[model[i], "AICc"], 0),
						", BIC = ",
						round(m[model[i], "BIC"], 0),
						sep="")
					flush.console()
				}
			}
			
			if (length(model) > 1) { # choose the best model
				w <- which.min(m[, informationCriterion])
				optProcessors <- which.min(resTime[w,])
				m <- m[w,, drop=FALSE]
				model <- model[w]
				if (verbose)
					cat("\n\nThe selected model was:  ",
						model,
						ifelse(method == 3 || method == 8, "\n\n", "\n"),
						sep="")
			} else {
				optProcessors <- which.min(resTime[1,])
				if (verbose)
					cat("\n\n")
			}
			
			model_params <- .giveParams(m,
				rownames(m),
				.rates)
		}
		
		if (nrow(myClusters) > 2 &&
			(method == 3 || method == 7)) {
			# initialize functions requiring the local environment
			if (method == 3) {
				.minimize <- function(x, branches=integer(), lengths=numeric()) {
					# given branch lengths return -LnL
					myClusters[, 4:5] <- x
					LnL <- .Call("clusterML",
						myClusters,
						myXStringSet,
						model_params,
						branches,
						lengths,
						0,
						typeX,
						weights_ML,
						optProcessors,
						PACKAGE="DECIPHER")
					
					w <- which.min(LnL)
					if (LnL[w] < .best - epsilon) {
						.best <<- LnL[w]
						if (verbose) {
							if (is.infinite(.overall))
								.overall <<- LnL[1]
							.printLine(LnL[w])
						}
					}
					
					LnL
				}
				
				.NNI <- function(myClusters) {
					out <- .localBranches(myClusters,
						myXStringSet,
						model_params,
						weights_ML,
						optProcessors)
					delta <- out[[4]] - epsilon - out[[2L]]
					o <- which(delta > 0)
					o <- o[order(out[[2L]][o])]
					myClustersTemp <- myClusters
					NNIs <- 0L
					for (i in o) {
						center <- out[[1L]][i, 1L]
						left <- myClusters[center, 7L]
						right <- myClusters[center, 8L]
						if (myClustersTemp[center, 7L] != left ||
							myClustersTemp[center, 8L] != right ||
							(left > 0 &&
							(myClustersTemp[left, 7L] != myClusters[left, 7L] ||
							myClustersTemp[left, 8L] != myClusters[left, 8L])) ||
							(right > 0 &&
							(myClustersTemp[right, 7L] != myClusters[right, 7L] ||
							myClustersTemp[right, 8L] != myClusters[right, 8L])))
							next # incompatible NNI
						
						NNIs <- NNIs + 1L
						side <- out[[1L]][i, 2L]
						flip <- out[[1L]][i, 3L]
						down <- myClusters[center, 6 + side]
						if (center == nrow(myClusters)) { # root
							opposite <- myClusters[center, 9 - side] # opposite-left
							up <- opposite + nrow(myClusters) # opposite-right
						} else {
							up <- match(center, myClusters[, 7:8])
						}
						myClusters[, 4:5][up] <- out[[3L]][5L, i]
						
						if (side == 1L) { # center branch is left
							if (flip > 0) { # swap opposite with down-left
								myClusters[center, 4L] <- out[[3L]][1L, i]
								myClusters[down, 4L] <- out[[3L]][2L, i]
								myClusters[down, 5L] <- out[[3L]][4L, i]
								if (center == nrow(myClusters)) { # root
									myClusters[opposite, 4L] <- out[[3L]][3L, i]
									myClusters <- .swapBranches(myClusters,
										opposite, 7, # opposite-left
										down, 7) # down-left
								} else {
									myClusters[center, 5L] <- out[[3L]][3L, i]
									myClusters <- .swapBranches(myClusters,
										center, 8, # opposite
										down, 7) # down-left
								}
							} else { # swap opposite with down-right
								myClusters[center, 4L] <- out[[3L]][1L, i]
								myClusters[down, 4L] <- out[[3L]][3L, i]
								myClusters[down, 5L] <- out[[3L]][2L, i]
								if (center == nrow(myClusters)) { # root
									myClusters[opposite, 4L] <- out[[3L]][4L, i]
									myClusters <- .swapBranches(myClusters,
										opposite, 7, # opposite-left
										down, 8) # down-right
								} else {
									myClusters[center, 5L] <- out[[3L]][4L, i]
									myClusters <- .swapBranches(myClusters,
										center, 8, # opposite
										down, 8) # down-right
								}
							}
						} else { # center branch is right
							if (flip > 0) { # swap opposite with down-left
								myClusters[center, 5L] <- out[[3L]][1L, i]
								myClusters[down, 4L] <- out[[3L]][2L, i]
								myClusters[down, 5L] <- out[[3L]][4L, i]
								if (center == nrow(myClusters)) { # root
									myClusters[opposite, 4L] <- out[[3L]][3L, i]
									myClusters <- .swapBranches(myClusters,
										opposite, 7, # opposite-left
										down, 7) # down-left
								} else {
									myClusters[center, 4L] <- out[[3L]][3L, i]
									myClusters <- .swapBranches(myClusters,
										center, 7, # opposite
										down, 7) # down-left
								}
							} else { # swap opposite with down-right
								myClusters[center, 5L] <- out[[3L]][1L, i]
								myClusters[down, 4L] <- out[[3L]][3L, i]
								myClusters[down, 5L] <- out[[3L]][2L, i]
								if (center == nrow(myClusters)) { # root
									myClusters[opposite, 4L] <- out[[3L]][4L, i]
									myClusters <- .swapBranches(myClusters,
										opposite, 7, # opposite-left
										down, 8) # down-right
								} else {
									myClusters[center, 4L] <- out[[3L]][4L, i]
									myClusters <- .swapBranches(myClusters,
										center, 7, # opposite
										down, 8) # down-right
								}
							}
						}
					}
					
					list(myClusters,
						NNIs,
						cbind(out[[1L]], out[[2L]]),
						out[[4]])
				}
			} else {
				.minimize <- function(C) {
					score <- .Sankoff(C,
						Z, # integer encoded XStringSet
						S, # substitution matrix
						weights_MP,
						scoreOnly=TRUE,
						processors=processors)
					
					if (score < .best - epsilon) {
						.best <<- score
						if (verbose) {
							if (is.infinite(.overall))
								.overall <<- score
							.printLine(score)
						}
					}
					
					score
				}
				
				.NNI <- function(myClusters) {
					C <- myClusters[, 7:8]
					res <- .Sankoff(C,
						Z,
						S,
						weights_MP,
						TRUE,
						-1, # add
						processors)
					
					NNIs <- 0L
					count <- 1L
					Corg <- C
					for (j in rev(seq_len(nrow(C)))) {
						if (Corg[j, 1] > 0) {
							count <- count + 1L
							if (res[1] > res[count] &&
								C[j, 1] > 0 &&
								all(C[c(j, C[j, 1]),] == Corg[c(j, C[j, 1]),])) {
								# swap left-left with right
								NNIs <- NNIs + 1L
								temp <- C[C[j, 1], 1]
								C[C[j, 1], 1] <- C[j, 2]
								C[j, 2] <- temp
								C <- .reorder2(C)
							}
							
							count <- count + 1L
							if (res[1] > res[count] &&
								C[j, 1] > 0 &&
								all(C[c(j, C[j, 1]),] == Corg[c(j, C[j, 1]),])) {
								# swap left-right with right
								NNIs <- NNIs + 1L
								temp <- C[C[j, 1], 2]
								C[C[j, 1], 2] <- C[j, 2]
								C[j, 2] <- temp
								C <- .reorder2(C)
							}
						}
						
						if (Corg[j, 2] > 0) {
							count <- count + 1L
							if (res[1] > res[count] &&
								C[j, 2] > 0 &&
								all(C[c(j, C[j, 2]),] == Corg[c(j, C[j, 2]),])) {
								# swap right-left with left
								NNIs <- NNIs + 1L
								temp <- C[C[j, 2], 1]
								C[C[j, 2], 1] <- C[j, 1]
								C[j, 1] <- temp
								C <- .reorder2(C)
							}
							
							count <- count + 1L
							if (res[1] > res[count] &&
								C[j, 2] > 0 &&
								all(C[c(j, C[j, 2]),] == Corg[c(j, C[j, 2]),])) {
								# swap right-right with left
								NNIs <- NNIs + 1L
								temp <- C[C[j, 2], 2]
								C[C[j, 2], 2] <- C[j, 1]
								C[j, 1] <- temp
								C <- .reorder2(C)
							}
						}
					}
					
					res <- .Sankoff(C,
						Z,
						S,
						weights_MP,
						FALSE,
						processors=processors)
					
					myClusters[, 7:8] <- C
					myClusters[, 4:5] <- res[[3L]]/N # changes per site
					list(myClusters,
						NNIs,
						NULL,
						res[[1L]])
				}
			}
			
			# print progress of likelihood maximization
			if (verbose) {
				.printLine <- function(value, print=interactive()) {
					if (!print)
						return(value)
					change <- (value - .overall)/value
					if (isTRUE(all.equal(change, 0))) {
						sign <- ""
					} else {
						change <- round(100*change, 3)
						if (change > 0) {
							sign <- "+"
						} else if (change == 0) {
							sign <- "~"
						} else {
							sign <- "-"
						}
					}
					change <- abs(change)
					cat(ifelse(interactive(),
						"\r",
						""),
						ifelse(method == 3,
							"-ln(L) = ",
							"score = "),
						formatC(value,
							digits=1,
							format="f"),
						" (",
						sign,
						formatC(change,
							digits=3,
							format="f"),
						"%)",
						ifelse(is.na(.Grafts), # don't report Climbs if Grafts
							ifelse(is.na(.Climbs),
								"",
								paste(",",
									.Climbs,
									ifelse(.Climbs == 1, "Climb", "Climbs"))),
							paste(", ",
								.Grafts,
								ifelse(.Grafts == 1, " Graft", " Grafts"),
								" of ",
								.totGrafts,
								sep="")),
						"  ",
						sep="")
					flush.console()
					value
				}
			}
			
			# initialize overall parameters
			relTol <- 0.0001 # relative convergence tolerance (0, 1]
			absTol <- 0.1 # absolute convergence tolerance (> 0)
			epsilon <- 1e-5 # threshold to accept changes
			fracParams <- 1.001 # attempt parameter optimization within fracParams*best score
			startingTrees <- c(1L, 5L) # number of initial trees to compare per iteration
			generationSize <- 100L # maximum number of regrowths per generation (also used to set maxGenerations)
			max_iterations <- c(100L, maxGenerations*generationSize, 1000L) # maximum number of iterations per phase
			max_time <- maxTime*max_iterations/sum(max_iterations) # maximum number of hours per phase
			
			# initialize initial tree parameters
			if (method == 3) {
				noiseLevel <- c(0, 0.5) # Gaussian noise added to substitution matrix (>= 0)
			} else {
				noiseLevel <- c(0, 0) # substitution matrix is fixed
			}
			nSamps <- Inf # number of resamples of parsimony weights (> 0 or Inf for zero noise)
			observations <- 10L # number of repeated observations of best score to skip phase
			
			# initialize climb parameters
			doClimbs <- TRUE
			
			# initialize regrowth parameters
			waitAttempts <- 10L # number of attempts (per generation) to wait before sampling all trees for PCA
			maxSD <- 0.001 # normalized standard deviation of scores to keep iterating after observations (>= 0)
			noiseEstimation <- c(1, 2.5) # Gaussian noise added to estimates (multiple of standard error >= 0)
			base <- 2 # base of weight function
			fracSamp <- 0.25 # fraction of previous trees to sample for PCA (0, 1]
			numSamp <- ceiling(fracSamp*max_iterations[1L]) # fixed number of trees to sample
			tolPCA <- seq(0.5, 0.5, length.out=generationSize) # tolerances for PCA each generation (0, 1]
			addNoise <- 0.5 # amount of noise to add when converging to the same tree (>= 0)
			popSize <- startingTrees[2L] # number of projected cophenetic matrices to store
			
			# initialize grafting parameters
			fracGraft <- 0.1 # fraction of high scoring trees for first grafts [0, 1]
			doGrafts <- FALSE
			.Grafts <- NA # number of tree recombinations
			
			# initialize shake parameters
			.Shakes <- NA # number of shakes (or NA for none)
			doShakes <- FALSE
			repsShakes <- 50L # number of shake replicates
			fracRandomNNIs <- 0.7 # fraction of random NNIs (0, 1]
			
			# initialize variables for maximum parsimony
			if (typeX == 1) {
				states <- DNA_BASES
			} else if (typeX == 2) {
				states <- RNA_BASES
			} else {
				states <- AA_STANDARD
			}
			if (method == 3 && indels)
				states <- c(states, "-")
			Z <- matrix(match(Z, states), # integer encode
				nrow=nrow(Z),
				ncol=ncol(Z))
			r <- apply(Z,
				2L,
				function(x) {
					x <- x[!is.na(x)]
					if (length(x) > 0) {
						max(x) - min(x)
					} else {
						0
					}
				})
			r <- which(r > 0) # informative sites
			weights_MP <- tabulate(r[selfmatch(as.data.frame(Z[, r, drop=FALSE]))], ncol(Z))
			if (method == 7) {
				S <- eval(sexpr,
					envir=list(x=states),
					enclos=parent.frame())
				if (!is(S, "matrix"))
					stop("costMatrix must evaluate to a matrix.")
				if (nrow(S) != ncol(S))
					stop("costMatrix must evaluate to a square matrix.")
				if (!is.numeric(S))
					stop("costMatrix must evaluate to a numeric matrix.")
				mode(S) <- "numeric"
				if (any(is.na(S)))
					stop("costMatrix evaluates to a matrix containing NA values.")
				if (any(S < 0))
					stop("costMatrix evaluates to a matrix containing negative values.")
				N <- sum(!is.na(Z))/nrow(Z) # denominator of parsimony branch length
			}
			
			# maximize likelihood of tree
			final <- FALSE # final pass
			it <- 0L # iteration number
			I <- 0L # estimate number
			cand <- 0L # current candidate
			t <- list() # tested topologies
			gen <- 0L # current generation
			allow <- TRUE # allow parameter optimization at end of iteration
			maxTrees <- sum(max_iterations)
			graft <- logical(sum(max_iterations))
			Scores <- numeric(maxTrees)
			Trees <- vector("list", maxTrees)
			Cophenetic <- vector("list", max_iterations[1L] + max_iterations[2L])
			if (method == 3)
				Models <- vector("list", maxTrees)
			currentTime <- Sys.time()
			if (verbose)
				cat("PHASE 1 OF 3: INITIAL TREES\n")
			while (!final) {
				if (it >= max_iterations[1L]) {
					if (gen == maxGenerations &&
						cand == generationSize) {
						currentTime <- Sys.time()
						if (verbose)
							cat("\n\nPHASE 3 OF 3: SHAKEN TREES\n")
						doGrafts <- TRUE # note: grafting sets best tree for shakes
						doClimbs <- FALSE # note: used to delineate iteration type
						doShakes <- TRUE
						.Shakes <- 0L # number of shakes (or NA for none)
						.totShakes <- 0L
						allow <- FALSE
						
						w <- (it - generationSize + 1L):it
						w <- w[Scores[w] > 0]
						graft[w[which.min(Scores[w])]] <- TRUE
						
						w <- which(Scores[seq_len(it)] > 0)
						graft[w[Scores[w] < quantile(Scores[w], fracGraft)]] <- TRUE
						w <- w[which.min(Scores[w])]
						.best <- Scores[w]
						
						outgroups <- Cophenetic[[w]]
						outgroups <- .rowSums(outgroups)
						outgroups <- order(outgroups)
						outgroups <- outgroups[c(1L, length(outgroups))]
						Cophenetic[] <- list(NULL)
					} else {
						if (cand == generationSize) {
							w <- (it - generationSize + 1L):it
							w <- w[Scores[w] > 0]
							W <- which.min(Scores[w])
							graft[w[W]] <- TRUE
							Cophenetic[w[-W]] <- list(NULL)
							observed <- Scores[graft]
							observed <- sum((observed - min(observed))/min(observed) < relTol)
							if (observed >= observations) {
								w <- (max_iterations[1L] + 1L):it
								w <- w[Scores[w] > 0]
								if (sd(Scores[w])/mean(Scores[w]) < maxSD) {
									gen <- maxGenerations
									it <- sum(max_iterations[1:2])
									next # skip remaining generations
								} else {
									observed <- observations - 1L
								}
							}
							attempt <- 1L
							gen <- gen + 1L
							cand <- 0L
							W <- integer()
							offset <- 0
						} else if (it > max_iterations[1L]) {
							W <- (it - cand + 1L):it
							offset <- min(Scores[W])
							offset <- sum((Scores[W] - offset)/offset < relTol)
							offset <- addNoise*(offset - 1) # add noise as needed
						} else {
							currentTime <- Sys.time()
							observed <- 0L
							attempt <- 1L
							gen <- gen + 1L
							W <- integer()
							offset <- 0
						}
						
						W <- c(seq_len(max_iterations[1L]), W)
						w <- W[which.min(Scores[W])]
						if (w == it && # last iteration was best
							sum(Scores[W] <= Scores[w] + absTol) == 1L) {
							attempt <- 1L
						} else if (cand > 0L) {
							attempt <- attempt + 1L
							if (attempt > waitAttempts) {
								w <- W[which.min(Scores[W])]
								W <- seq_len(it) # use all trees
								W <- W[lengths(Cophenetic)[W] > 0]
								if (min(Scores[W]) >= Scores[w] - absTol) { # best overall occurred this generation
									it <- it + generationSize - cand
									cand <- generationSize
									next
								}
							}
						}
						
						if (it > max_iterations[1L] &&
							difftime(Sys.time(), currentTime, units="hours") > max_time[2L]) {
							it <- it + generationSize - cand
							cand <- generationSize
							next
						}
						.overall <- Scores[w]
						if (method == 3) {
							m <- Models[[w]][[1]]
							model_params <- Models[[w]][[2]]
						}
						
						it <- it + 1L
						.Climbs <- 0L # number of Climbs (or NA for none)
						cand <- cand + 1L
						I <- I + 1L
						.best <- Inf
						doClimbs <- TRUE
						
						if (verbose) {
							cat(ifelse(cand == 1L,
									paste("\n\nPHASE 2 OF 3: REGROW GENERATION ",
										gen,
										" OF ",
										ifelse(gen + observations - observed <= maxGenerations,
											paste(gen + observations - observed - 1L, "TO "),
											""),
										maxGenerations,
										"\n\n",
										sep=""),
									"\n"),
								"2/3. Optimizing regrown tree #",
								cand,
								" of ",
								if (min(Scores[Scores > 0]) >= Scores[w] - absTol) { # best overall
									ifelse(cand + waitAttempts - attempt >= generationSize,
										generationSize,
										paste(cand + waitAttempts - attempt,
											"to",
											generationSize))
								} else {
									ifelse(cand + waitAttempts >= generationSize,
										generationSize,
										paste(cand + waitAttempts,
											"to",
											generationSize))
								},
								":\n",
								sep="")
							flush.console()
						}
						
						if (w == it - 1L ||
							cand == 1L ||
							I + startingTrees[2L] - 1L > popSize) {
							w <- W[sample(seq_along(W),
								numSamp,
								prob=base^(-1000*(Scores[W]/Scores[w] - 1)))]
							for (i in seq_along(w)) {
								d <- Cophenetic[[w[i]]]
								if (i == 1L)
									D <- matrix(NA_real_,
										nrow=length(w),
										ncol=length(d))
								D[i,] <- d
							}
							
							I <- 1L
							E <- .estimate(Scores[w],
								D,
								base,
								tol=tolPCA[cand],
								mult=noiseEstimation + offset,
								N=min((generationSize - cand + 1L)*startingTrees[2L], popSize))
						}
						
						candScore <- Inf
						for (I in I:(I + startingTrees[2L] - 1L)) {
							# perform regrowth
							myClusters <- .Call("cluster",
								E[I,],
								-Inf, # cutoff
								1L, # NJ
								-dim,
								FALSE,
								NULL,
								processors,
								PACKAGE="DECIPHER")
							
							if (method == 3) { # optimize all branch lengths
								params <- as.vector(myClusters[, 4:5])
								params <- .globalBranches(.minimize, params)
								if (.best < candScore) {
									myClusters[, 4:5] <- params
									myClustersTemp <- myClusters
									candScore <- .best
								}
							} else { # method == 7
								.best <- .Sankoff(myClusters[, 7:8],
									Z, # integer encoded XStringSet
									S, # substitution matrix
									weights_MP,
									scoreOnly=TRUE,
									processors=processors)
								if (.best < candScore) {
									myClustersTemp <- myClusters
									candScore <- .best
								}
							}
						}
						
						myClusters <- myClustersTemp
						
						if (verbose)
							.printLine(.best)
					}
				} else {
					.Climbs <- 0L # number of Climbs (or NA for none)
					allow <- TRUE
					
					if (it == 0L) {
						if (method == 3) {
							.best <- m[, "-LnL"]
						} else {
							.best <- .Sankoff(myClusters[, 7:8],
								Z, # integer encoded XStringSet
								S, # substitution matrix
								weights_MP,
								scoreOnly=TRUE,
								processors=processors)
						}
						.overall <- .best
						observed <- 0L
					} else {
						Score <- Scores[seq_len(it)]
						w <- which.min(Score)
						observed <- sum((Score - Score[w])/Score[w] < relTol)
						if (observed >= observations ||
							(it > 2L &&
							difftime(Sys.time(), currentTime, units="hours") > max_time[1L])) {
							max_iterations[1L] <- it
							if (repsShakes > it)
								repsShakes <- it
							if (numSamp > it)
								numSamp <- it
							next
						}
						.overall <- Scores[w]
						if (method == 3) {
							m <- Models[[w]][[1]]
							model_params <- Models[[w]][[2]]
						}
					}
					it <- it + 1L
					
					if (verbose) {
						cat("\n1/3. Optimizing initial tree #",
							it,
							" of ",
							ifelse(it + observations - observed <= max_iterations[1L],
								paste(it + observations - observed - 1L, "to "),
								""),
							max_iterations[1L],
							":\n",
							sep="")
						flush.console()
					}
					
					if (it > 1L) {
						if (method == 3) {
							.best <- Inf
							S <- .makeS(model_params,
								length(states),
								v=mean(Trees[[w]][, 4:5]),
								noise=noiseLevel[2L]*(it - 2)/max_iterations[1L] + noiseLevel[1L])
						}
						candScore <- Inf
						for (j in seq_len(startingTrees[1L])) {
							if (nSamps < Inf) {
								weights <- sample(seq_along(weights_MP),
									nSamps*length(weights_MP),
									replace=TRUE,
									prob=weights_MP)
								weights <- tabulate(weights, length(weights_MP))
							} else {
								weights <- weights_MP
							}
							
							temp <- .clusterMP(Z,
								S,
								sample(dim), # seed
								weights=weights,
								processors=processors)
							temp[[4L]] <- temp[[4L]]/N # changes per site
							
							if (method == 3) { # optimize all branch lengths
								myClustersTemp <- myClusters
								myClusters[, 4:5] <- temp[[4]]
								myClusters[, 7:8] <- temp[[1]]
								params <- as.vector(temp[[4]])
								params <- .globalBranches(.minimize, params)
								if (.best < candScore) {
									myClusters[, 4:5] <- params
									candScore <- .best
								} else {
									myClusters <- myClustersTemp
								}
							} else { # method == 7
								.best <- temp[[2]]
								if (.best < candScore) {
									myClusters[, 4:5] <- temp[[4]]
									myClusters[, 7:8] <- temp[[1]]
									candScore <- .best
									if (verbose)
										.printLine(.best)
								}
							}
						}
					} else if (verbose) {
						.printLine(.best)
					}
				}
				
				repeat {
					# perform tree fusions
					if (doGrafts) {
						.Grafts <- 0L
						.totGrafts <- 0L
						
						w <- which(Scores[seq_len(it)] > 0)
						w <- w[which.min(Scores[w])]
						myClusters <- Trees[[w]]
						.overall <- .best <- best <- Scores[w]
						if (method == 3) {
							m <- Models[[w]][[1]]
							model_params <- Models[[w]][[2]]
						}
						
						s <- which(graft)
						graft[s] <- FALSE
						s <- s[s != w]
						
						if (verbose) {
							cat("\nGrafting",
								length(s),
								ifelse(length(s) > 1,
									"trees",
									"tree"),
								"to the best tree:\n")
							.printLine(.best)
							flush.console()
						}
						
						for (outgroup in outgroups) {
							myClusters <- .root(myClusters, outgroup)
							
							t1 <- .extractClades(myClusters)
							s1 <- lapply(t1, sort)
							myClustersTemp <- myClusters
							
							for (i in s) {
								# find untested topologies
								myClustersAlt <- Trees[[i]]
								myClustersAlt <- .root(myClustersAlt, outgroup)
								
								t2 <- .extractClades(myClustersAlt)
								s2 <- lapply(t2, sort)
								x <- match(s1, s2)
								u <- which(!is.na(x)) # same leaves
								u <- u[!tail(duplicated(c(t2, t1[u])), length(u))] # different branching
								u <- u[!tail(duplicated(c(t, t2[x[u]])), length(u))] # untried
								rows <- ifelse(x[u] > nrow(myClustersAlt), x[u] - nrow(myClustersAlt), x[u])
								u <- u[!(myClustersAlt[, 7:8][x[u]] %in% rows)] # only keep the smallest subtrees
								
								# test alternative topologies
								while (length(u) > 0) {
									.totGrafts <- .totGrafts + 1L
									r1 <- myClusters[, 7:8][u[1]]
									r2 <- myClustersAlt[, 7:8][x[u[1]]]
									
									i1 <- .getClusters(r1, Inf, myClusters)
									i2 <- .getClusters(r2, Inf, myClustersAlt)
									i1 <- sort(i1)
									i2 <- sort(i2)
									subTree <- myClustersAlt[i2,, drop=FALSE]
									node <- subTree[, 7:8] > 0
									subTree[, 7:8][node] <- i1[match(subTree[, 7:8][node], i2)]
									myClusters[i1,] <- subTree
									
									if (method == 3) { # optimize all branch lengths
										params <- as.vector(myClusters[, 4:5])
										params <- .globalBranches(.minimize, params)
										myClusters[, 4:5] <- params
									} else {
										.minimize(myClusters[, 7:8])
									}
									
									if (.best < best - epsilon) {
										.Grafts <- .Grafts + 1L
										best <- .best
										
										if (method == 7) {
											params <- .Sankoff(myClusters[, 7:8],
												Z, # integer encoded XStringSet
												S, # substitution matrix
												weights_MP,
												scoreOnly=FALSE,
												processors=processors)[[3]]
											myClusters[, 4:5] <- params/N # changes per site
										}
										myClustersTemp <- myClusters
										
										last <- t1
										t1 <- .extractClades(myClusters)
										t <- unique(c(t, last[!(last %in% t1)]))
										s1 <- lapply(t1, sort)
									} else {
										myClusters <- myClustersAlt
										subTree <- myClustersTemp[i1,, drop=FALSE]
										node <- subTree[, 7:8] > 0
										subTree[, 7:8][node] <- i2[match(subTree[, 7:8][node], i1)]
										myClusters[i2,] <- subTree
										myClustersAlt <- myClusters
										myClusters <- myClustersTemp
										
										last <- t2
										t2 <- .extractClades(myClustersAlt)
										t <- unique(c(t, last[!(last %in% t2)]))
										s2 <- lapply(t2, sort)
									}
									x <- match(s1, s2)
									u <- which(!is.na(x)) # same leaves
									u <- u[!tail(duplicated(c(t2, t1[u])), length(u))] # different branching
									u <- u[!tail(duplicated(c(t, t2[x[u]])), length(u))] # untried
									rows <- ifelse(x[u] > nrow(myClustersAlt), x[u] - nrow(myClustersAlt), x[u])
									u <- u[!(myClustersAlt[, 7:8][x[u]] %in% rows)] # only keep the smallest subtrees
									
									if (verbose) {
										.printLine(.best)
										flush.console()
									}
								}
							}
						}
						
						if (method == 3 &&
							.Grafts > 0L &&
							((typeX != 3 && model != "JC69") ||
							(typeX == 3 && !(model %in% colnames(ProtModels)))) &&
							.best < .overall - epsilon) {
							temp <- .optimizeModel(myClusters,
								rownames(m),
								myXStringSet,
								N,
								TRUE, # scaleTree
								.rates,
								m,
								weights_ML,
								factr=min(1e10, min(absTol/.best, relTol)/.Machine$double.eps/10),
								processors=optProcessors)
							if (temp[15L] < .best - epsilon) { # improvement
								m[1,] <- temp[-1:-2]
								myClusters[, 4:5] <- temp[1L]*myClusters[, 4:5] # scale tree
								model_params <- .giveParams(m,
									rownames(m),
									.rates)
								.best <- m[1, "-LnL"]
							}
							
							if (verbose) {
								.printLine(.best)
								flush.console()
							}
						}
						
						# replace best tree
						Trees[[w]] <- myClusters
						Scores[w] <- .best
						if (method == 3)
							Models[[it]] <- list(m, model_params)
						
						if (.totShakes > 0L &&
							difftime(Sys.time(), currentTime, units="hours") > max_time[3L]) {
							doShakes <- FALSE
							doGrafts <- FALSE
							final <- TRUE
						} else if (!doShakes) {
							if (.totShakes >= max_iterations[3L] ||
								((currentScore - .best)/.best < relTol &&
								currentScore - .best < absTol)) {
								doGrafts <- FALSE
								final <- TRUE
							} else {
								doShakes <- TRUE
								if (verbose) {
									if (!interactive())
										.printLine(.best, TRUE)
									cat("\n")
								}
							}
						} else {
							if (verbose) {
								if (!interactive())
									.printLine(.best, TRUE)
								cat("\n")
							}
						}
					}
					
					currentScore <- .best
					currentClimbs <- .Climbs
					
					if (final || (!is.na(.Climbs) && doClimbs)) { # perform nearest neighbor interchanges
						if (final) {
							myClusters <- .reorderClusters(myClusters, all=TRUE)
							myClusters <- .adjustTreeHeights(myClusters)
							myClusters <- .root(myClusters, root)
						}
						out <- .NNI(myClusters)
						if (out[[2L]] > 0) {
							.Climbs <- .Climbs + 1L
							myClusters <- out[[1L]]
							if (method == 3) {
								if (out[[2L]] > 1L) { # optimize all branch lengths
									params <- as.vector(myClusters[, 4:5])
									params <- .globalBranches(.minimize, params)
									myClusters[, 4:5] <- params
								} else if (verbose) {
									.printLine(.best)
								}
							} else if (method == 7) {
								.best <- out[[4L]]
								if (verbose)
									.printLine(.best)
							}
							if (final) {
								myClusters <- .reorderClusters(myClusters, all=TRUE)
								myClusters <- .adjustTreeHeights(myClusters)
							}
						}
					}
					
					if (!is.na(.Shakes) && doShakes) {
						.overall <- best <- .best
						.Grafts <- NA
						myClustersTemp <- myClusters
						new <- TRUE
						if (verbose)
							.lastShake <- .totShakes + repsShakes
						for (i in seq_len(repsShakes)) {
							it <- it + 1L
							graft[it] <- TRUE
							.totShakes <- .totShakes + 1L
							.Shakes <- .Shakes + 1L
							.Climbs <- 0L
							.best <- Inf
							
							if (verbose) {
								cat("\n3/3. Optimizing shaken tree #",
									.totShakes,
									" of ",
									ifelse(.lastShake >= max_iterations[3L],
										max_iterations[3L],
										paste(.lastShake, "to", max_iterations[3L])),
									":\n",
									sep="")
								flush.console()
							}
							
							if (new) {
								new <- FALSE
								support <- .support(myClusters, Trees)
							}
							
							myClusters <- .rNNIs(myClusters,
								min(fracRandomNNIs, mean(support < 1)),
								1/support)
							if (method == 3) { # optimize all branch lengths
								params <- as.vector(myClusters[, 4:5])
								params <- .globalBranches(.minimize, params)
								myClusters[, 4:5] <- params
								if (verbose)
									.printLine(.best)
							}
							
							lastNNI <- Inf
							repeat {
								outTemp <- .NNI(myClusters)
								
								if (outTemp[[2L]] == 0L) {
									if (.Climbs == 0L &&
										method == 7) {
										.best <- out[[4L]]
										if (verbose)
											.printLine(.best)
									}
									break
								}
								.Climbs <- .Climbs + 1L
								myClusters <- outTemp[[1L]]
								if (method == 3 &&
									outTemp[[2L]] > 1L) { # optimize all branch lengths
									params <- as.vector(myClusters[, 4:5])
									params <- .globalBranches(.minimize, params)
									myClusters[, 4:5] <- params
								} else if (method == 7) {
									.best <- outTemp[[4L]]
									if (verbose)
										.printLine(.best)
								}
								if (outTemp[[4L]] > lastNNI - absTol)
									break
								lastNNI <- outTemp[[4L]]
							}
							
							if (method == 3 &&
								((typeX != 3 && model != "JC69") ||
								(typeX == 3 && !(model %in% colnames(ProtModels)))) &&
								.best < best - epsilon) {
								temp <- .optimizeModel(myClusters,
									rownames(m),
									myXStringSet,
									N,
									TRUE, # scaleTree
									.rates,
									m,
									weights_ML,
									factr=min(1e10, min(absTol/.best, relTol)/.Machine$double.eps/10),
									processors=optProcessors)
								if (temp[15L] < .best - epsilon) { # improvement
									m[1,] <- temp[-1:-2]
									myClusters[, 4:5] <- temp[1L]*myClusters[, 4:5] # scale tree
									model_params <- .giveParams(m,
										rownames(m),
										.rates)
									.best <- m[1, "-LnL"]
								}
								
								if (verbose) {
									.printLine(.best)
									flush.console()
								}
							}
							
							if (verbose && !interactive())
								.printLine(.best, TRUE)
							
							Scores[it] <- .best
							Trees[[it]] <- myClusters
							if (method == 3)
								Models[[it]] <- list(m, model_params)
							
							if (.best < best - epsilon) {
								myClustersTemp <- myClusters
								new <- TRUE
								best <- .best
								.overall <- .best
								out <- outTemp
							} else {
								.Shakes <- .Shakes - 1L
								myClusters <- myClustersTemp
							}
						}
						
						.best <- best
						if (.totShakes >= max_iterations[3L] ||
							((currentScore - .best)/.best < relTol &&
							currentScore - .best < absTol) ||
							.Shakes == 0L)
							doShakes <- FALSE
					}
					
					exit <- !doGrafts &&
						(((currentScore - .best)/.best < relTol &&
						currentScore - .best < absTol) ||
						((is.na(.Climbs) || currentClimbs == .Climbs) &&
						(is.na(.Shakes) || .Shakes == 0L)))
					
					if (method == 3) {
						if (!final &&
							allow &&
							((typeX != 3 && model != "JC69") ||
							(typeX == 3 && !(model %in% colnames(ProtModels)))) &&
							((it == 1 && !exit) ||
							(it > 1 && (.best < .overall || (exit && .best < .overall*fracParams))))) {
							temp <- .optimizeModel(myClusters,
								rownames(m),
								myXStringSet,
								N,
								TRUE, # scaleTree
								.rates,
								m,
								weights_ML,
								factr=min(1e10, min(absTol/.best, relTol)/.Machine$double.eps/10),
								processors=optProcessors)
							m[1,] <- temp[-1:-2]
							myClusters[, 4:5] <- temp[1L]*myClusters[, 4:5] # scale tree
							temp_params <- .giveParams(m,
								rownames(m),
								.rates)
							if (m[1, "-LnL"] < .best - epsilon) {
								.best <- m[1, "-LnL"]
								if (.best > currentScore - epsilon)
									allow <- FALSE
								if ((currentScore - .best)/.best < relTol &&
									currentScore - .best < absTol) {
									model_params <- temp_params
								} else {
									exit <- FALSE
								}
								if (verbose)
									.printLine(.best)
							}
						} else {
							temp_params <- model_params
						}
					}
					
					if (exit) {
						if (verbose && !interactive())
							.printLine(.best, TRUE)
						break # convergence
					}
					
					if (method == 3)
						model_params <- temp_params
				}
				
				if (doClimbs) { # not shaking
					Scores[it] <- .best
					Trees[[it]] <- myClusters
					Cophenetic[[it]] <- .cophenetic(myClusters)
					if (method == 3)
						Models[[it]] <- list(m, model_params)
				}
			}
			
			if (method == 7) {
				params <- .Sankoff(myClusters[, 7:8],
					Z, # integer encoded XStringSet
					S, # substitution matrix
					weights_MP,
					scoreOnly=FALSE)
				myClusters[, 4:5] <- params[[3]]/N # changes per site
				if (reconstruct && type > 1) {
					f <- function(x) {
						x <- states[x]
						x[is.na(x)] <- "."
						paste(x, collapse="")
					}
					orgXStringSet <- apply(Z, 1, f)
					states <- list(apply(params[[2]], 1, f))
				}
			}
			
			myClusters <- .Call("reclusterNJ",
				myClusters,
				cutoff[1],
				PACKAGE="DECIPHER")
			
			if (verbose) {
				.printLine(.best)
				cat("\n")
				flush.console()
			}
		} else {
			.best <- m[, "-LnL"]
		}
		
		if (verbose && method != 7) {
			params <- formatC(round(model_params, 3),
				digits=3,
				format="f")
			if (typeX == 3L) {
				cat(ifelse(all(is.na(m[c(211, 213)])),
						"",
						"\nModel parameters:"),
					ifelse(is.na(m[211]),
						"",
						paste("\nFrequency(-) =", params[211],
							"\nRate indels =", params[212])),
					ifelse(is.na(m[213]),
						"",
						paste("\nAlpha = ",
							formatC(round(m[213], 3),
								digits=3,
								format="f"),
							sep="")),
					ifelse(all(is.na(m[c(211, 213)])),
							"",
							"\n"),
					sep="")
			} else {
				cat(ifelse(all(is.na(m[c(1, 5, 11)])),
						"",
						"\nModel parameters:"),
					ifelse(is.na(m[1]),
						"",
						ifelse(grepl("T92", rownames(m), fixed=TRUE),
							paste("\nFrequency(A) = Frequency(T) = ", params[1],
								"\nFrequency(C) = Frequency(G) = ", params[2],
								sep=""),
							paste("\nFrequency(A) = ", params[1],
								"\nFrequency(C) = ", params[2],
								"\nFrequency(G) = ", params[3],
								"\nFrequency(T) = ", params[4],
								sep=""))),
					ifelse(is.na(m[5]),
						"",
						paste("\nFrequency(-) =", params[5],
							"\nRate indels =", params[11])),
					ifelse(is.na(m[6]),
						"",
						ifelse(is.na(m[7]),
							paste("\nTransition rates = ", params[6],
								"\nTransversion rates = 1",
								sep=""),
							ifelse(is.na(m[8]),
								paste("\nRate A <-> G = ", params[6],
									"\nRate C <-> T = ", params[7],
									"\nTransversion rates = 1",
									sep=""),
								paste("\nRate A <-> C = ", params[8],
									"\nRate A <-> G = ", params[6],
									"\nRate A <-> T = ", params[9],
									"\nRate C <-> G = ", params[10],
									"\nRate C <-> T = ", params[7],
									"\nRate G <-> T = 1.000",
									sep="")))),
					ifelse(is.na(m[12]),
						"",
						paste("\nAlpha = ",
							formatC(round(m[12], 3),
								digits=3,
								format="f"),
							sep="")),
					ifelse(all(is.na(m[c(1, 5, 12)])),
							"",
							"\n"),
					sep="")
			}
		}
		if (method == 3) {
			if (typeX == 3L) {
				params <- m[, 1:213]
			} else {
				params <- m[, 1:12]
			}
		} else {
			model <- NULL
			params <- NULL
		}
	} else {
		model <- NULL
		params <- NULL
		.best <- NULL
	}
	
	if (showPlot || type > 1) {
		if (method != 3 &&
			method != 7 &&
			(method == 1 ||
			root > 0))
			myClusters <- .root(myClusters, root)
		
		# create a dendrogram object
		myClustersList <- list()
		dNames <- labels(myDistMatrix)
		if (is.null(dNames)) {
			myClustersList$labels <- 1:(dim(myClusters)[1] + 1)
		} else {
			if (is(dNames, "list")) {
				myClustersList$labels <- dNames[[1]]
			} else {
				myClustersList$labels <- dNames
			}
			
			w <- which(duplicated(myClustersList$labels))
			if (length(w) > 0) {
				warning("Duplicated labels in myDistMatrix appended with index.")
				myClustersList$labels[w] <- paste(myClustersList$labels[w],
					w,
					sep="_")
			}
		}
		
		myClustersList$merge <- myClusters[, 7:8, drop=FALSE]
		myClustersList$height <- myClusters[, 6, drop=FALSE]
		myClustersList$lengths <- myClusters[, 4:5, drop=FALSE]
		if (dim > 100) {
			fontSize <- 0.6
		} else if (dim > 70) {
			fontSize <- 0.7
		} else if (dim > 40) {
			fontSize <- 0.8
		} else {
			fontSize <- 0.9
		}
		if (dim > 300) {
			leaves <- "none"
		} else {
			leaves <- "perpendicular"
		}
		
		if (method == 3) {
			if (nrow(myClusters) > 2) {
				if (out[[2L]] == 0L) { # no remaining NNIs
					scores <- out[[3L]][, 4L]
					w <- out[[3L]][, 1:2]
					LnL <- out[[4]]
				} else { # convergence despite possible NNIs
					out <- .localBranches(myClusters,
						myXStringSet,
						model_params,
						weights_ML,
						optProcessors)
					scores <- out[[2L]]
					w <- out[[1L]][, 1:2]
					LnL <- out[[4]]
				}
				
				p <- tapply(scores,
					myClusters[, 7:8][w],
					function(x)
						1/sum(1, exp(LnL - x)))
				probs <- rep(NA_real_, nrow(myClusters))
				probs[as.numeric(names(p))] <- p
			} else {
				probs <- rep(1, nrow(myClusters))
			}
		} else {
			probs <- NULL
		}
		
		if (method == 3 || method == 7) {
			support <- rep(NA_real_, nrow(myClusters))
			w <- which(myClusters[, 7:8] > 0)
			support[myClusters[, 7:8][w]] <- .support(myClusters, Trees)
		} else {
			support <- NULL
		}
		
		if (reconstruct && type > 1) {
			# perform ancestral state reconstruction
			if (method != 7) {
				states <- .Call("clusterML",
					myClusters,
					orgXStringSet,
					model_params,
					integer(),
					numeric(),
					reconstruct,
					typeX,
					weights_ML,
					optProcessors,
					PACKAGE="DECIPHER")
			}
			d <- to.dendrogram(myClustersList,
				states[[1]],
				probs,
				support)
			orgXStringSet <- unname(as.character(orgXStringSet))
			d <- rapply(d,
				function(x) {
					attr(x, "state") <- orgXStringSet[x]
					x
				},
				how="replace")
			if (method == 3)
				attr(d, "siteLnLs") <- states[[2]]
		} else {
			d <- to.dendrogram(myClustersList,
				p=probs,
				s=support)
		}
		
		# convert bifurcating tree to multifurcating
		if (collapse >= 0)
			d <- .collapse(d, collapse, dim)
		
		attr(d, "method") <- METHODS[method]
		attr(d, "model") <- model
		attr(d, "parameters") <- params
		attr(d, "score") <- .best
		
		# specify the order of clusters that
		# will match the plotted dendrogram
		orderDendrogram <- order.dendrogram(d)
		
		c <- .organizeClusters(myClusters, myClustersList$labels, orderDendrogram)
		
		if (is.finite(cutoff[1])) {
			# create a visibily different vector of colors
			cl <- colors()
			v1 <- c(117,254,73,69,152,51,26,450,503,596,652,610,563,552,97)
			r <- cl[v1]
			
			d <- .colEdge(d, dim, myClustersList$labels, r, c)
			d <- .reorder(d, dim, c)
		}
		
		# add midpoints to the tree
		d <- .applyMidpoints(d, dim)
	}
	if (type==1 || type==3) {
		dNames <- labels(myDistMatrix)
		if (is.null(dNames)) {
			dNames <- 1:(dim(myClusters)[1] + 1)
		} else {
			if (is(dNames, "list"))
				dNames <- dNames[[1]]
			
			w <- which(duplicated(dNames))
			if (length(w) > 0) {
				if (type == 1 && !showPlot)
					warning("Duplicated labels in myDistMatrix appended with index.")
				dNames[w] <- paste(dNames[w],
					w,
					sep="_")
			}
		}
		
		if (type == 1) # do not number clusters by order of tree
			c <- .organizeClustersFast(myClusters, dNames)
		
		if (length(cutoff) > 1) {
			names(c) <- paste("cluster",
				gsub("\\.", "_", cutoff[1]),
				METHODS[method],
				sep="")
			for (i in 2:length(cutoff)) {
				if (method == 1 ||
					method == 3 ||
					method == 7 ||
					root > 0) {
					myClusters <- .Call("reclusterNJ",
						myClusters,
						cutoff[i],
						PACKAGE="DECIPHER")
				} else { # ultrametric
					myClusters <- .Call("reclusterUPGMA",
						myClusters,
						cutoff[i],
						PACKAGE="DECIPHER")
				}
				x <- .organizeClustersFast(myClusters, dNames)
				if ((method == 1 ||
					method == 3 ||
					method == 7 ||
					root > 0) &&
					!ASC) # ensure clusters are subsets
					x[, 1] <- .splitClusters(x[, 1], c[, dim(c)[2]])
				names(x) <- paste("cluster",
					gsub("\\.", "_", cutoff[i]),
					METHODS[method],
					sep="")
				c <- cbind(c, x)
			}
		}
		myClusters <- c
	}
	
	if (showPlot)
		plot(d,
			horiz=FALSE,
			leaflab=leaves,	
			nodePar=list(lab.cex=fontSize, pch = NA))
	
	if (verbose) {
		time.2 <- Sys.time()
		cat("\n")
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	
	if (type==1) {
		return(myClusters)
	} else if (type==2) {
		return(d)
	} else {
		return(list(myClusters, d))
	}
}
