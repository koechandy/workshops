
# source("S:/etude40/4093_formation_ouvrage_SN_Corte/contenu_formation_ouvrage/Pratical/function/lexis.R", echo=T)

# see also http://staff.pubhealth.ku.dk/~bxc/Lexis/
#-----------------------------------------------


lexis <-function(entry = 0, exit, fail, origin = 0, scale = 1, breaks,
  include = NULL, data = NULL)
{
# If there is a dataframe argument, attach it so that arguments are matched
	if(!is.null(data)){
    attach(data, 2)
    on.exit(detach(pos=2))
	}
  a.en <- attributes(entry)
	a.ex <- attributes(exit)
	a.fa <- attributes(fail)
	a.en$names <- NULL
	a.ex$names <- NULL
	a.fa$names <- NULL
	nexit <- length(exit)
	if(length(entry) == 1)
		entry <- rep(entry, nexit)
	if(length(entry) != nexit)
		stop("wrong argument length (entry)")
	if(length(fail) != nexit)
		stop("wrong argument length (fail)")
	if(any(exit < entry)) stop("exit time < entry time")
	# Is there more than one time scale?
	if(!is.list(origin)) {
		if(length(origin) != 1 && length(origin) != nexit)
			stop("illegal argument (entry)")
		if(!is.numeric(breaks))
			stop("illegal argument (breaks)")
		if(!is.numeric(scale) || length(scale) != 1)
			stop("illegal argument (scale)")
		origin <- list(.time = origin)
		breaks <- list(.time = breaks)
		names(scale) <- ".time"
	}
	if(!is.list(breaks)) stop("breaks argument should be a list")
	# Find time scale names in breaks argument and check, otherwise assign
	tnames <- names(origin)
	nt <- length(tnames)
	bnames <- names(breaks)
	if(is.null(bnames)) {
		if(length(breaks) != nt)
			stop("wrong argument length (breaks)")
		names(breaks) <- tnames
	}
	else {
		if(any(sort(bnames) != sort(tnames)))
			stop("breaks and origin names don't match")
	}
	bnames <- names(scale)
	if(is.null(bnames)) {
		if(length(scale) == 1 && nt != 1)
			scale <- rep(scale, nt)
		names(scale) <- tnames
	}
# Covariates
	covexp <- substitute(include)
	if(mode(covexp) == "name") {
		include <- list(include)
		names(include) <- as.character(covexp)
	}
	else if(mode(covexp) == "call" && covexp[[1]] == "list") {
		nc <- length(covexp)
		cnames <- names(covexp)
		if(!is.null(cnames)) {
			for(i in 2:nc) {
				if(cnames[i] == "") {
				  if(mode(covexp[[i]]) == "name")
				    cnames[i] <- as.character(covexp[[i]])
				  else stop(
				      "illegal argument (include)")
				}
			}
		}
		else {
			for(i in 2:nc) {
				if(mode(covexp[[i]]) == "name")
				  cnames[i] <- as.character(covexp[[i]])
				else stop("illegal argument (include)")
			}
			cnames[1] <- ""
		}
		names(covexp) <- cnames
		include <- eval(covexp)
	}
	else {
		stop("illegal argument (include)")
	}
	nc <- length(include)
	cnames <- names(include)
	for(i in 1:nc) {
		if(length(include[[i]]) != nexit) {
			stop("incorrect length for included variable")
		}
	}
# List to hold level names
	lnames <- vector("list", nt)
	# Exclude follow-up before Lexis diagram begins
	for(s in tnames) {
		stt <- scale[s] * as.numeric(breaks[[s]][1]) + origin[[s]]
		entry <- ifelse(stt > entry, stt, entry)
	}
# Exclude follow-up after Lexis diagram ends
	for(s in tnames) {
		br <- breaks[[s]]
		lastbr <- as.numeric(br[length(br)])
		stp <- scale[s] * lastbr + origin[[s]]
		fail <- ifelse(stp < exit, 0, fail)
		exit <- ifelse(stp < exit, stp, exit)
	}
# Matrix to hold cell index
	cell <- matrix(nrow = nexit, ncol = nt, dimnames = list(NULL,	tnames))
# Calculate expanded vector length
	ni <- rep(1, nexit)
	for(s in tnames) {
		t <- as.numeric(entry - origin[[s]])/scale[s]
		first <- cut(t, breaks[[s]], include.lowest = T)
		cell[, s] <- first
		t <- as.numeric(exit - origin[[s]])/scale[s]
		last <- cut(t, breaks[[s]], include.lowest = T)
    ni <- ni + as.integer(last) - as.integer(first)
		lnames[[s]] <- levels(first)
	}
# Identify individuals with follow-up entirely outside the diagram
	surv <- exit >= entry
	ni <- ifelse(surv, ni, 0)
	nint <- sum(ni)
# Cycle thru cells of Lexis diagram
# At each cycle, increment entry times
# Build up expanded arrays
	x.en <- vector(mode = mode(entry), nint)
	x.ex <- vector(mode = mode(exit), nint)
	x.fa <- vector(mode = mode(fail), nint)
	attributes(x.en) <- a.en
	attributes(x.ex) <- a.ex
	attributes(x.fa) <- a.fa
	x.in <- numeric(nint)
	x.ce <- matrix(nrow = nint, ncol = nt, dimnames = list(NULL, tnames))
	no <- 0
	while(any(surv)) {
		ex <- exit
		for(s in tnames) {
			hi <- scale[s] * as.numeric(breaks[[s]][1 +	cell[, s]]) + origin[[s]]
			ex <- ifelse(hi < ex, hi, ex)
		}
		st <- no + 1
		no <- no + sum(surv)
		x.en[st:no] <- entry[surv]
		x.ex[st:no] <- ex[surv]
		x.fa[st:no] <- ifelse(ex == exit, fail, 0)[surv]
		x.in[st:no] <- (1:nexit)[surv]
		x.ce[st:no,  ] <- cell[surv,  ]
		surv <- ex < exit
		entry <- ex
		for(s in tnames) {
			hi <- scale[s] * as.numeric(breaks[[s]][1 +	cell[, s]]) + origin[[s]]
			cell[, s] <- cell[, s] + (ex == hi & surv)
		}
	}
# Build result list
	x.or <- ifelse(is.na(x.in), 0, order(x.in))
	x.in <- x.in[x.or]
	res <- vector("list", 4 + nt + nc)
#	names(res) <- c(".expand", ".entry", ".exit", ".fail", tnames, cnames)
	names(res) <- c("Expand", "Entry", "Exit", "Fail", tnames, cnames)
	res[[1]] <- x.in
	res[[2]] <- x.en[x.or]
	res[[3]] <- x.ex[x.or]
	res[[4]] <- x.fa[x.or]
	for(s in tnames) {
		res[[s]] <- factor(x.ce[x.or, s], levels =
                  1:length(lnames[[s]]), labels = lnames[[s]])
	}
	for(s in cnames) {
		res[[s]] <- include[[s]][x.in]
	}
	as.data.frame(res)
}

