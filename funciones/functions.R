

lisread <- function( fname,quiet=TRUE )
{
  # lisread: Function to read a list of data objects from a file.
  # The initial characters "##" denote a comment line (ignored).
  # The initial characters "# " denote a variable name.
  # All other lines must contain scalars or vectors of numbers.
  # Furthermore, all rows in a matrix must contain the same number of
  # columns. Row and column vectors are not converted to matrices.
  #
  # fname  : File name.
  # quiet  : If true, shut up about reporting progress.
  # result : List object with components named in the file.

  # Original functions courtesy of Jon Schnute.
  # Modifications by A.R. Kronlund.

  lis2var <- function( x )
  {
    # lis2var: Makes global variables from the components of a list
    # x      : list object with named components.
    # result : global variables with names and contents extracted from x.

    namx <- names( x )
    nx <- length( namx )
    if (nx > 0) for (i in 1:nx)
    {
      if (namx[i] != "")
      {
        cmd <- paste( namx[i],"<<- x[[i]]" )
        eval( parse(text=cmd) )
      }
    }
    namx[namx != ""]
  }

  # numvecX functions:
  #
  # Function to convert a single string with white spaces into a numeric
  # vector. The string is parsed into separate components and converted
  # to char and numeric. A direct conversion to numeric fails.

  numvec <- function( x )
  {
    # Deprecated.
    xp <- parse( text=x,white=T )
    xc <- as.character( xp )
    as.numeric( xc )
  }

  numvec2 <- function( x )
  {
    # Patch required for S6.0 where white=T option is defunct in parse.
    # Deprecated:  xp <- parse( text=x,white=T )
    # ARK 30-Oct-03: R text connections get used up, must open/close.
    tc <- textConnection( x )
    xp <- scan( tc )
    close( tc )
    xc <- as.character( xp )
    as.numeric( xc )
  }

  numvec3 <- function( x,quiet )
  {
    # ARK 12-Jan-06: Patch to read characters because Rashmi asked nicely.
    # This is a largely untested hack, no expressed or implied warrantee.

    tmpwarn <- options( "warn" )
    options( warn=-1 )
    tc <- textConnection( x )
    xp <- scan( tc, what="character",quiet=quiet )
    close( tc )
    xc <- as.character( xp )
    if ( !all(is.na(as.numeric(xc))) )
      xc <- as.numeric( xc )

    options( tmpwarn )
    xc
  }

  #------------------------------------------------------------------#

  file <- scan( fname, what=character(), sep="\n", quiet=quiet )

  f2 <- file[ regexpr("##",file)==-1 ]           # remove comments
  nf2 <- length( f2 )                            # number of lines
  llab <- regexpr( "#",f2 )==1                   # identifies label lines
  vlab <- substring( f2[llab],3 )                # variable labels

  # ARK 30-Oct-03 R does not coerce logical to character for grep.
  ilab <- grep( "TRUE",as.character(llab) )      # label indices

  nvar <- length( vlab )                         # number of variables

  # ARK 19-Jan-10 When there is only one varaible in a file, the original
  # code does not work, namely:
  #    nrow <- c( ilab[2:nvar],nf2+1) - ilab - 1
  # returns an NA because the ilab vector is of length 1.
  #
  # Calculate the number of line for each variable.
  if ( nvar == 1 )
    nrow <- (nf2+1) - ilab - 1
  else
    nrow <- c( ilab[2:nvar],nf2+1 ) - ilab - 1

  zout <- list( NULL )

  for ( i in 1:nvar )
  {
    i1 <- ilab[i] + 1                            # line of first var element
    i2 <- i1 + nrow[i] - 1                       # line of last  var element
    zstr <- paste(f2[i1:i2],collapse=" ")
#    zvec <- numvec2(zstr)                        # numeric vector
    zvec <- numvec3(zstr,quiet)                  # numeric or character vector

    nz <- length(zvec)
    zrow <- nrow[i]
    zcol <- nz / zrow                            # dimensions
    if ( (zrow>1) & (zcol>1) )                   # a true matrix
    {
      zvec <- matrix( zvec,nrow=zrow,ncol=zcol,byrow=T )
#      print( vlab[i] )
#      print( zvec )
#      scan()
    }

    zout[[i]] <- zvec
    if ( !quiet )
      cat( "vlab = ", vlab[i], "\n" )
  }
  names(zout) <- vlab
  zout
}




writeData <-function (name, L, append = FALSE) 
{
	n <- nchar(name)
	file_name <- 
		
		#if (tools::file_ext(name) == ".dat") {
		name
#	}
#	else paste(name, "dat", sep = ".")
	#cat("# \"", file_name, "\" produced by dat_write() from R2admb ", 
	#	date(), "\n", file = file_name, sep = "", append = append)
	for (i in 1:length(L)) {
		x <- L[[i]]
		dc <- data.class(x)
		if (dc == "numeric") {
			cat("#", names(L)[i], "\n", L[[i]], "\n\n", file = file_name, 
				append = TRUE)
		}
		else {
			if (dc == "matrix") {
				cat("#", names(L)[i], "\n", file = file_name, 
					append = TRUE)
				write.table(L[[i]], , col.names = FALSE, row.names = FALSE, 
					quote = FALSE, file = file_name, append = TRUE)
				cat("\n", file = file_name, append = TRUE)
			}
			else {
				stop(paste("can't handle data type '", dc, "' (variable ", 
					names(L)[i], ")", sep = ""))
			}
		}
	}
}



reptoRlist <- function(fn){
  ifile <- scan(fn,what="character",flush=T,blank.lines.skip=F,quiet=T)
  idx <- sapply(as.double(ifile),is.na)
  vnam <- ifile[idx] #list names
  nv <- length(vnam) #number of objects
  A <- list()
  ir <- 0
  for(i in 1:nv){
    ir <- match(vnam[i],ifile)
    if(i!=nv){
      irr <- match(vnam[i+1],ifile)
    }else{
      irr <- length(ifile)+1 #next row
    }
    dum <- NA
    if(irr-ir==2){
      dum <- as.double(scan(fn,skip=ir,nlines=1,quiet=T,what=""))
    }
    if(irr-ir>2){
      dum <- as.matrix(read.table(fn,skip=ir,nrow=irr-ir-1,fill=T))
    }
    if(is.numeric(dum)) #Logical test to ensure dealing with numbers
      {
        A[[vnam[i]]] <- dum
      }
  }
  return(A)
}

run_admb<-function (fn, verbose = FALSE, mcmc = FALSE, mcmc.opts = mcmc.control(), 
	profile = FALSE, extra.args = "", admb_errors = c("stop", 
		"warn", "ignore")) 
{
	admb_errors <- match.arg(admb_errors)
	args <- ""
	if (mcmc) {
		if (is.null(mcmc.opts$mcmcpars)) 
			stop("you must specify at least one parameter in 'mcmc.opts$mcmcpars' (see ?mcmc.control)")
		args <- paste(args, mcmc.args(mcmc.opts))
	}
	if (profile) 
		args <- paste(args, "-lprof")
	if (!missing(extra.args)) {
		args <- paste(args, extra.args)
	}
	if (verbose) 
		cat("running compiled executable with args: '", args, 
			"'...\n")
	outfn <- paste(fn, "out", sep = ".")
	if (.Platform$OS.type == "windows") {
		cmdname <- paste(fn, ".exe", sep = "")
		shellcmd <- shell
	}
	else {
		cmdname <- paste("./", fn, sep = "")
		shellcmd <- system
	}
	if (!file.exists(cmdname)) 
		stop("executable ", cmdname, " not found: did you forget to compile it?")
	res <- shellcmd(paste(cmdname, args, ">", outfn), intern = TRUE)
	outfile <- readLines(paste(fn, ".out", sep = ""))
	if (mcmc) {
		mcinfofile <- file(paste(fn, "mcinfo", sep = "."), "w")
		mctab <- unlist(mapply(function(x, y) {
			c(paste("# ", x), if (is.null(y)) "" else paste(y, 
				collapse = " "))
		}, names(mcmc.opts), mcmc.opts))
		writeLines(mctab, paste(fn, "mcinfo", sep = "."))
	}
	if (verbose) {
		cat("Run output:\n", res, "\n", sep = "\n")
		cat(outfile, "\n", sep = "\n")
	}
	if (length(grep("^Error", outfile) > 0)) {
		runerrmsg <- "errors detected in ADMB run: run with verbose=TRUE to view"
		if (admb_errors == "stop") 
			stop(runerrmsg)
		else if (admb_errors == "warn") 
			warning(runerrmsg)
	}
	invisible(res)
}