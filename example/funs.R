mseBD <- function (OM, start, sr, srRsdl = FLQuant(1, dimnames = dimnames(window(rec(OM),
    start = start))), CV = 0.3, Ftar = 0.75, Btrig = 0.75, Fmin = Ftar *
    0.1, Blim = Btrig * 1e-04, Bpct = 0.5, Fpct = 0.5, jk = FALSE,
    bounds = NULL, intYr=1)
{
	# deal with iterations
    nits <- c(OM = dims(OM)$iter, sr = dims(params(sr))$iter, rsdl = dims(srRsdl)$iter)
    if (length(unique(nits)) >= 2 & !(1 %in% nits)) ("Stop, iters not '1 or n' in OM")
    nits <- max(nits)
	#stock(OM) = propagate(stock(OM), nits)
 	if(dims(OM)$iter!=nits) stock(OM) <- propagate(stock(OM), nits)

    bd <- as(OM, "FLBioDym")
    if (!is.null(bounds)) bd@bounds <- bounds
    #bd = propagate(bd, nits)
    index(bd) <- index(bd) * rlnorm(prod(dim(index(bd))), 0, CV)

	# object to register TACs, first intYr equal to last observed
	TACrec <- FLQuant(NA, dimnames=dimnames(window(catch(bd), start=start)))
	TACrec[, ac(start+1)] <- catch(bd)[, ac(start)]

	# Go !!
    for (iYr in start:(range(OM, "maxyear") - 2)) {
        cat("===================", iYr, "===================\n")
        bd <- window(bd, end = iYr)
        # OEM ???
        index(bd)[, ac(iYr)] <- stock(OM)[, ac(iYr)] * rlnorm(prod(dim(index(bd)[, ac(iYr)])), 0, CV)
        catch(bd)[, ac(iYr)] <- computeCatch(OM)[, ac(iYr)]
		# HCR
        if (jk) {
            hv <- hcrJK(bd, Ftar, Btrig, Fmin, Blim, Fpct, Bpct)
        }
        else {
            bd <- admbBD(bd)
			intCat <- window(catch(bd), start=iYr, end=iYr+intYr)
			# TAC where caught, good place to introduce implementation error
			intCat[,ac((iYr+1):(iYr+intYr))] <- TACrec[,ac((iYr+1):(iYr+intYr))] 
			bd <- fwd(bd, catch=intCat)
            hv <- hcr(bd, FLPar(Ftar = Ftar, Btrig = Btrig, Fmin = Fmin, Blim = Blim), lag=intYr)
        }
        tac <- TAC(bd, hv)
        ctrl <- fwdControl(data.frame(year = iYr + 2, max = c(NA, 2), quantity = c("catch", "f")))
        dms <- dimnames(ctrl@trgtArray)
        dms$iter <- 1:nits
        ctrl@trgtArray <- array(NA, lapply(dms, length), dms)
        ctrl@trgtArray[1, "val", ] <- tac
        ctrl@trgtArray[2, "max", ] <- 0.35
        OM <- fwd(OM, ctrl = ctrl, sr = sr, sr.residuals = srRsdl)
    }
    return(OM)
}


