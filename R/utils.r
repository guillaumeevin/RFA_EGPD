#==============================================================================
#' EGPD.pGI, EGPD.dGI, EGPD.qGI
#'
#' First parametric family for G(v) = v^kappa: distribution, density and quantile function
#'
#' @param v probability
#' @param kappa transformation parameter greater than 0
#' @param p probability
#'
#' @return distribution, density and quantile of EGPD
#'
#' @author Guillaume Evin
#' @name functions.EGPD.GI
NULL

#' @rdname functions.EGPD.GI
EGPD.pGI = function(v,kappa) return(v^kappa)
#' @rdname functions.EGPD.GI
EGPD.dGI = function(v,kappa) return(kappa*v^(kappa-1))
#' @rdname functions.EGPD.GI
EGPD.qGI = function(p,kappa) return(p^(1/kappa))


#==============================================================================
#' dEGPD.GI, pEGPD.GI, qEGPD.GI, rEGPD.GI
#'
#' Density function, distribution function, quantile function, random generation
#' for the unified EGPD distribution
#'
#' @param x Vector of quantiles
#' @param p Vector of probabilities
#' @param n Number of observations
#' @param kappa transformation parameter greater than 0
#' @param sig Scale parameter
#' @param xi Shape parameter
#'
#' @return dEGPD.GI gives the density function, pEGPD.GI gives the distribution function, 
#' qEGPD.GI gives the quantile function, and rEGPD.GI generates random deviates.
#'
#' @author Guillaume Evin
#' @name dist.functions.EGPD.GI
NULL

#' @rdname dist.functions.EGPD.GI
dEGPD.GI = function(x,kappa,sig,xi){
  pH = Renext::pGPD(q=x,scale=sig,shape=xi)
  dH = Renext::dGPD(x=x,scale=sig,shape=xi)
  dEGPD = dH*EGPD.dGI(pH,kappa)
  return(dEGPD)
}
#' @rdname dist.functions.EGPD.GI
pEGPD.GI = function(x,kappa,sig,xi){
  pH = Renext::pGPD(q=x,scale=sig,shape=xi)
  return(EGPD.pGI(pH,kappa))
}
#' @rdname dist.functions.EGPD.GI
qEGPD.GI = function(p,kappa,sig,xi){
  qG = EGPD.qGI(p,kappa)
  qH = Renext::qGPD(p=qG,scale=sig,shape=xi)
  return(qH)
}
#' @rdname dist.functions.EGPD.GI
rEGPD.GI = function(n,kappa,sig,xi){
  u = stats::runif(n=n)
  rEGPD = qEGPD.GI(u,kappa,sig,xi)
  return(rEGPD)
}


#==============================================================================
#' EGPD.GI.mu0, EGPD.GI.mu1, EGPD.GI.mu2
#'
#' Probability Weighted Moments of order 0, 1 and 2 of the unified EGPD distribution
#'
#' @param kappa transformation parameter greater than 0
#' @param sig Scale parameter
#' @param xi Shape parameter
#'
#' @return Probability Weighted Moments
#'
#' @author Guillaume Evin
#' @name PWM.EGPD.GI
NULL

#' @rdname PWM.EGPD.GI
EGPD.GI.mu0 = function(kappa,sig,xi){
  mu0 = (sig/xi)*(kappa*beta(kappa,1-xi)-1)
  return(mu0)
}

#' @rdname PWM.EGPD.GI
EGPD.GI.mu1 = function(kappa,sig,xi){
  mu1 = (sig/xi)*(kappa*(beta(kappa,1-xi)-beta(2*kappa,1-xi))-1/2)
  return(mu1)
}

#' @rdname PWM.EGPD.GI
EGPD.GI.mu2 = function(kappa,sig,xi){
  mu2 = (sig/xi)*(kappa*(beta(kappa,1-xi)-2*beta(2*kappa,1-xi)+beta(3*kappa,1-xi))-1/3)
  return(mu2)
}


#==============================================================================
#' EGPD.GI.fPWM
#'
#' Parameter estimation of the unified EGPD distribution with the PWM method.
#' Set of equations which have to be equal to zero
#'
#' @param par vector of parameters kappa,sig,xi
#' @param pwm set of probability weighted moments of order 0, 1 and 2
#' @param xi shape parameter
#'
#' @return differences between expected and target weighted moments
#'
#' @author Guillaume Evin
# set of equations which have to be equal to zero
EGPD.GI.fPWM =  function(par,pwm,xi){
  kappa = par[1]
  sig = par[2]
  
  y = numeric(2)
  y[1] = EGPD.GI.mu0(kappa,sig,xi) - pwm[1]
  y[2] = EGPD.GI.mu1(kappa,sig,xi) - pwm[2]
  
  return(y)
}

#==============================================================================
#' EGPD.GI.fit.PWM
#'
#' Parameter estimation of the unified EGPD distribution with the PWM method.
#' Numerical solver of the system of nonlinear equations
#'
#' @param x vector of parameters kappa,sig
#' @param xi shape parameter
#'
#' @return estimated parameters kappa, sig, xi
#'
#' @author Guillaume Evin
EGPD.GI.fit.PWM = function(x,xi=0.05){
  sam.pwm = c(EnvStats::pwMoment(x,k=0),
              EnvStats::pwMoment(x,k=1))
  fit.out = nleqslv::nleqslv(c(2,stats::sd(x)),EGPD.GI.fPWM,jac=NULL,pwm=sam.pwm,xi=xi)
  fit.out$x = c(fit.out$x,xi)
  return(fit.out)
}


#==============================================================================
#' egpd_fit_main
#'
#' Main function executing the fit with auto adjustment of threshold
#'
#' @param list_scaled_POT list with one element by month. Each Element is a matrix 
#' nStation x nDate of scaled POT (precipitation exceeding a high threshold)
#' @param list_precip_byMonth  list with one element by month. Each Element is a matrix 
#' nStation x nDate of precipitation 
#' @param mat_thres matrix 12 x nStation of scaled POT thresholds
#' @param scale_byStation vector of length nStation of scale factors: mean of all positive rainfalls
#' @param thresh.loc threshold that defines zero precipitation in mm (e.g. 0.1 mm)
#' @param ROI ROI provided as a list with two levels (month/station). For each station,
#' the ROI is provided as a vector of index corr. to the stations in the vector 1:nStation
#' @param nCluster number of clusters for the parallelization
#'
#' @return list of fitted parameters with one list per month
#'
#' @author Guillaume Evin
#' 
#' @export
egpd_fit_main <- function(list_scaled_POT, 
                          list_precip_byMonth,
                          mat_thres,
                          scale_byStation,
                          thresh.loc=0.1,
                          ROI,
                          nCluster=1){
  
  # parallel computation and monitor progress
  cl <- snow::makeSOCKcluster(nCluster)
  doSNOW::registerDoSNOW(cl)
  
  # nStation
  nStation = length(scale_byStation)
  
  # define objects
  fitted_params = list() # to store the fitted parmas
  
  # fit  egpd only if the precipitation data provides 20 months of daily data (30*20)
  thNbDaysFit = 600
  
  # vector of parameters set to NA when a fit is not possible
  par.NA = list("kappa"= NA,"sigma" = NA,"xi" = NA, "cens_thres" = NA,
                "best_rmse" = NA)
  
  # prepare output
  list.par = list()
  
  # monitor progress
  pb <- txtProgressBar(min=1, max=nStation, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  
  #loop over the months
  for (im in 1:12){
    # extract the scaled POT of the current month
    scaled_POT = list_scaled_POT[[im]]
    
    # extract the observations in the current month
    precip_byMonth = list_precip_byMonth[[im]]
    
    
    # loop on stations
    st = NULL
    par.all.st = foreach(st=1:nStation, .options.snow=opts,
                         .packages=c('regEgpd','stats')) %dopar% {
      # retrieve data for this station
      neighbors = ROI[[im]][[st]]
      
      #==========================================================================
      # CASE 1: ROI identified -> calls fit_egdp_roi
      #==========================================================================
      # if at least 1 neighbor is identified, find the regional estimates
      if (length(neighbors) > 1){
        nbDays4fit = sum(!is.na(precip_byMonth[neighbors,]))
        
        # assume the scaled values of all the positive rainfall in the region follows EGPD
        if (nbDays4fit >= thNbDaysFit){
          par.fit = fit_egdp_roi(neighbors = neighbors,
                                 thresh.loc = thresh.loc,
                                 precip_byMonth = precip_byMonth,
                                 scale_byStation = scale_byStation,
                                 st = st)
        }else{
          par.fit = par.NA
        }
        
        #==========================================================================
        # CASE 2: no ROI -> calls fit.egpd.local
        #==========================================================================
      }else if (length(neighbors) ==1 & sum(!is.na(precip_byMonth[st,])) >= thNbDaysFit) {
        # station has no neighbor, fit local EGPD
        st_precip = stats::na.omit(precip_byMonth[st,]) 
        par.fit = fit.egpd.local(st_precip = st_precip,thresh.loc = thresh.loc)
      }else{
        # no ROI and not enough local data to provide a local fit
        par.fit = par.NA
      }
      
      return(par.fit)
    }
    
    # provide a data.frame for this month
    list.par[[im]] = par.all.st
  }
  
  return(list.par)
}



#==============================================================================
#' fit_egdp_roi
#'
#' automatization of ROI egdp censoring threshold. similar to above, but here we 
#' get the censoring threshold that minimizes the NRMSE
#'
#' @param neighbors vector of  the index of the identified neighbors
#' @param thresh.loc threshold that defines zero precipitation in mm (e.g. 0.1 mm)
#' @param precip_byMonth matrix nStation x nDate of the observed precip
#' @param scale_byStation vector of length nStation of scale factors: mean of all positive rainfalls
#' @param st the index of the current station
#'
#' @return regional estimate of k, sigma and xi (sigma multiplied by the current station scale)
#'
#' @author Guillaume Evin
fit_egdp_roi <- function(neighbors,thresh.loc,precip_byMonth,scale_byStation,st){
  # find the lower censoring that minimizes the nrmse
  optim_rmse = stats::optim(par = 1, fn = get_rmse_egdp_roi_MAXLK, gr = NULL,  
                            neighbors = neighbors,thresh.loc=thresh.loc,
                            precip_byMonth =  precip_byMonth, 
                            scale_byStation = scale_byStation, st = st, method = "Brent", 
                            lower = 1, upper = 9)
  best_cens_thre = optim_rmse$par
  best_rmse = optim_rmse$value
  
  # repeat the fit with this censoring threshold
  par.fit = fit_egdp_roi_MAXLK(neighbors = neighbors, precip_byMonth =  precip_byMonth, 
                               censor_value = best_cens_thre, thresh.loc=thresh.loc,
                               scale_byStation = scale_byStation, st= st)
  
  #finally returns the parameters
  return(list("kappa"= par.fit[1],"sigma" = par.fit[2],"xi" = par.fit[3],  
              "cens_thres" = best_cens_thre,
              "best_rmse" = best_rmse))
}


#==============================================================================
#' fit_egdp_roi_MAXLK
#'
#' fitting of the parameters of the EGPD using the scaled data of the ROI
#'
#' @param neighbors vector of  the index of the identified neighbors
#' @param precip_byMonth matrix nStation x nDate of the observed precip in the region
#' @param censor_value lower threshold value for censoring in mm (e.g. 2 mm)
#' @param thresh.loc threshold that defines zero precipitation in mm (e.g. 0.1 mm)
#' @param scale_byStation vector of length nStation of scale factors: mean of all positive rainfalls
#' @param st the index of the current station
#'
#' @return vector: regional estimate of k, sigma and xi (sigma multiplied by the current
#'  station scale) and the censoring threshold
#'
#' @author Guillaume Evin
fit_egdp_roi_MAXLK <- function(neighbors, precip_byMonth, 
                               censor_value, thresh.loc,
                               scale_byStation, st){
  n_cens= vec_p = nei.F = NULL
  for (nei in neighbors) {
    P.nei <- stats::na.omit(precip_byMonth[nei,]) # remove NA in ROI obs.
    if (length(P.nei)==0){
      nei.F = c(nei.F, nei)
    } else {
      # count no of obs in ROI less than censoring threshold and positive
      n_cens <- c(n_cens,length(P.nei[P.nei>0 & P.nei < censor_value])) 
      nei_scale = scale_byStation[nei] 
      
      # scale the obs in ROI by its scale
      scaled_p = P.nei/nei_scale 
      
      # select obs that are >= censoring value for ROI
      obs_hi_c = scaled_p[scaled_p>=censor_value/nei_scale] 
      vec_p <- c(vec_p,  obs_hi_c) 
    }
  }
  if(length(nei.F)>0){neighbors = neighbors[-(which(neighbors %in% nei.F))]}
  
  par.fit = EGPD.MAXLK(x = vec_p, n_cens = n_cens, censor_value=censor_value, 
                       thresh.loc = thresh.loc,
                       censor_vect = censor_value/scale_byStation[neighbors])
  
  return(c(par.fit[1], par.fit[2]*scale_byStation[st], par.fit[3], censor_value))
}


#==============================================================================
#' EGPD.MAXLK
#'
#' fitting of the parameters of the EGPD using the scaled data of the ROI
#'
#' @param x obs that are >= censoring value for ROI
#' @param n_cens no of obs in ROI less than censoring threshold and positive
#' @param censor_value lower threshold value for censoring in mm (e.g. 2 mm)
#' @param thresh.loc threshold that defines zero precipitation in mm (e.g. 0.1 mm)
#' @param censor_vect censor_value/scale[neighbors]
#'
#' @return regional estimate of k, sigma and xi
#'
#' @author Guillaume Evin
EGPD.MAXLK <- function(x, n_cens, censor_value,thresh.loc, censor_vect){
  # initial values for the optims
  inits = extRemes::fevd(x, method = "MLE", type="GP", threshold=thresh.loc, 
                         optim.args=list(method="L-BFGS-B", lower =c(0.1, 0.0001), 
                                         upper = c(100, 0.2))  )$results$par
  par_st <- c(0.9,inits)
  
  # otimizing of the likelihood function
  par.optim = stats::optim(par=par_st,fn=LK.EGPD,gr=NULL,x=x, n_cens = n_cens, 
                           censor_value= censor_value, censor_vect = censor_vect,
                           control = list(maxit = 1000), hessian = FALSE,
                           method="Nelder-Mead")$par
  
  # retourne les parametres estimes
  names(par.optim)=c("kappa","sigma","xi")
  return(par.optim)
}


#==============================================================================
#' LK.EGPD
#'
#' Likelihood of the egpd, taking into account lower censoring 
#'
#' @param par 3-parameter vector: kappa; sigma (scale); xi (shape)
#' @param x scaled data in the ROI that are >= censoring value for ROI
#' @param n_cens no of obs in ROI less than censoring threshold and positive
#' @param censor_value lower threshold value for censoring in mm (e.g. 2 mm)
#' @param censor_vect censor_value/scale[neighbors]
#'
#' @return Likelihood value
#'
#' @author Guillaume Evin
LK.EGPD <- function(par, x, n_cens, censor_value, censor_vect){
  kappa = par[1]
  sigma=par[2]
  xi=par[3]
  if(sigma>0 & abs(xi)<1 & xi > 10^(-6) & kappa >0){
    contrib.cens1 <- ifelse(censor_value > 0, 
                            sum(n_cens * log(mev::pextgp(censor_vect, 
                                                         kappa = kappa, 
                                                         sigma = sigma, 
                                                         xi = xi, type = 1))), 0)
    contrib.not.cens <- sum(mev::dextgp(x, kappa = kappa, sigma =sigma, xi = xi, 
                                        type = 1, log = TRUE))
    return(-(contrib.cens1+ contrib.not.cens))
  }else{
    return(Inf)
  }
}

#==============================================================================
#' get_rmse_egdp_roi_MAXLK
#'
#' Function linked to the automatic lower censoring threshold determination in 
#' the EGPD fit. This function computes the NRMSE corresponding to a particular 
#' choice of threshold. It is minimized by the optim function called in 
#' \code{fit_egdp_roi} function.
#'
#' @param censor_value lower threshold value for censoring in mm (e.g. 2 mm)
#' @param thresh.loc threshold that defines zero precipitation in mm (e.g. 0.1 mm)
#' @param neighbors vector of  the index of the identified neighbors
#' @param precip_byMonth matrix nStation x nDate of the observed precip in the region
#' @param scale_byStation vector of length nStation of scale factors: mean of all positive rainfalls
#' @param st the index of the current station
#'
#' @return rmse
#'
#' @author Guillaume Evin
get_rmse_egdp_roi_MAXLK = function(censor_value, thresh.loc, neighbors, 
                                   precip_byMonth,  scale_byStation, st){
  # fit the parameters
  par.fit = fit_egdp_roi_MAXLK(censor_value = censor_value, thresh.loc=thresh.loc,
                               neighbors = neighbors, precip_byMonth =  precip_byMonth, 
                               scale_byStation = scale_byStation, st)
  
  # rmse computation
  data = stats::na.omit(as.numeric(precip_byMonth[neighbors, ]))
  empirical <- data.frame(d = sort(data[data>0], decreasing = F) , 
                          probs = (1:length(data[data>0]))/(length(data[data>0])+1))
  
  # quantile of the edgpd
  mle = mev::qextgp(p=empirical$probs, kappa = par.fit[1], sigma = par.fit[2], 
                    xi = par.fit[3], type=1)
  
  # rmse
  rmse=sqrt(mean((mle - empirical$d)^2))/mean(data[data>0], na.rm=T)
  return(rmse)
}


#==============================================================================
#' fit.egpd.local
#'
#' Fit the egpd with the data from one station
#'
#' @param st_precip vector of precipitation for the station
#' @param thresh.loc threshold that defines zero precipitation in mm (e.g. 0.1 mm)
#'
#' @return list with the following fields: kappa, sigma, xi, cens_thres
#'
#' @author Guillaume Evin
fit.egpd.local  <- function(st_precip, thresh.loc){
  # matrix with two columns for the censoring threshold: 
  # - col. 1: lower bound ("censoring")
  # - col. 2: step ("rounded")
  cens_thresholds = expand.grid(censoring = c(1,2,5), rounded = c(0.1, 1, 2, 2.5))
  
  # vector of precipitation
  data = stats::na.omit(st_precip)
  
  # apply a generalized pareto as a first initialisation
  inits0 = extRemes::fevd(data, method = "MLE", type="GP", threshold=thresh.loc,
                          optim.args=list(method="L-BFGS-B", lower =c(0.1, 0.0001), 
                                          upper = c(100, 0.3))  )$results$par
  
  
  # fit an EGPD with the function fit.extgp of the mev package
  # provide a finer initialisation
  pot = data[data>thresh.loc]
  extgp <- mev::fit.extgp(data = pot, model=1, method = "mle", 
                          init =  c(0.9, inits0), censoring=c(2,Inf), rounded = 0.1,
                          confint = F,  plots = F)
  inits1 = extgp$fit$mle
  
  # empirical estimations of cdf values
  empirical <-data.frame(d = sort(pot, decreasing = F), 
                         probs = (1:length(pot))/(length(pot)+1))
  
  # find the lower censoring that minimizes the nrmse
  optim_rmse = stats::optim(par = 1, fn = get.rmse.egpd.local, gr = NULL,  
                            pot = pot,init=inits1,empirical=empirical,
                            method = "Brent", 
                            lower = 1, upper = 9)
  best_cens_thre = optim_rmse$par
  best_rmse = optim_rmse$value
  final_params <- mev::fit.extgp(data = pot, model=1, method = "mle", 
                                 init =  inits1, 
                                 censoring=c(best_cens_thre,Inf), 
                                 rounded = 0.1,
                                 confint = F,  plots = F)$fit$mle
  
  
  
  # choose the set of parameters that minises the rmse
  kappa <- final_params[1]
  sigma <- final_params[2]
  xi <- final_params[3]
  
  return(list("kappa"= kappa,"sigma" = sigma, "xi" = xi,  
              "cens_thres" = best_cens_thre,"best_rmse"=best_rmse))
}

#==============================================================================
#' get.rmse.egpd.local
#'
#'
#' Function linked to the automatic lower censoring threshold determination in 
#' the EGPD fit. This function computes the NRMSE corresponding to a particular 
#' choice of threshold.
#' 
#' @param censor_value lower threshold value for censoring in mm (e.g. 2 mm)
#' @param pot peak over threshold
#' @param init vector of EGPD parameters provided as initial values
#' @param empirical list of empirical estimates of cdf values with tho fields:
#' probs (empirical probabilies) and d (corresponding quantiles from pot)
#'
#' @return rmse
#'
#' @author Guillaume Evin
get.rmse.egpd.local = function(censor_value,pot,init,empirical){
  # fit with a different censoring
  fit.extgp.out = mev::fit.extgp(data = pot, model=1, method = "mle", init =  init ,
                                 censoring=c(censor_value,Inf), rounded = 0.1, 
                                 confint = F,  plots = F)
  par=fit.extgp.out$fit$mle
  
  # rmse computation
  mle = mev::qextgp(p=empirical$probs, prob=NA, kappa= par[1], 
                    delta=NA, sigma=par[2], xi= par[3], type=1)
  rmse = sqrt(mean((mle - empirical$d)^2))/mean(pot)
  
  return(rmse)
}


#==============================================================================
#' get_POT_matPrec_byMonth
#'
#' @param P.dates vector of dates
#' @param P.mat matrix of raw precipitation
#' @param thresh.loc threshold that defines a zero precipitation (e.g. noise)
#' @param Qseuil probability corr. to the quantile that defines the POT values
#' @param scale_byStation vector of length nStation of scale factors: mean of all positive rainfalls
#'
#' @return list with the following fields: 
#' \itemize{
#'  \item{"list_scaled_POT"}{list with one element per month: matrix nS x nP of POT values,
#'  values lower than the quantile defined by \code{Qseuil} are set to NA}
#'  \item{"list_precip_byMonth"}{list with one element per month: matrix nS x nP 
#'  of precipitation values, values lower than \code{thresh.loc} are set to 0}
#'  \item{"mat_thres"}{matrix of scaled POT thresholds}
#' }
#'
#' @author Guillaume Evin
#' 
#' @export
get_POT_matPrec_byMonth = function(P.dates,P.mat,thresh.loc=0.1,Qseuil=0.95,scale_byStation){
  # Monthly fit using a 3-month rolling window
  vec_month.series = as.numeric(format(P.dates, "%m"))
  vec.month.cycl = c(12,1:12,1)
  
  # number of stations / dates
  nS = ncol(P.mat)
  
  # initialize objects
  mat_thres = matrix(NA, 12, nS)
  list_scaled_POT = list()
  list_precip_byMonth = list()
  
  # loop on months
  for (i.m in 1:12){
    # matrix for the month
    months.i = vec.month.cycl[i.m+(0:2)]
    ind_month = vec_month.series %in% months.i
    mat_prec = t(P.mat[ind_month,])
    mat_prec[mat_prec <= thresh.loc] = 0
    # matrix of exceedances
    mat_scaled_POT = matrix(NA,nS,ncol(mat_prec))
    #for each station, find the POT and scale them
    for (i.s in 1:nS){
      vec_prec_s = mat_prec[i.s,]
      # threshold on only positive values
      vec_prec_p = vec_prec_s[vec_prec_s>thresh.loc] 
      thres_s = stats::quantile(vec_prec_p,Qseuil,na.rm=TRUE)
      mat_thres[i.m,i.s] = thres_s/scale_byStation[i.s]
      #keep only the POT and set non_extreme values to NA
      vec_prec_s[vec_prec_s<= thres_s] <- NA
      #scale the POT
      mat_scaled_POT[i.s,] = vec_prec_s/scale_byStation[i.s]
    }
    list_scaled_POT[[i.m]] = mat_scaled_POT
    list_precip_byMonth[[i.m]] = mat_prec
  }
  
  return(list(mat_thres=mat_thres,
              list_scaled_POT=list_scaled_POT,
              list_precip_byMonth=list_precip_byMonth))
}


#==============================================================================
#' get_scale_byStation
#'
#' scaling factor for the regionalization: one value per station (mean of positive
#' precipitations)
#'
#' @param P.mat matrix of raw precipitation
#' @param thresh.loc threshold that defines a zero precipitation (e.g. noise)
#'
#' @return scale_byStation vector of length nStation
#'
#' @author Guillaume Evin
#' 
#' @export
get_scale_byStation = function(P.mat,thresh.loc=0.1){
  # number of stations
  nStation = ncol(P.mat)
  
  # scale factor, mean of all positive rainfalls
  scale_byStation = vector(length = nStation)
  for(i in 1:nStation){
    P = P.mat[,i]
    scale_byStation[i] = mean(P[P>thresh.loc], na.rm=T)
  }
  
  return(scale_byStation)
}

