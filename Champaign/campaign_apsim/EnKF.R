# function to compute the analysis distribution from the forecasted and observed states 
DA <- function(Forecast, Observed, R = NULL, H, theta = 1, extraArg = NULL, ...) {
  
  # Forecast inputs
  Q <- Forecast$Q # process error
  X <- Forecast$X # states from all ensembles 
  
  # Observed inputs
  if (is.null(R)) R <- Observed$R
  Y <- Observed$Y
  
  # EnKF ---------------------------------------------
  
  # forecast prior info 
  mu.f <- as.numeric(apply(X, 2, mean, na.rm = TRUE))
  Pf <- cov(X)
  
  # MK: We might want to rethink this hack in the context of SDA with soil moisture
  # 0.1 is a pretty high variance value for the forecast
  diag(Pf)[which(diag(Pf) == 0)] <- 0.1 # hack for zero variance
  
  # remove covariance values from Pf is ensemble size is too small to allow for good estimation 
  if (nrow(X) < 30) Pf <- diag(diag(Pf))
  
  # process error
  if (!is.null(Q)) {
    Pf <- Pf + Q
  }
  
  # inflate forecast variance (just the diagonals) 
  Pf_theta <- Pf
  diag(Pf_theta) <- diag(Pf) * diag(theta)
  #diag(Pf_theta) <- diag(Pf) * theta
  
  # Kalman gain calculation with inflation factor
  K <- Pf_theta %*% t(H) %*% solve((R + H %*% Pf_theta %*% t(H)))
  
  # compute analysis distribution
  mu.a <- mu.f + K %*% (Y - H %*% mu.f)
  Pa   <- (diag(ncol(X)) - K %*% H) %*% Pf_theta
  
  return(list(
    mu.f = mu.f,
    Pf = Pf,
    mu.a = as.numeric(mu.a),
    Pa = Pa,
    Pf_theta = Pf_theta,
    H=H
  ))
}

# adjust analysis using original likelihood of forecast (not inflated) 
adj.ens<-function(Pf, X, mu.f, mu.a, Pa){
  
  vars = colnames(X)
  X = as.matrix(X)
  
  S_f  <- svd(Pf)
  L_f  <- S_f$d
  V_f  <- S_f$v
  
  # normalize
  Z <- X*0
  
  for(i in seq_len(nrow(X))){
    
    Z[i,] <- 1/sqrt(L_f) * t(V_f)%*%(X[i,]-mu.f)
    
  }
  Z[is.na(Z)]<-0
  Z[is.infinite(Z)] <- 0
  
  # analysis
  S_a  <- svd(Pa)
  
  L_a  <- S_a$d
  V_a  <- S_a$v
  
  # analysis ensemble 
  X_a <- X*0
  for(i in seq_len(nrow(X))){
    # she decomposed Pa - then it's putting it back together but with a different Z which comes from the likelihood of that ens    
    X_a[i,] <- V_a %*%diag(sqrt(L_a))%*%Z[i,] + mu.a
  }
  
  
  #if(sum(mu.a - colMeans(X_a)) > 1 | sum(mu.a - colMeans(X_a)) < -1) logger.warn('Problem with ensemble adjustment (1)')
  #if(sum(diag(Pa) - diag(cov(X_a))) > 5 | sum(diag(Pa) - diag(cov(X_a))) < -5) logger.warn('Problem with ensemble adjustment (2)')
  
  analysis <- as.matrix(X_a)
  
  # The following is a hack to ensure that soil moisture does not go below zero
  # also need to remember that soil moisture cannot go above 1, though this doesn't typically happen 
  SWS = grep('sw',vars)
  if (length(SWS > 0)){
    analysis[,SWS] <- apply(analysis[,SWS], 2,
                            function(x) {
                              #x[x <= 0] <- median(x)
                              x[x <= 0] <- 0.08 # MK: the air-dry values for layers 2, 3, and 4 are 0.089, 0.087, 0.087 
                              x[x >= 0.8] <- 0.44 # MK: the saturated values for layers 2, 3, and 4 are 0.437, 0.435, 0.435
                              return(x)
                            })
  }
  
  # Ks also cannot go below zero
  Ks = grep('ks',vars)
  if (length(Ks > 0)){
    analysis[,Ks] <- apply(analysis[,Ks], 2,
                           function(x) {
                             x[x <= 0] <- 0.1 # MK: from the "reasonable" range of values for this parameter in real soils
                             return(x)
                           })
  }
  
  return(analysis)
}

# function to apply method by Miyoshi et. al (2013) for inflation factor and R matrix 
inflate <- function(Observed, X, Xa, H, R, theta, p = 0.5){
  
  #scalar = FALSE
  
  # determine innovations for analysis and forecast
  obs.mat = matrix(rep(Observed$Y, nrow(X)), nrow = nrow(X), byrow = TRUE)
  
  # check to see if we are adjusting variables outside of obs data 
  # i.e. if we have more forecast variables than observed 
  flag = FALSE
  if (ncol(X) > ncol(obs.mat)){
    # which forecast are in obs?
    inds = which(apply(H, 2, function(col){any(col == 1)}))
    fn = ncol(X)
    
    # adjust all data to the correct dimensionos
    X <- X[,inds]
    Xa <- Xa[,inds]
    theta <- theta[inds,inds]
    H <- H[,inds]
    flag = TRUE
  }
  
  dof = obs.mat - X 
  daf = obs.mat - Xa
  
  # update R using E(daf x dof) = cov(daf, dof) + mean(daf)*mean(dof)
  Ru = cov(daf, dof) + (apply(daf, 2, mean) %*% t(apply(dof, 2, mean)))
  
  # then make sure R has no covariance values
  diags = diag(Ru)
  Ru = Ru * 0 
  diag(Ru) = diags
  
  # estimate updated inflation factor
  exp = cov(dof) + (apply(dof, 2, mean) %*% t(apply(dof, 2, mean)))
  
  # scalar
  #if (scalar){
  #  denom = H %*% (cov(X) * theta) %*% t(H) * solve(Ru)
  #  theta_update = (sum(diag(exp * solve(Ru))) - ncol(Ru)) / sum(diag(denom))
  #}
  #lse{ 
  # non-scalar
  denom = solve(H %*% (cov(X) * theta) %*% t(H))
  theta_update = (exp - Ru) * denom
  #}
  
  # smooth inflation factor using past estimated and newly estimated
  theta = ((1-p) * theta) + (p * theta_update)
  R = ((1-p) * R) + (p * Ru)
  
  if (flag){
    theta_bigger = diag(fn)
    theta_bigger[inds, inds] = theta
    theta = theta_bigger
  }
  
  # set minimum R variance values to be around accuracy ranges for sensors (at the very least)
  vars = colnames(R)
  diags = diag(R)
  SMs = grep('sw',vars)
  if (length(SMs) > 0) {
    toChange = which(diags[SMs] < 0.01^2)
    diags[SMs[toChange]] = 0.01^2
  }
  
  STs = grep('st',vars)
  if (length(STs) > 0){
    toChange = which(diags[STs] < 0.1^2)
    diags[STs[toChange]] = 0.1^2
  }
  diag(R) = diags
  
  # TO DO: we probably should add a lower limit for other variables in R 
  
  # make sure theta isn't shrinking the forecast uncertainty and remove all inflation of covariance values
  #if (scalar){
  #  if (theta < 1) theta = 1
  #}else{
  tempDiag = diag(theta)
  tempDiag[which(tempDiag < 1)] = 1
  theta = diag(length(tempDiag))
  diag(theta) = tempDiag
  #}
  
  return(list(theta = theta, R = R))
}

# data assimilation workflow with switches
data_assimilation <- function(Forecast, Observed, R, theta, miyoshi, p = 0.5, Rcov){
  
  
  if (!is.null(Forecast$X)){
    if (nrow(Forecast$X) > 0){
      
      # form H based on available data and forecast 
      H <- Observed$H %>% as.matrix()
      
      Forecast$X <- Forecast$X %>% select(Observed$Forecast.name)
      
      # first data assimilation without adjustment or later with Miyoshi algorithm
      if (miyoshi == 0){
        
        # should R covariance values be included? if no, remove them 
        if (Rcov == 0){
          Rnames = colnames(Observed$R)
          Observed$R = diag(diag(Observed$R))
          colnames(Observed$R) = Rnames
          rownames(Observed$R) = Rnames
        }
        
        DAR <- DA(Forecast, Observed, H = H)
        DAResult <- adj.ens(DAR$Pf, Forecast$X, DAR$mu.f, DAR$mu.a, DAR$Pa)
        
      }else{
        
        # in the case that this is the first day of assimilation, R and theta will be NULL
        if (is.null(R)) R <- Observed$R
        if (is.null(theta)) theta <- matrix(1,ncol(Forecast$X),ncol(Forecast$X))
        #if (is.null(theta)) theta <- 1
        
        # should R covariance values be included? if no, remove them 
        if (Rcov == 0){
          Rnames = colnames(R)
          R = diag(diag(R))
          colnames(R) = Rnames
          rownames(R) = Rnames
        }
        
        DAR <- DA(Forecast, Observed, R = R, H = H, theta = theta)
        
        DAResult <- adj.ens(DAR$Pf, Forecast$X, DAR$mu.f, DAR$mu.a, DAR$Pa)
        
        infl <- inflate(Observed = Observed, X = as.matrix(Forecast$X),
                        Xa = DAResult, H = H, R = R, theta = theta, p = p)
        
        Observed$R <- R 
        DAR$Theta <- theta
        
        # set up values for next time step
        R <- infl$R
        theta <- infl$theta
      }
      
      # gather info on SDA for this time step
      DAR$Adj <- DAResult
      DAR$Obs <- Observed
      DAR$Forecast <- Forecast
      
      return(list(Xa = DAResult, DAR = DAR, R = R, theta = theta))
    }
  }else{
    return(NULL)
  }
}

# function to get the obs.mean and obs.cov for the date when there is known data available
get_Observed <- function(today, obs.list, dateinfo){
  
  ind <- which(names(obs.list) == today)
  
  Observed <- obs.list[[ind]]
  return(Observed)
}

# function to get the forecast ensemble set up 
get_Forecast <- function(data, date){
  
  X = data %>%
    filter(day == date)%>%
    distinct(ensemble,.keep_all = TRUE) %>% 
    dplyr::select(-ensemble) 
  
  return(X)
}

# function to get the maintained ensemble set up 
#get_Maintained <- function(data, date){

#  X = data %>%
#    filter(day == date)%>%
#    distinct(ensemble,.keep_all = TRUE) %>% 
#    dplyr::select(-ensemble) 

#  return(X)
#}

# function to set up values for next time step 
nextStep <- function(Xa, Maintained){
  
  results <- colnames(Xa) %>%
    purrr::map(function(col.name) {
      #I remove the numbers from the name and find the state variables with the 
      #exact same name with no numbers - for example for sw1 this gives me all
      # sw1-sw8
      tmp <- Maintained [, grep(paste0(gsub('[[:digit:]]+', '', col.name),"$"),
                                gsub('[[:digit:]]+', '', colnames(Maintained))
      )]
      
      #original names
      orgi.name <- colnames(tmp)
      
      #those that will be replaced
      rep.col.xa <- grep(paste0(gsub('[[:digit:]]+', '', col.name),"$"),
                         gsub('[[:digit:]]+', '', colnames(Xa))
      )
      
      if (is.null(ncol(tmp))){
        tmp <- data.frame(x = tmp)
        colnames(tmp) <- col.name
      } 
      
      # Where would the rep.col.xa will go
      rep.col.main <- purrr::map_dbl(colnames(Xa)[rep.col.xa], ~grep(.x, colnames(tmp)))
        
      #replace the column from Xa with the the one from maintained
      tmp[, rep.col.main] <- Xa[, rep.col.xa]

      tmp%>%
        `colnames<-`(c(orgi.name))
      
      
    }) %>%
    setNames(colnames(Xa))
  
  
  return(results)
}