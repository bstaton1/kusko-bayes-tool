
##### FUNCTION TO PREPARE DATA FOR FITTING REGRESSION RELATIONSHIP #####
prepare_fit_data = function(dt, index, N) {
  # extract this day from the index for historical years
  index_sub = filter(index, date == dt)
  
  # merge in total run with 
  dat = merge(index_sub, N, by = "year")
  
  # add a q period
  dat$q = factor(ifelse(dat$year < 2008, 1, 2))
  
  dat[,c("year", "N", "q", "ccpue")]
}

##### FUNCTION TO SAMPLE FROM PREDICTED RUN SIZE DISTRIBUTION #####
sample_likelihood = function(fit, pred_ccpue, pred_q, n_mc = 500) {
  
  # extract the mean coefficients in the model
  coefs = coef(fit)
  
  # extract the covariance matrix on the coefficients
  Sigma = vcov(fit)
  
  # extract the residual standard error
  sigma = summary(fit)$sigma
  
  # turn q to be binary dummy variable
  pred_q = ifelse(pred_q == 1, 0, 1)
  
  if (n_mc > 0) {  # if doing monte carlo simulation
    # create random coefficients
    b_rand = rmvnorm(n = n_mc, mean = coefs, sigma = Sigma)
    colnames(b_rand) = names(coefs)
    
    # create random residual errors
    e_rand = rnorm(n_mc, -0.5 * sigma^2, sigma)
    
    # create the log prediction
    pred_log_N = 
      b_rand[,"(Intercept)"] + 
      b_rand[,"q2"] * pred_q +
      b_rand[,"ccpue"] * pred_ccpue +
      b_rand[,"q2:ccpue"] * pred_q * pred_ccpue
    
    # add residual error and bring to natural scale
    pred_N = exp(pred_log_N + e_rand)
  } else {  # if doing point prediction only
    pred_log_N = 
        coefs["(Intercept)"] + 
        coefs["q2"] * pred_q + 
        coefs["ccpue"] * pred_ccpue + 
        coefs["q2:ccpue"] * pred_q * pred_ccpue
    pred_N = exp(pred_log_N)
  }
  
  cbind(log_N = pred_log_N, N = pred_N)
}

##### FUNCTION TO GENERATE RANDOM HARVEST SAMPLES #####
gen_harv_samps = function(n_mc, mean, cv) {
  
  if (is.na(mean) | mean < 0) 
    mean = 0
  if (is.na(cv) | cv < 0) 
    cv = 0
  
  if (mean == 0) {
    x = rep(0, n_mc)
  } else {
    if (cv == 0) {
      x = rep(mean, n_mc)
    } else {
      x = exp(rnorm(n_mc, log(mean) - 0.5 * cv2sig(cv)^2, cv2sig(cv)))
    }
  }
  
  return(x)
}

##### FUNCTION TO GENERATE PROPOSAL FROM LOGNORMAL JUMPING DIST #####
generate.proposal = function(x, sig) {
  rlnorm(1, log(x) - 0.5 * sig^2, sig)
}

##### FUNCTION TO CALCULATE THE ASYMETRY CALCULATION #####
calc.correction = function(current, proposal, sig) {
  # find the mu for dist centered on current
  mu.current = log(current) - 0.5 * sig^2
  
  # find the mu for dist centered on proposal
  mu.proposal = log(proposal) - 0.5 * sig^2
  
  # calculate density of drawing proposal from dist centered on current
  prop.given.curr = dlnorm(x = proposal, meanlog = mu.current, sdlog = sig)
  
  # calculate density of drawing current from dist centered on current
  curr.given.prop = dlnorm(x = current, meanlog = mu.proposal, sdlog = sig)
  
  # return the output
  return(curr.given.prop/prop.given.curr)
}

##### FUNCTION TO PERFORM METROPOLIS-HASTINGS MCMC SAMPLING #####
## ARGUMENT LIST
# init: the initial value of the chain
# ni: the number of MCMC iterations (total)
# nb: the number of burnin interations
# nt: the thinning interval
# prop.sig: sd of the lognormal proposal dist
# like.fun: pdf of likelihood
# prior.fun: pdf of prior
# prior.params: vector with mean and sd of the prior distributuion

MH = function(init, ni, nb, nt, prop.sig, like.fun, prior.fun) {
  # containers to store values
  current = numeric(ni)
  proposal = numeric(ni - 1)
  accept = numeric(ni)
  prob.current = numeric(ni)
  prob.proposal = numeric(ni)
  prob.accept = numeric(ni)
  
  # step 1.) propose an initial value: is the current value for i = 1
  current[1] = init
  
  sml = 1e-300 # prevent dividing by zero in p(prop)/p(current)
  # repeat steps 2-4
  for (i in 1:(ni - 1)) {
    
    # try(progress_bar(i, ni - 1))
    
    # step 2.) calculate the posterior prob of current value with bayes' theorem
    prob.current[i] = (like.fun(current[i]) + sml) * (prior.fun(current[i]) + sml)
    
    # step 3.) generate proposal and calculate posterior prob with bayes' theorem
    proposal[i] = generate.proposal(current[i], sig = prop.sig)
    prob.proposal[i] = (like.fun(proposal[i]) + sml) * (prior.fun(proposal[i]) + sml)
    
    # step 4.) pick one to keep and one to discard
    correct = calc.correction(current = current[i], proposal = proposal[i], sig = prop.sig)
    prob.accept[i] = min(c((prob.proposal[i]/prob.current[i]) * correct, 1))
    if (is.na(prob.accept[i])) {
      current[i+1] = current[i]
    } else{
      if(prob.accept[i] > runif(1, 0, 1)) current[i+1] = proposal[i] else current[i+1] = current[i]
    }
    
    # record if the proposal was kept or if it was rejected
    accept[i] = if (current[i+1] == current[i]) 0 else 1
  }
  # acceptance rate after burnin and before thinning
  # should be between 0.2 - 0.6
  accept.rate = mean(accept[(nb+1):ni])
  
  # burnin and thin chain
  keep.iter = numeric(ni)
  if(nt > 1) keep.iter[seq(nb + 1, ni, by = nt)] = 1 else keep.iter = rep(1, ni) 
  if(nb > 0) keep.iter[1:nb] = 0 else keep.iter[1] = 1
  post.samp = current[which(keep.iter == 1)]
  
  # create and return output
  output = list(post.samp = post.samp, accept.rate = accept.rate)
  return(output)
}
