
##### FUNCTION TO CREATE DISTRIBUTION PLOT #####
create_dist_plot = function(out) {
  if (!is.null(out$EstOut$boot_out)) {
    
    ### calculate quantities ###
    # bin dimensions
    # bin_dims = get_bin_dim(0, 8e5, 160)
    bin_dims = get_bin_dim(0, 5e6, 1001)
    
    # extract non-mcmc samples and calculate bin probabilities
    boot_out = out$EstOut$boot_out
    like_samp = boot_out[,"N"]
    prior_samp = out$EstOut$prior_samp
    
    prior_probs = get_bin_probs(x = prior_samp, bin_dim = bin_dims)
    like_probs = get_bin_probs(x = like_samp, bin_dim = bin_dims)
    
    if (out$EstOut$did_mcmc) {
      post_samp = as.numeric(out$EstOut$post_samp)
      post_probs = get_bin_probs(post_samp, bin_dims)
      ymax = max(c(prior_probs, like_probs[1:(length(like_probs) - 1)], post_probs))
    } else {
      post_samp = NULL
      ymax = max(c(prior_probs, like_probs[1:(length(like_probs) - 1)]))
    }
    x_var = bin_dims$mp
    n = bin_dims$n
    
    ### make plot ###
    par(cex.lab = 1.5, cex.axis = 1.4, xaxs = "i", yaxs = "i", mar = c(5,3,1,1))
    plot(1, 1 , type = "n",
         col = "black", xlim = c(0, 400000), ylim = c(0, ymax * 1.05), lwd = 2,
         xaxt = "n", yaxt = "n", ylab = "",
         xlab = "Run Size (1000s)")
    mtext(side = 2, line = 1, "Relative Probability", cex = 1.5)
    
    # text(200000, 0.5, labels = ymax, cex = 1.5)
    
    polygon(x = c(x_var, rev(x_var)), y = c(prior_probs, rep(0, n)), col = alpha("red", 0.5), border = NA)
    polygon(x = c(x_var, rev(x_var)), y = c(like_probs, rep(0, n)), col = alpha("blue", 0.5), border = NA)
    
    # if (!is.null(post_samp)) polygon(x = c(x_var, rev(x_var)), y = c(post_probs, rep(0, n)), col = alpha("yellow", 0.5), border = NA)
    # if (!is.null(post_samp)) lines(post_probs ~ bin_mp, col = "black", lwd = 3)
    
    if (!is.null(post_samp)) {
      polygon(x = c(x_var, rev(x_var)), y = c(post_probs, rep(0, n)), col = alpha("yellow", 0.5), border = NA)
      lines(post_probs ~ x_var, col = "black", lwd = 3)
      legend("topright", legend = c("Forecast", "BTF Only", "Updated"),
             pch = 22, col = "black", pt.cex = 2.5,
             pt.bg = c(alpha("red", 0.5), alpha("blue", 0.5), alpha("yellow", 0.5)), bty = "n", cex = 1.2)
    } else {
      legend("topright", legend = c("Forecast", "BTF Only"),
             pch = 22, col = "black", pt.cex = 2.5,
             pt.bg = c(alpha("red", 0.5), alpha("blue", 0.5)), bty = "n", cex = 1.2)
      
    }
    lines(like_probs ~ x_var, col = "black", lwd = 3)
    lines(prior_probs ~ x_var, col = "black", lwd = 3)
    axis(side = 1, at = seq(0, 400000, 50000), labels = seq(0, 400, 50), lwd = 4)
    box(lwd = 4)
  } else {
    par(cex.lab = 1.5, cex.axis = 1.4, xaxs = "i", yaxs = "i", mar = c(5,3,1,1))
    plot(1,1, type = "n",
         col = "black", xlim = c(0, 400000), ylim = c(0, 1), lwd = 2,
         xaxt = "n", yaxt = "n", ylab = "",
         xlab = "Run Size (1000s)")
    mtext(side = 2, line = 1, "Relative Probability", cex = 1.5)
    axis(side = 1, at = seq(0, 400000, 50000), labels = seq(0, 400, 50), lwd = 4)
    box(lwd = 4)
  }
}

##### FUNCTION TO CREATE ESTIMATE TABLE #####
create_est_table = function(out) {
  pretty_prior = prettyNum(round(summ(out$EstOut$prior_samp), -3), big.mark = ",", scientific = F)
  pretty_like = prettyNum(round(summ(out$EstOut$boot_out[,"N"]), -3), big.mark = ",", scientific = F)
  pretty_p = paste(round(summ(out$EstOut$boot_out[,"p"]), 2) * 100, "%", sep = "")
  pretty_eos = prettyNum(round(summ(out$EstOut$boot_out[,"eos"])), big.mark = ",", scientific = F)
  rnames = names(summ(rnorm(10)))
  rnames[1:2] = c("Mean", "SD")
  if (!is.null(out$EstOut$post_samp)) {
    cnames = c("Statistic", "Forecast", "BTF Only", "Updated", "% Complete", "EOS")
    pretty_post = prettyNum(round(summ(as.numeric(out$EstOut$post_samp)), -3), big.mark = ",", scientific = F)
    tab = data.frame(rnames, pretty_prior, pretty_like, pretty_post, pretty_p, pretty_eos)
    colnames(tab) = cnames
  } else {
    cnames = c("Statistic", "Forecast", "BTF Only", "% Complete", "EOS")
    tab = data.frame(rnames, pretty_prior, pretty_like, pretty_p, pretty_eos)
    colnames(tab) = cnames
  }
  # tab[,!(colnames(tab) %in% c("mean", "sd"))]
  tab
}
##### FUNCTION TO CREATE PLCY PLOT (P S< P.STAR) #####
create_plcy_plot = function(out, input) {
  if (!is.null(out$PolicyOut$know.ind) & length(out$PolicyOut$know.ind) > 0 & !is.null(out$EstOut)) {
    par(xaxs = "i", yaxs = "i", mar = c(4,4,2,2), cex.axis = 1.3, cex.lab = 1.3)
    plot(1,1, type = "n", xlim = range(out$PolicyOut$Hcand), ylim = c(0,1), 
         xlab = "Additional Harvest Candidates (1000s)", xaxt = "n", yaxt = "n",
         ylab = out$PolicyOut$ylab, las = 1)
    
    # windows()
    # plot(1,1, type = "n", ylim = c(0,1))
    usr = par("usr")
    
    if (out$PolicyOut$direction == "less") {
      rect(xleft = usr[1], ybottom = usr[3], xright = usr[2], ytop = out$PolicyOut$p.star, col = "grey90", border = NA)
    } else {
      rect(xleft = usr[1], ybottom = out$PolicyOut$p.star, xright = usr[2], ytop = usr[4], col = "grey90", border = NA)
    }
    
    abline(h = out$PolicyOut$p.star, col = "grey", lwd = 2, lty = 1)
    
    if (1 %in% out$PolicyOut$know.ind) {
      lines(out$PolicyOut$risk_p_prior ~ out$PolicyOut$Hcand, col = "red", lwd = 4)
      segments(out$PolicyOut$Htarget_prior, usr[3], out$PolicyOut$Htarget_prior, 
               out$PolicyOut$p.star, col = "red", lwd = 4, lty = 2)
    }
    
    if (2 %in% out$PolicyOut$know.ind) {
      lines(out$PolicyOut$risk_p_like ~ out$PolicyOut$Hcand, col = "blue", lwd = 4)
      segments(out$PolicyOut$Htarget_like, usr[3], out$PolicyOut$Htarget_like, 
               out$PolicyOut$p.star, col = "blue", lwd = 4, lty = 2)
    }
    
    if (3 %in% out$PolicyOut$know.ind & !is.null(out$PolicyOut$risk_p_post)) {
      lines(out$PolicyOut$risk_p_post ~ out$PolicyOut$Hcand, col = "black", lwd = 4)
      segments(out$PolicyOut$Htarget_post, usr[3], out$PolicyOut$Htarget_post, 
               out$PolicyOut$p.star, col = "black", lwd = 4, lty = 2)
    } 
    
    at.x = axisTicks(usr = usr[1:2], log = F)
    axis(side = 1, at = at.x, labels = at.x/1000, lwd = 4)
    axis(side = 2, lwd = 4, las = 2)
    mod.names = c("Forecast", "BTF Only", "Updated")
    mod.cols = c("red", "blue", "black")
    legend(ifelse(out$PolicyOut$direction == "less", "topleft", "topright"),
           legend = mod.names[out$PolicyOut$know.ind],
           col = mod.cols[out$PolicyOut$know.ind], bty = "n", pch = 15, pt.cex = 2.5, cex = 1.3)
    box(lwd = 4)
  } else {
    if (input$risk_direction == 1) {
      ylab = paste("Probability of Escapement Less than", prettyNum(input$Sobj, big.mark = ",", scientific = F))
    } else {
      ylab = paste("Probability of Escapement Greater than", prettyNum(input$Sobj, big.mark = ",", scientific = F))
    }
    # par(xaxs = "i", yaxs = "i", mar = c(4,4,2,2), cex.axis = 1.3, cex.lab = 1.3)
    # plot(1,1, type = "n", xlim = input$Hcand_lim, ylim = c(0,1), 
    #      xlab = "Additional Harvest Candidates (1000s)", xaxt = "n",
    #      ylab = ylab,
    #      las = 1)
    # at.x = axisTicks(usr = par("usr")[1:2], log = F)
    # axis(side = 1, at = at.x, labels = at.x/1000, lwd = 4)
    # axis(side = 2, lwd = 4, las = 2)
    # 
    # box(lwd = 4)
    NULL
  }
}



##### FUNCTION TO CREATE LEARNING PLOT #####
create_learn_plot = function(out) {
  par(mar = c(4,5,2,2))
  
  dat = out$LearnOut$plot_dat$dat
  dat = dat[order(dat$day, decreasing = T),]
  est = unique(dat$source)
  x.025 = dat[dat$stat == "2.5%","day"]; y.025 = dat[dat$stat == "2.5%","value"]
  x.10 = dat[dat$stat == "10%","day"]; y.10 = dat[dat$stat == "10%","value"]
  x.25 = dat[dat$stat == "25%","day"]; y.25 = dat[dat$stat == "25%","value"]
  x.50 = dat[dat$stat == "50%","day"]; y.50 = dat[dat$stat == "50%","value"]
  x.mean = dat[dat$stat == "mean","day"]; y.mean = dat[dat$stat == "mean","value"]
  x.75 = dat[dat$stat == "75%","day"]; y.75 = dat[dat$stat == "75%","value"]
  x.90 = dat[dat$stat == "90%","day"]; y.90 = dat[dat$stat == "90%","value"]
  x.975 = dat[dat$stat == "97.5%","day"]; y.975 = dat[dat$stat == "97.5%","value"]
  
  plot(1,1,type = "n",
       xlim = out$LearnOut$plot_dat$xlim,
       ylim = out$LearnOut$plot_dat$ylim, axes = F, ann = F
  )
  
  polygon(x = c(x.025, rev(x.975)), y = c(y.025, rev(y.975)), col = "grey90", border = NA)
  polygon(x = c(x.10, rev(x.90)), y = c(y.10, rev(y.90)), col = "grey60", border = NA)
  polygon(x = c(x.25, rev(x.75)), y = c(y.25, rev(y.75)), col = "grey40", border = NA)
  lines(y.025 ~ x.025, col = "grey80", lwd = 2)
  lines(y.10 ~ x.10, col = "grey50", lwd = 2)
  lines(y.25 ~ x.25, col = "grey10", lwd = 2)
  lines(y.75 ~ x.75, col = "grey10", lwd = 2)
  lines(y.90 ~ x.90, col = "grey50", lwd = 2)
  lines(y.975 ~ x.975, col = "grey80", lwd = 2)
  lines(y.50 ~ x.50, type = "o", lwd = 4, pch = 16, cex = 2)
  
  pd = seq(out$LearnOut$plot_dat$xlim[1], out$LearnOut$plot_dat$xlim[2], 1)
  npd = length(pd)
  inc = ifelse(npd <= 7, 1, 
               ifelse(npd <= 10, 2,
                      ifelse(npd < 20, 3, 5)))
  
  at.x = pd[seq(1,npd, inc)]
  lab.x = dates$date[dates$day %in% at.x]
  
  axis(side = 1, at = at.x, labels = lab.x, las = 1, lwd = 4, cex.axis = 1.4)
  
  if (est %in% c("prior", "likelihood", "posterior")) {
    axis(side = 2, at = seq(0, 800000, 50000), labels = seq(0, 800, 50), las = 2, cex.axis = 1.4, lwd = 4)
  }
  if (est == "p") {
    axis(side = 2, at = seq(0, 1, 0.1), labels = seq(0, 1, 0.1), las = 2, cex.axis = 1.4, lwd = 4)
  }
  if (est == "eos") {
    axis(side = 2, at = seq(0, 2000, 200), labels = seq(0, 2000, 200), las = 2, cex.axis = 1.4, lwd = 4)
  }
  
  mtext(side = 2, line = 3.8, out$LearnOut$plot_dat$ylab, cex = 1.4)
  box(lwd = 4)
  
  if(out$LearnOut$plot_dat$legend_loc != "none") {
    legend(out$LearnOut$plot_dat$legend_loc,
           legend = c("Median", "Central 50%", "Central 80%", "Central 95%"),
           pch = c(16, 22, 22, 22), lty = c(1, NA, NA, NA), lwd = c(2, NA, NA, NA),
           pt.bg = c(NA, "grey40", "grey60", "grey90"), pt.cex = c(1.5, 3, 3, 3), bty = "n", cex = 1.2,
           col = c("black", "grey10", "grey50", "grey80"))
  }
}

##### FUNCTION TO CREATE PLCY TABLE 2 #####
create_plcy2_table = function(out) {
  prettyES1 = prettyNum(round(out$PlcyOut2$ES1,-3), big.mark = ",", scientific = F)
  prettyES2 = prettyNum(round(out$PlcyOut2$ES2,-3), big.mark = ",", scientific = F)
  prettyES3 = prettyNum(round(out$PlcyOut2$ES3,-3), big.mark = ",", scientific = F)
  pretty_p1 = round(out$PlcyOut2$p1, 2); pretty_p1 = ifelse(pretty_p1 < 0.01, "<0.01", pretty_p1)
  pretty_p2 = round(out$PlcyOut2$p2, 2); pretty_p2 = ifelse(pretty_p2 < 0.01, "<0.01", pretty_p2)
  pretty_p3 = round(out$PlcyOut2$p3, 2); pretty_p3 = ifelse(pretty_p3 < 0.01, "<0.01", pretty_p3)
  # pretty_odds = round(out$PlcyOut2$odds_ratio,2)
  
  prettyH1 = prettyNum(out$PlcyOut2$Hobj1, big.mark = ",", scientific = F)
  prettyH2 = prettyNum(out$PlcyOut2$Hobj2, big.mark = ",", scientific = F)
  prettyH3 = prettyNum(out$PlcyOut2$Hobj3, big.mark = ",", scientific = F)
  
  q = c("Expected Escapement", "Pr(S < 65,000)", "Pr(65,000 < S < 120,000)", "Pr(S > 120,000)", "Pr(S < 95,000)", "Pr(S < 110,000)")
  df = data.frame(q, c(prettyES1, pretty_p1), c(prettyES2, pretty_p2), c(prettyES3, pretty_p3))#, c("---", pretty_odds))
  colnames(df) = c("Quantity", paste("H =", prettyH1), paste("H =", prettyH2), paste("H =", prettyH3))#, "H2:H1 Odds Ratio")
  
  # keep = ifelse(prettyH3 == prettyH2 | prettyH3 == prettyH1, c(1,2), c(1,2,3))
  
  # df[,!duplicated(c(out$PlcyOut2$Hobj1, out$PlcyOut2$Hobj2, out$PlcyOut2$Hobj3))]
  if(!is.null(out$EstOut)) {
    df
  } else {
    NULL
  }
}

##### FUNCTION TO PRETTIFY NUMBERS #####
prettify = function(x, rnd = NULL) {
  if (!is.null(rnd)) x = round(x, rnd)
  out = prettyNum(x = x, big.mark = ",", scientific = F)
  
  return(out)
}

##### FUNCITON TO CREATE POLICY TABLE 1 #####
create_plcy_table = function(out) {
  if (1 %in% out$PolicyOut$know.ind) {
    pretty_H_prior = prettify(out$PolicyOut$Htarget_prior, -3)
    pretty_ES_prior = prettify(out$PolicyOut$ES_prior,-3)
    pretty_less65_prior = round(out$PolicyOut$Pcrit_prior["less65"], 2)
    pretty_less95_prior = round(out$PolicyOut$Pcrit_prior["less95"], 2)
    pretty_less120_prior = round(out$PolicyOut$Pcrit_prior["less120"], 2)
    pretty_less110_prior = round(out$PolicyOut$Pcrit_prior["less110"], 2)
    pretty_in_EG_prior = pretty_less120_prior - pretty_less65_prior
    pretty_above120_prior = 1 - pretty_less120_prior
  } else {
    pretty_H_prior = NA
    pretty_ES_prior = NA
    pretty_less65_prior = NA
    pretty_less95_prior = NA
    pretty_less120_prior = NA
    pretty_less110_prior = NA
    pretty_in_EG_prior = NA
    pretty_above120_prior = NA
  }
  
  if (2 %in% out$PolicyOut$know.ind) {
    pretty_H_like = prettify(out$PolicyOut$Htarget_like, -3)
    pretty_ES_like = prettify(out$PolicyOut$ES_like,-3)
    pretty_less65_like = round(out$PolicyOut$Pcrit_like["less65"], 2)
    pretty_less95_like = round(out$PolicyOut$Pcrit_like["less95"], 2)
    pretty_less120_like = round(out$PolicyOut$Pcrit_like["less120"], 2)
    pretty_less110_like = round(out$PolicyOut$Pcrit_like["less110"], 2)
    pretty_in_EG_like = pretty_less120_like - pretty_less65_like
    pretty_above120_like = 1 - pretty_less120_like
  } else {
    pretty_H_like = NA
    pretty_ES_like = NA
    pretty_less65_like = NA
    pretty_less95_like = NA
    pretty_less120_like = NA
    pretty_less110_like = NA
    pretty_in_EG_like = NA
    pretty_above120_like = NA
  }
  
  if (3 %in% out$PolicyOut$know.ind) {
    pretty_H_post = prettify(out$PolicyOut$Htarget_post, -3)
    pretty_ES_post = prettify(out$PolicyOut$ES_post,-3)
    pretty_less65_post = round(out$PolicyOut$Pcrit_post["less65"], 2)
    pretty_less95_post = round(out$PolicyOut$Pcrit_post["less95"], 2)
    pretty_less120_post = round(out$PolicyOut$Pcrit_post["less120"], 2)
    pretty_less110_post = round(out$PolicyOut$Pcrit_post["less110"], 2)
    pretty_in_EG_post = pretty_less120_post - pretty_less65_post
    pretty_above120_post = 1 - pretty_less120_post
  } else {
    pretty_H_post = NA
    pretty_ES_post = NA
    pretty_less65_post = NA
    pretty_less95_post = NA
    pretty_less120_post = NA
    pretty_less110_post = NA
    pretty_in_EG_post = NA
    pretty_above120_post = NA
  }
  
  knowledge = c("Forecast", "BTF Only", "Updated")
  df = data.frame(c(pretty_H_prior, pretty_H_like, pretty_H_post),
                  c(pretty_ES_prior, pretty_ES_like, pretty_ES_post),
                  c(pretty_less65_prior, pretty_less65_like, pretty_less65_post),
                  c(pretty_in_EG_prior, pretty_in_EG_like, pretty_in_EG_post),
                  c(pretty_above120_prior, pretty_above120_like, pretty_above120_post),
                  c(pretty_less95_prior, pretty_less95_like, pretty_less95_post),
                  c(pretty_less110_prior, pretty_less110_like, pretty_less110_post))
  
  df = t(df)
  df = cbind(c("Additional Harvest", "Expected Escapement", "Pr(S < 65,000)", "Pr(65,000 < S < 120,000)", "Pr(S > 120,000)", "Pr(S < 95,000)", "Pr(S < 110,000)"), df)
  colnames(df) = c("Quantity", knowledge)
  
  if (length(out$PolicyOut$know.ind) == 0 | is.null(out$EstOut)) {
    NULL
  } else {
    df[,c(1,out$PolicyOut$know.ind+1)]
  }
}
##### FUNCTION TO CALCULATION COMPUTATION TIMES #####
calc_mcmc_time = function(){
  prior_fun = approxfun(density(rnorm(1e6, 150000, 25000), from = 0, to = 5e6))
  like_fun = approxfun(density(rnorm(1e6, 250000, 25000), from = 0, to = 5e6))
  starttime = Sys.time()
  samps = MH(init = 150000, ni = 10000, nb = 1000, nt = 1,
             prop.sig = 0.7, prior.fun = prior_fun, like.fun = like_fun)
  as.numeric(Sys.time() - starttime)/10000
}

calc_prior_like_time = function(){
  starttime = Sys.time()
  fit.eos = fit.reg.mod(dat = N.eos.dat)
  
  hist.p = filter(hist.btf, date == "6/12" & rt.type %in% c(1,2,3)) %>% select(p.ccpue) %>% unlist %>% unname
  rt.mu = estBetaParams(mean(hist.p), var(hist.p))
  
  prior.samp = exp(rnorm(1e6, log(150000) - 0.5 * cv2sig(0.27)^2, cv2sig(0.27)))
  
  boot_out = sample.likelihood(n.boot = 1e6,
                               obs.ccpue = 21,
                               fit.eos = fit.eos,
                               rt.method = "beta",
                               rt.mu = rt.mu, 
                               sub.harv = 0, 
                               sub.cv = 0,
                               com.harv = 0
  )
  as.numeric(Sys.time() - starttime)
}


##### FUNCTION TO CALCULATE BIN DIMENSIONS #####
get_bin_dim = function(min = 0, max = 8e5, n) {
  bins = seq(min, max, length = n)
  n = length(bins) - 1
  mp = numeric(n)
  for (i in 1:n) {
    mp[i] = (bins[i+1] + bins[i])/2
  }
  
  list(
    min = min,
    max = max,
    bins = bins,
    n = n,
    mp = mp
  )
}

##### FUNCTION TO CALCULATE BIN PROBS #####
get_bin_probs = function(x, bin_dim, type = "include") {
  # use exlclude for plotting ONLY
  if (type == "include") {
    x[x <= bin_dim$min] = bin_dim$min
    x[x >= bin_dim$max] = bin_dim$max
  }
  
  if (type == "exclude") {
    x = x[x >= bin_dim$min & x <= bin_dim$max]
  }
  
  counts = hist(x, breaks = bin_dim$bins, plot = F)$counts
  probs = counts/sum(counts)
  probs
}

# bin_dims = get_bin_dim(min = 0, max = 8e5, n = 160)
# 
# p = get_bin_probs(x = rlnorm(1e6, log(45000), 0.27) - 100000, bin_dim = bin_dims)
# 
# windows()
# plot(cumsum(p) ~ bin_dims$mp, type = "l", ylim = c(0,1), xlim = c(0,300000))
# lines(cumsum(p) ~ bin_dims$mp)

##### LOGNORMAL SD FUNCTIONS #####
cv2sig = function(x) {
  sqrt(log(x^2 + 1))
}

sig2cv = function(x) {
  sqrt(exp(x^2) - 1)
}

##### FIND.RISK.P() #####
# uses estimated CDF to get P(S < Sobj) if Hcand fish are harvested (direction = "less")
# uses estimated CDF to get P(S > Sobj) if Hcand fish are harvestedv(direction = "greater")
find.risk.p = function(Sobj, Hcand, cdf, direction) {
  
  if (!(direction %in% c("less", "greater"))) {
    stop("direction must be one of 'less' or 'greater'")
  }
  
  if (direction == "less") {
    p = cdf(Sobj + Hcand)
  } else {
    p = 1 - cdf(Sobj + Hcand)
  }
  return(p)
}

##### TURN.TO.JDAY #####
turn.to.jday = function(dates) {
  # first, split up dates
  date.split = matrix(unlist(strsplit(as.character(dates), "/")), nrow = length(dates), ncol = 3, byrow = T)
  date.split = data.frame(day = as.numeric(date.split[,2]), month = as.numeric(date.split[,1]), year = as.numeric(date.split[,3]))
  
  # error handling
  if (length(unique(date.split$year)) != 1) stop("This function can only deal with one year at a time; consider using a loop through years or use dplyr")
  
  # next, decide if the year is a leap year or not
  is.whole.number = function(x) {
    x == round(x)
  }
  
  year = unique(date.split$year)
  yr.type = ifelse(!is.whole.number(year/4), "common", 
                   ifelse(!is.whole.number(year/100), "leap", 
                          ifelse(!is.whole.number(year/400), "common", "leap")))
  yr.type = unique(yr.type)
  
  # make keys
  # leap year key
  dates = seq(as.Date("2016-01-01"), as.Date("2016-12-31"), 1)
  dates = gsub(x = dates, pattern = "2016-", replacement = "")
  dates = gsub(x = dates, pattern = "-", replacement = "/")
  month.key = as.integer(substr(dates, 1, 2))
  day.key = as.integer(substr(dates, 4, 5))
  leap.key = data.frame(day = day.key, month = month.key, jday = 1:366)
  
  # common year key
  dates = seq(as.Date("2015-01-01"), as.Date("2015-12-31"), 1)
  dates = gsub(x = dates, pattern = "2015-", replacement = "")
  dates = gsub(x = dates, pattern = "-", replacement = "/")
  month.key = as.integer(substr(dates, 1, 2))
  day.key = as.integer(substr(dates, 4, 5))
  common.key = data.frame(day = day.key, month = month.key, jday = 1:365)
  
  if (unique(yr.type) == "leap") key = leap.key else key = common.key
  
  n = nrow(date.split)
  jdays = numeric(n)
  for (i in 1:n) {
    jdays[i] = key$jday[key$month == date.split$month[i] & key$day == date.split$day[i]]
  }
  
  return(jdays)
}

##### TURN.FROM.JDAY #####
turn.from.jday = function(jday, year) {
  # error handling
  if(missing(year)) stop("Must specify the year, otherwise you can't account for leap years!")
  if(missing(jday)) stop("Must specify the jday, otherwise you can't get a date!")
  if(length(year) > 1 | length(jday) > 1) stop("jday and year must be vectors of length 1")
  
  # decide if the year is a leap year
  is.whole.number = function(x) {
    x == round(x)
  }
  yr.type = ifelse(!is.whole.number(year/4), "common", 
                   ifelse(!is.whole.number(year/100), "leap", 
                          ifelse(!is.whole.number(year/400), "common", "leap")))
  yr.type = unique(yr.type)
  
  # make keys
  # leap year key
  dates = seq(as.Date("2016-01-01"), as.Date("2016-12-31"), 1)
  dates = gsub(x = dates, pattern = "2016-", replacement = "")
  dates = gsub(x = dates, pattern = "-", replacement = "/")
  month.key = as.integer(substr(dates, 1, 2))
  day.key = as.integer(substr(dates, 4, 5))
  leap.key = data.frame(day = day.key, month = month.key, jday = 1:366)
  
  # common year key
  dates = seq(as.Date("2015-01-01"), as.Date("2015-12-31"), 1)
  dates = gsub(x = dates, pattern = "2015-", replacement = "")
  dates = gsub(x = dates, pattern = "-", replacement = "/")
  month.key = as.integer(substr(dates, 1, 2))
  day.key = as.integer(substr(dates, 4, 5))
  common.key = data.frame(day = day.key, month = month.key, jday = 1:365)
  
  if (jday > 365 & yr.type == "common") stop("There is no 366th day in a common year!")
  if (yr.type == "leap") key = leap.key else key = common.key
  
  date.info = key[key$jday == round(jday),]
  
  day = as.character(date.info$day); if(nchar(day) < 2) day = paste("0", day, sep = "")
  month = as.character(date.info$month); if(nchar(month) < 2) month = paste("0", month, sep = "")
  
  date = paste(month, day, year, sep = "/")
  return(date)
  
}

##### FIX.DATES() #####
fix.dates = function(x) {
  patterns = paste("0", 1:9, sep = "")
  replace = as.character(1:9)
  
  new.x1 = character(length(x))
  new.x2 = character(length(x))
  for (i in 1:length(x)) {
    for (j in 1:length(patterns)) {
      if (str_detect(x[i], patterns[j])) new.x1[i] = str_replace_all(x[i], patterns[j], replace[j])
      if (new.x1[i] == "") new.x1[i] = x[i]
    }
    
    for (j in 1:length(patterns)) {
      if (str_detect(new.x1[i], patterns[j])) new.x2[i] = str_replace_all(new.x1[i], patterns[j], replace[j])
      if (new.x2[i] == "") new.x2[i] = new.x1[i]
    }
  }
  
  return(new.x2)
}

##### HIST.INDEX.PREP() #####
hist.index.prep = function(f.date = "6/1", l.date = "8/24", chinook, chum, sockeye, coho) {
  
  if (!missing(chinook)) {
    keep = seq(which(chinook$Date == f.date), which(chinook$Date == l.date),1)
    chinook = chinook[keep,]
    chinook[is.na(chinook)] = 0
    chinook$species = "Chinook"
  } else {
    chinook = NULL
  }
  
  if (!missing(chum)) {
    keep = seq(which(chum$Date == f.date), which(chum$Date == l.date),1)
    chum = chum[keep,]
    chum[is.na(chum)] = 0
    chum$species = "Chum"
  } else {
    chum = NULL
  }
  
  if (!missing(sockeye)) {
    keep = seq(which(sockeye$Date == f.date), which(sockeye$Date == l.date),1)
    sockeye = sockeye[keep,]
    sockeye[is.na(sockeye)] = 0
    sockeye$species = "Sockeye"
  } else {
    sockeye = NULL
  }
  
  if (!missing(coho)) {
    keep = seq(which(coho$Date == f.date), which(coho$Date == l.date),1)
    coho = coho[keep,]
    coho[is.na(coho)] = 0
    coho$species = "Coho"
  } else {
    coho = NULL
  }
  
  index = rbind(chinook, chum, sockeye, coho)
  
  # manipulate dates
  colnames(index) = c("date", substr(colnames(index)[2:(ncol(index)-1)], 2, 5), "species")
  years = as.numeric(colnames(index)[2:(ncol(index)-1)])
  nyrs = length(years)
  
  index = cbind(day = seq(dates$day[dates$date == f.date], dates$day[dates$date == l.date]), index)
  index = melt(index, id.vars = c("day", "date", "species"), value.name = "cpue", variable.name = "year")
  index$year = as.numeric(as.character(index$year))
  index$date = as.character(index$date)
  
  # index = index %>% filter(year < c.year)
  
  index = index %>% group_by(year) %>% mutate(jday = turn.to.jday(paste(date, year, sep = "/")))
  index = index %>% group_by(year, species) %>%
    mutate(ccpue = cumsum(cpue), p.ccpue = ccpue/sum(cpue)) %>%
    ungroup %>%
    select(day, jday, date, year, species, cpue, ccpue, p.ccpue) %>%
    arrange(year, species, day)
  
  return(index)
}

##### LOGISTIC() #####
logistic = function(x, d50, h) {
  1/(1 + exp(-h * (x - d50)))
}

##### FIT.TIMING() #####
fit.timing = function(hist, spp, init = NULL) {
  
  if(is.null(init)) {
    init = numeric(2)
    init[1] = ifelse(spp == "Chinook", 167, 
                     ifelse(spp == "Chum", 172, 
                            ifelse(spp == "Sockeye", 172, 
                                   ifelse(spp == "Coho", 212, NA))))
    init[2] = ifelse(spp == "Chinook", 0.1, 
                     ifelse(spp == "Chum", 0.1, 
                            ifelse(spp == "Sockeye", 0.3, 
                                   ifelse(spp == "Coho", 0.3, NA))))
  }
  
  names(init) = c("d50", "h")
  
  # model to fit run timing curves
  mod = function(theta, df) {
    d50 = theta["d50"]
    h = theta["h"]
    
    pred.p.ccpue = logistic(x = df$jday, d50 = d50, h = h)
    
    resid = df$p.ccpue - pred.p.ccpue
    sqr.resid = resid^2
    ssq = sum(sqr.resid)
    out = list(resid = resid, sqr.resid = sqr.resid, ssq = ssq)
    return(out)
  }
  
  # gives the sum of sqrs
  mod.ssq = function(theta, df) {
    mod(theta, df)$ssq
  }
  
  years = unique(hist$year)
  nyrs = length(years)
  
  ests = data.frame(year = years, d50 = NA, h = NA, converge = NA)
  y = 1
  for (y in 1:nyrs) {
    temp.fit = optim(par = init, fn = mod.ssq, df = filter(hist, species == spp & year == years[y]), control = list(maxit = 1000))
    ests[y,"d50"] = temp.fit$par["d50"]
    ests[y,"h"] = temp.fit$par["h"]
    ests[y,"converge"] = temp.fit$convergence
  }
  
  if (any(ests$converge == 1)) {
    n.bad = sum(ests$converge == 1)
    warning (paste(n.bad, "years did not converge. Returning all years that did"))
    return(ests[ests$converge != 1,c("year", "d50", "h")])
  } else {
    return(ests[,c("year", "d50", "h")])
  }
}

##### DATES.PREP() #####
dates.prep = function(f.date = "6/1", l.date = "8/24", year) {
  
  # julian days
  f.jday = turn.to.jday(paste(f.date, year, sep = "/"))
  l.jday = turn.to.jday(paste(l.date, year, sep = "/"))
  jday.seq = f.jday:l.jday
  ndays = length(jday.seq)
  
  # days past may 31st
  days = jday.seq - (min(jday.seq) - 1)
  
  # create character dates
  dates = character(ndays)
  for (i in 1:ndays) {
    dates[i] = turn.from.jday(jday.seq[i], year = year)
  }
  dates = fix.dates(substr(dates, 1, nchar(dates) - 5))
  
  # create data frame with all days, jdays, and dates
  dates = data.frame(day = days, jday = jday.seq, date = dates, stringsAsFactors = F)
  
  return(dates)
}

############################################
#### FUNCTIONS FOR LIKELIHOOD MACHINERY ####
############################################

##### FUNCTION TO GENERATE RANDOM HARVEST SAMPLES #####
gen_harv_samps = function(n.boot, mean, cv) {
  
  if (is.na(mean) | mean < 0) 
    mean = 0
  if (is.na(cv) | cv < 0) 
    cv = 0
  
  if (mean == 0) {
    x = rep(0, n.boot)
  } else {
    if (cv == 0) {
      x = rep(mean, n.boot)
    } else {
      x = exp(rnorm(n.boot, log(mean) - 0.5 * cv2sig(cv)^2, cv2sig(cv)))
    }
  }
  
  return(x)
}

##### FUNCTION TO GET BETA PARAMETERS #####
estBetaParams <- function(mu, var) {
  alpha = mu * ((mu * (1 - mu))/var - 1)
  beta = (1 - mu) * ((mu * (1 - mu))/var - 1)
  
  return(c(alpha, beta))
}

##### FUNCITON TO FIT THE N VS. EOS REGRESSION MODEL #####
fit.reg.mod = function(dat) {
  
  fit = lm(log(N.btf) ~ log(eos), data = dat)
  
  list(
    coef = coef(fit),
    sig = summary(fit)$sigma,
    vcov = vcov(fit)
  )
}

##### FUNCTIONS FOR CREATING FILE NAMES #####

date_for_filename = function(dates) {
  unname(sapply(dates, function(z) paste(unlist(str_split(z, pattern = "/")), collapse = "_")))
}

create_export_filename = function(date, suffix, ext, base = "Ests_") {
  paste(base, paste(date_for_filename(date), collapse = "_thru_"),
        ifelse(nchar(suffix) == 0, "", paste("_", suffix, sep = "")), ext, sep = "")
}

##### FUNCTION TO PREDICT FROM THE REGRESSION MODEL #####
pred.reg.mod = function(fit, rand.eos, n.boot) {
  rand.coefs = rmvnorm(n.boot, fit$coef, fit$vcov)
  rand.e = rnorm(n.boot, -0.5 * fit$sig^2, fit$sig)
  
  rand.log.N = rand.coefs[,1] + rand.coefs[,2] * log(rand.eos) + rand.e
  rand.N = unname(exp(rand.log.N))
  rand.N
}

##### FUNCTION TO SAMPLE FROM LIKELIHOOD DISTRIBUTION #####
sample.likelihood = function(n.boot = 10000, obs.ccpue, fit.eos, rt.mu, 
                             sub.harv = 0, sub.cv = 0, com.harv = 0, rt.method = "beta") {
  
  if (rt.method == "beta") {
    rand.p = rbeta(n.boot, rt.mu[1], rt.mu[2])
  }
  
  if (rt.method == "known") {
    rand.p = rep(rt.mu, n.boot)
  }
  
  # STEP 3: OBTAIN THE END OF SEASON CCPUE SUGGESTED BY EACH COMPLETED PROPORTION
  rand.eos = obs.ccpue/rand.p
  rand.log.eos = log(rand.eos)
  
  # STEP 4: GENERATE RANDOM RUN SIZES PAST THE BTF
  rand.N.btf = pred.reg.mod(fit = fit.eos, rand.eos = rand.eos, n.boot = n.boot)
  
  # STEP 5: ADD ON HARVEST BELOW BETHEL TO DATE
  rand.sub.harv = gen_harv_samps(n.boot, mean = sub.harv, cv = sub.cv)
  
  rand.N = rand.N.btf + rand.sub.harv + com.harv
  rand.log.N = log(rand.N)
  
  # STEP 6: SUMMARIZE
  rand.out = cbind(p = rand.p, eos = rand.eos, sub.harv = rand.sub.harv, N.btf = rand.N.btf, N = rand.N)
  # rand.out = rand.out[rand.out[,"N"] <= maxN,]
  
  # STEP 7: RETURN OUTPUT
  return(rand.out)
}

############################################
##### FUNCTIONS FOR BAYESIAN MACHINERY #####
############################################

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

##### FUNCTION TO PLOT DIAGNOSTICS: TRACE PLOT AND DENSITY PLOTS #####
diag.plots = function(post.samp, n.trace = 20000) {
  
  n.chain = ncol(post.samp)
  n.samp = nrow(post.samp)
  
  cols = rep(c("blue", "red", "green", "lightskyblue1", "yellow", "pink", "purple"), 10)
  
  x.lims = matrix(NA, n.chain, 2)
  y.max = numeric(n.chain)
  post.samp = post.samp/1000
  
  for (c in 1:n.chain) {
    temp.dens = density(post.samp[,c])
    
    x.lims[c,] = range(temp.dens$x)
    y.max[c] = max(temp.dens$y)
  }
  
  xlim = c(min(x.lims[,1]), max(x.lims[,2]))
  ylim = c(0, max(y.max)) * 1.05
  
  par(mfrow = c(1,2), mar = c(4,4,2,0.5), xaxs = "i", yaxs = "i", cex.main = 1.4, cex.lab = 1.4, cex.axis = 1.3)
  plot(1, 1, type = "n", xlim = xlim, ylim = ylim, main = "Density Plot", xaxs = "i", las = 1,
       xlab = "", ylab = "", yaxt = "n", xaxt = "n")
  mtext(side = 1, "Value of N (1000s)", line = 3, cex = 1.4)
  mtext(side = 2, "Density", line = 2, cex = 1.4)
  for (c in 1:n.chain) {
    lines(density(post.samp[,c]), col = cols[c], lwd = 4)
  }
  axis(side = 1, lwd = 4)
  box(lwd = 4)
  
  # thin the traceplot so doesn't take forever to render
  if (n.samp < n.trace) keep.i = 1:n.samp else keep.i = sort(sample(x = 1:n.samp, size = n.trace, replace = F))
  post.samp.all = post.samp
  post.samp = post.samp[keep.i,]
  plot(1, 1, type = "n", ann = F, axes = F, xlim = c(0,n.samp))
  at.x = axisTicks(par("usr")[1:2], log = F)
  lab.x = at.x/1000
  
  ylim = range(post.samp.all) * c(0.95, 1.05)
  par(mar = c(4,0.5,2,4), xaxs = "r", new = T)
  plot(1,1, type = "n", xlim = c(0, min(c(n.trace, n.samp))), ylim = ylim, xaxt = "n",
       yaxt = "n", xlab = "Iteration (1000s)", main = "Trace Plot")
  mtext(side = 4, "Value of N (1000s)", line = 2, cex = 1.4)
  for (c in 1:n.chain) {
    lines(post.samp[,c], col = cols[c])
  }
  axis(side = 4, lwd = 4)
  if (n.samp > n.trace) {
    axis(side = 1, at = (at.x/n.samp) * n.trace, lab = lab.x, lwd = 4)
  } else {
    axis(side = 1, at = at.x, lab = lab.x, lwd = 4)
  }
  
  box(lwd = 4)
  
}
# 
# windows()
# ni = 100000; nc = 6
# post.samp = matrix(rnorm(ni * nc, 500, 10), ni, nc)
# diag.plots(post.samp)
##### FUNCTION TO CALCULATE POSTERIOR SUMMARY STATS #####
summ = function(x) {
  c(mean = mean(x), sd = sd(x), quantile(x, c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)))
}

##### FUNCTION TO CREATE PROGRESS BAR #####
progress_bar = function(i, ni) {
  print.i = round(seq(ni/100, ni, length = 100))
  if (i %in% c(1, print.i)) {
    j = which(i == c(1, print.i))
    j = j - 1
    
    if (j == 0) {
      cat("|-", rep(" ", 100), " -|", paste(" ", 0, "%", sep = ""), sep = "")
    } else {
      cat("\r", "|- ", rep("*", j), rep(" ", 100-j), " -|", paste(" ", j, "%", sep = ""), sep = "")
    }
  }
}

##### FUNCTIONS TO WRAP PDF FORMATION, MCMC, AND POST PROCESSING INTO ONE #####
run.mcmc = function(init, ni = 100000, nb = 1000, nt = 2, nc = 2, prop.sig, prior.fun, like.fun) {
  
  # containers
  post.samp = matrix(NA, (ni - nb)/nt, nc)
  accept.rate = numeric(nc)
  
  # call MH separately for each chain
  for (c in 1:nc) {
    # cat("\n Chain ", c, " Sampling:\n ", sep = "")
    incProgress(1/nc, detail = paste("Chain #", c, sep = ""))
    MH.temp = MH(init = init[c], 
                 ni = ni, nb = nb, nt = nt,
                 prop.sig = prop.sig, 
                 like.fun = like.fun, 
                 prior.fun = prior.fun)
    
    post.samp[,c] = MH.temp$post.samp
    accept.rate[c] = MH.temp$accept.rate
  }
  
  output = list(
    post.samp = post.samp,
    accept.rate = accept.rate
  )
  
  return(output)
}

##### FUNCTION TO FIND INDEX OF CLOSEST MATCHING ELEMENT #####
closest.index = function(x, p) {
  which.min(abs(x - p))
}

