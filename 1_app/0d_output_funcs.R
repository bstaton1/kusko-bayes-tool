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
    # par(cex.lab = 1.5, cex.axis = 1.4, xaxs = "i", yaxs = "i", mar = c(5,3,1,1))
    # plot(1,1, type = "n",
    #      col = "black", xlim = c(0, 400000), ylim = c(0, 1), lwd = 2,
    #      xaxt = "n", yaxt = "n", ylab = "",
    #      xlab = "Run Size (1000s)")
    # mtext(side = 2, line = 1, "Relative Probability", cex = 1.5)
    # axis(side = 1, at = seq(0, 400000, 50000), labels = seq(0, 400, 50), lwd = 4)
    # box(lwd = 4)
    NULL
  }
}

##### FUNCTION TO CREATE ESTIMATE TABLE #####
create_est_table = function(out) {
  pretty_prior = prettyNum(round(summ(out$EstOut$prior_samp), -3), big.mark = ",", scientific = F)
  pretty_like = prettyNum(round(summ(out$EstOut$boot_out[,"N"]), -3), big.mark = ",", scientific = F)
  rnames = names(summ(rnorm(10)))
  rnames[1:2] = c("Mean", "SD")
  if (!is.null(out$EstOut$post_samp)) {
    cnames = c("Statistic", "Forecast", "BTF Only", "Updated")
    pretty_post = prettyNum(round(summ(as.numeric(out$EstOut$post_samp)), -3), big.mark = ",", scientific = F)
    tab = data.frame(rnames, pretty_prior, pretty_like, pretty_post)
    colnames(tab) = cnames
  } else {
    cnames = c("Statistic", "Forecast", "BTF Only")
    tab = data.frame(rnames, pretty_prior, pretty_like)
    colnames(tab) = cnames
  }
  # tab[,!(colnames(tab) %in% c("mean", "sd"))]
  tab
}
##### FUNCTION TO CREATE ESTPLOT2 (RELATIONSHIP) #####
create_EstPlot2 = function(input, out) {
  if (!is.null(out$EstOut$boot_out)) {
    layout(matrix(c(1,2), 1, 2), widths = c(1,0.3))
    par(mar = c(2,3,1,0), xaxs = "i", yaxs = "i", oma = c(2,2,0,0))
    
    fit_data = prepare_fit_data(dt = out$EstOut$cdate, index = btf_data, N = N_data)
    fit = lm(log(N) ~ q * ccpue, data = fit_data)
    line_x = with(filter(fit_data, q == 2), seq(floor(min(ccpue)), ceiling(max(ccpue)), length = 100))
    line_y = sample_likelihood(fit, pred_ccpue = line_x, pred_q = 2, n_mc = 0)[,"N"]
    point_y = sample_likelihood(fit, pred_ccpue = out$EstOut$ccpue, pred_q = 2, n_mc = 0)[,"N"]
    rand_y = out$EstOut$boot_out[,"N"]
    
    # scatterplot relationship
    plot(N ~ ccpue, data = filter(fit_data, q == 2),
         type = "n", xaxt = "n", yaxt = "n",
         ylim = c(input$EstPlot2_ylim[1], input$EstPlot2_ylim[2]),
         # xlim = c(0, max(filter(btf_data, date == "8/24" & year >= 2008)$ccpue) * 1.05))
         xlim = c(input$EstPlot2_xlim[1], input$EstPlot2_xlim[2]))
    usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
    lines(line_y ~ line_x, lwd = 2)
    segments(out$EstOut$ccpue, usr[3], out$EstOut$ccpue, point_y, lwd = 2, lty = 2, col = "blue", lend = 2)
    segments(out$EstOut$ccpue, point_y, usr[2], point_y, lwd = 2, lty = 2, col = "blue", lend = 2)
    with(filter(fit_data, q == 2),
         text(x = ccpue, y = N, labels = paste("'", substr(year, 3, 4), sep = ""), cex = 1.2, font = 2))
    axis(side = 1, at = seq(0, 1200, 100), labels = seq(0, 1200, 100), lwd = 4, cex.axis = 1.2)
    axis(side = 2, at = seq(0, 450000, 50000), labels = seq(0, 450, 50), las = 2, lwd = 4, cex.axis = 1.2)
    
    if (input$EstPlot2_predata) {
      with(filter(fit_data, q == 1),
           text(x = ccpue, y = N, labels = paste("'", substr(year, 3, 4), sep = ""), cex = 1.2, font = 2, col = "grey"))
      
      line_x = with(filter(fit_data, q == 1), seq(floor(min(ccpue)), ceiling(max(ccpue)), length = 100))
      line_y = sample_likelihood(fit, pred_ccpue = line_x, pred_q = 1, n_mc = 0)[,"N"]
      lines(line_y ~ line_x, lwd = 2, col = "grey")
    }
    
    
    box(lwd = 4, ljoin = 1)
    
    par(mar = c(2,0,1,0.5))
    dens = density(rand_y, from = 0, to = 400000)
    plot(dens$x ~ dens$y, 
         # ylim = c(0, 400000), 
         ylim = c(input$EstPlot2_ylim[1], input$EstPlot2_ylim[2]),
         type = "l", axes = F,
         xlim = c(0, max(dens$y) * 1.1), xlab = "", ylab = "")
    polygon(x = c(dens$y, rev(dens$y)),
            y = c(dens$x, rep(0, length(dens$x))),
            col = "skyblue2", border = NA)
    lines(dens$x ~ dens$y, lwd = 2)
    usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
    segments(usr[1], mean(point_y),
             approxfun(density(rand_y, from = 0, to = 400000))(point_y),
             point_y, lwd = 2, lty = 2, col = "blue")
    box(lwd = 4)
    
    mtext(side = 1, outer = T, "BTF CCPUE", cex = 1.2, line = 0.5)
    mtext(side = 2, outer = T, "Run Size (1000s)", cex = 1.2, line = 0.5)
  } else {
    NULL
  }
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

##### FUNCTION TO CREATE POLICY TABLE 1 #####
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



