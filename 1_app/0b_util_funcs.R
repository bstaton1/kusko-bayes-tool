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
  
  # a date lookup key
  dates = dates.prep(year = 2019)
  
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

##### DATES.PREP()
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

##### FUNCTION TO PRETTIFY NUMBERS #####
prettify = function(x, rnd = NULL) {
  if (!is.null(rnd)) x = round(x, rnd)
  out = prettyNum(x = x, big.mark = ",", scientific = F)
  
  return(out)
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
##### LOGNORMAL SD FUNCTIONS #####
cv2sig = function(x) {
  sqrt(log(x^2 + 1))
}

sig2cv = function(x) {
  sqrt(exp(x^2) - 1)
}

##### FUNCTION TO CALCULATE COMPUTATION TIMES #####
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
  # subset historical data
  fit_data = prepare_fit_data(dt = "6/12", index = btf_data, N = N_data)
  
  # fit the historical regression
  fit = lm(log(N) ~ q * ccpue, data = fit_data)
  prior.samp = exp(rnorm(1e6, log(150000) - 0.5 * cv2sig(0.28)^2, cv2sig(0.28)))
  
  # generate N samples based on index
  boot_out = sample_likelihood(fit, pred_ccpue = 100, pred_q = 2, n_mc = 1e6)
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


##### FUNCTIONS FOR CREATING FILE NAMES #####

date_for_filename = function(dates) {
  unname(sapply(dates, function(z) paste(unlist(str_split(z, pattern = "/")), collapse = "_")))
}

create_export_filename = function(date, suffix, ext, base = "Ests_") {
  paste(base, paste(date_for_filename(date), collapse = "_thru_"),
        ifelse(nchar(suffix) == 0, "", paste("_", suffix, sep = "")), ext, sep = "")
}
##### FUNCTION TO CALCULATE SUMMARY STATS #####
summ = function(x) {
  c(mean = mean(x), sd = sd(x), quantile(x, c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)))
}

##### FUNCTION TO FIND INDEX OF CLOSEST MATCHING ELEMENT #####
closest.index = function(x, p) {
  which.min(abs(x - p))
}