# directory of historical data
data_dir = paste(getwd(), "Inputs", sep = "/")
# data_dir = "../Inputs"
# data_dir = "Inputs"
# prepare btf data for historical years
dates = dates.prep(year = 2017)
hist.btf = hist.index.prep(f.date = "6/1", l.date = "8/24", chinook = read.csv(paste(data_dir, "Chinook Daily BTF.csv", sep = "/"), stringsAsFactors = F))
eos = hist.btf %>% filter(date == "8/24" & year %in% 2008:2017) %>% select(ccpue) %>% unlist %>% unname

# prepare run size information: harvest information is only downstream of btf
# N.eos.dat = read.csv(paste(data_dir, "run_size_ests.csv", sep = "/"))
N.eos.dat = read.csv(paste(data_dir, "run_size_ests_new.csv", sep = "/"))
N.eos.dat$eos = eos

# add a run timing category to the data set based on the midpoint date
rt.ests = fit.timing(hist.btf, spp = "Chinook")
q = quantile(rt.ests$d50, c(0.33, 0.66))
# q = seq(min(rt.ests$d50), max(rt.ests$d50), length = 4)[2:3]
rt.types = data.frame(year = rt.ests$year, rt.type = ifelse(rt.ests$d50 <= q[1], 1, ifelse(rt.ests$d50 <= q[2], 2, 3)))
hist.btf = merge(hist.btf, rt.types, by = "year")
hist.btf = arrange(hist.btf, year, day)
