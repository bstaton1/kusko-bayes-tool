##### SESSION SETUP #####
library(stringr)
library(reshape2)
library(gifski)
library(dplyr)

# CLEAR THE WORKSPACE
rm(list = ls(all = TRUE))

# LOAD ALL NECESSARY FUNCTIONS 
# SET THE WORKING DIRECTORY TO THE LOCATION OF THIS FILE BEFORE RUNNING REST OF CODE
source("../../1_app/0b_util_funcs.R")
source("../../1_app/0c_update_funcs.R")
source("../../1_app/0d_output_funcs.R")

data_dir = "../../1_app/Inputs"

# PREPARE DATA FILES
# a date lookup key
dates = dates.prep(year = 2022)

# historical BTF data
btf_data = hist.index.prep(
  f.date = "6/1", l.date = "8/24",
  chinook = read.csv(file.path(data_dir, "Chinook Daily BTF.csv"), stringsAsFactors = F)
)

# historical total run abundance data
N_data = read.csv(file.path(data_dir, "run_size_ests.csv"), stringsAsFactors = F)

# graphics device resolution (pixels per inch)
ppi = 600

# function to fit regression and plot data for one day of the season
make_plot = function(cdate) {
  
  # format date for file name
  date_split = unlist(stringr::str_split(cdate, "/"))
  date_split = stringr::str_pad(date_split, 2, "left", "0")
  file = paste0(paste(date_split, collapse = "-"), ".png")
  
  # open a graphics device file
  png(file.path(fig_dir, file), height = 4 * ppi, width = 5 * ppi, res = ppi)

  # extract the data for regression on this day
  fit_data = prepare_fit_data(dt = cdate, index = btf_data, N = N_data)
  
  # fit this day's regression
  fit = lm(log(N) ~ q * ccpue, data = fit_data)
  
  # obtain values along fitted curve
  pred_x = with(subset(fit_data, q == 2), seq(min(floor(ccpue)), max(ceiling(ccpue)), length = 50))
  pred_y = exp(predict(fit, newdata = data.frame(ccpue = pred_x, q = factor(2))))
  
  # obtain predicted values of interest
  pred_x10 = with(subset(fit_data, q == 2), quantile(ccpue, 0.1))
  pred_y10 = exp(predict(fit, newdata = data.frame(ccpue = pred_x10, q = factor(2))))
  pred_x90 = with(subset(fit_data, q == 2), quantile(ccpue, 0.9))
  pred_y90 = exp(predict(fit, newdata = data.frame(ccpue = pred_x90, q = factor(2))))
  
  # graphical parameters
  par(mar = c(3,3,1,1), mgp = c(2,0.35,0), tcl = -0.15, cex.lab = 1.2, font.lab = 2,
      xaxs = "i", yaxs = "i", lend = "square", ljoin = "mitre")
  
  # create a blank plot with proper dimensions and labeling
  plot(N ~ ccpue, data = subset(fit_data, q == 2), type = "n",
       xlim = c(0, 1000),
       ylim = c(50000, 250000),
       xlab = paste0("BTF CCPUE as of ", cdate),
       ylab = "Total Run Abundance (1000s)", yaxt = "n")
  
  # draw better axes
  axis(side = 2, at = seq(0, 400000, 25000), labels = seq(0, 400, 25), las = 2, lwd = 2)
  axis(side = 1, at = seq(0, 1000, 100), lwd = 2)
  
  # add prediction points
  segments(pred_x10, 0, pred_x10, pred_y10, col = "blue")
  segments(0, pred_y10, pred_x10, pred_y10, col = "blue", lty = 2)
  segments(pred_x90, 0, pred_x90, pred_y90, col = "red")
  segments(0, pred_y90, pred_x90, pred_y90, col = "red", lty = 2)
  
  # add observed data
  text(N ~ ccpue, data = subset(fit_data, q == 2), 
       col = scales::alpha("grey20", 0.5),
       labels = substr(year, 3, 4), font = 2, xpd = TRUE)
  
  # add fitted line
  lines(pred_y ~ pred_x, lwd = 2)
  
  # add legend
  legend("bottomright", title = "BTF CCPUE\nQuantile", bty = "n",
         cex = 0.8, legend = c("10%", "90%"), lty = 1, col = c("blue", "red"))
  
  # draw boundary box
  box(lwd = 2)
  
  # close device
  junk = dev.off()
}

# a test plot
# make_plot("6/30")

# which dates to plot
fdate = "6/5"      # first one
ldate = "7/30"     # last one
date_interval = 1  # incremented by this many days

# extract the appropriate dates
fday = dates$day[dates$date == fdate]
lday = dates$day[dates$date == ldate]
day_range = seq(fday, lday, by = date_interval)
if (day_range[length(day_range)] != lday) day_range = c(day_range, lday)
date_range = dates$date[dates$day %in% day_range]

# loop through dates creating the plot file for each
fig_dir = "daily-figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir)
junk = sapply(date_range, function(d) {
  cat("\rMaking Plot: ", d, " (", which(date_range == d), " of ", length(date_range), ")", sep = "")
  make_plot(d)
})

# frames per second
fps = 4

# get all png file names
files = list.files(fig_dir, pattern = "*\\.png$", full.names = TRUE)

# number of seconds to view the gif
length(files)/fps

# create the gif
gifski(png_files = files, gif_file = "regression-progression.gif", width = 5 * ppi, height = 4 * ppi, delay = 1/fps, loop = TRUE, progress = TRUE)

# delete individual files
unlink(fig_dir, recursive = TRUE)
