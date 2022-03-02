
# set working directory to location of THIS file

# load necessary functions (just cv2sig() and sig2cv())
source("../../1_app/0b_util_funcs.R")

# location of data files
data_dir = "../../1_app/Inputs"

# historical total run abundance data
N_data = read.csv(file.path(data_dir, "run_size_ests.csv"), stringsAsFactors = F)

# calculate forecast mean and cv defaults
fcst_mean_def = round(N_data$N[nrow(N_data)], -3)
true_runs = N_data[2:nrow(N_data),"N"]; fcst_runs = N_data[1:(nrow(N_data)-1),"N"]
fcst_cv_def = round(sig2cv(sd(log(fcst_runs) - log(true_runs))),2)

# PRODUCE A FINE-SCALE VERSION OF THE FORECAST QUANTILES TABLE
p_seq = c(0.01, seq(0.05, 0.95, by = 0.05), 0.99)
q_seq = round(qlnorm(p_seq, log(fcst_mean_def) - 0.5 * cv2sig(fcst_cv_def)^2, cv2sig(fcst_cv_def)), -3)
table_fine = cbind(N = q_seq, p_below = p_seq, p_above = 1 - p_seq)
write.csv(table_fine, "quantiles.csv", row.names = FALSE)

# PRODUCE A FIGURE SHOWING THE CUMULATIVE PROBABILITY DENSITY FUNCTION
p_seq = seq(0.01, 0.99, by = 0.01)
q_seq = qlnorm(p_seq, log(fcst_mean_def) - 0.5 * cv2sig(fcst_cv_def)^2, cv2sig(fcst_cv_def))

# location of horizontal and vertical ticks and grid lines
h_grid = seq(0, 1, by = 0.1)
v_grid = seq(10000, 450000, by = 10000)

# resolution
ppi = 600

# open a graphics device
png("CDF.png", width = 6 * ppi, height = 4 * ppi, res = ppi)

# set graphics parameters
par(mar = c(2.5,3,0.5,0.5), xaxs = "i", yaxs = "i", mgp = c(2,0.25,0), tcl = -0.15, cex.axis = 0.8, lend = "square")

# create a blank plot with correct dimensions
plot(p_seq ~ q_seq, type = "n", ylim = c(0,1), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "")

# draw the polygon of the curve
polygon(x = c(q_seq, rev(q_seq)), y = c(p_seq, rep(0, length(q_seq))), col = "grey50", border = NA)

# add gridlines
abline(h = h_grid, col = "grey80", lty = 3)
abline(v = v_grid, col = "grey80", lty = 3)

# draw the curve itself
lines(p_seq ~ q_seq, lwd = 2)

# draw axis labels
mtext(side = 1, line = 1.5, "Run Size (1000s)")
mtext(side = 2, line = 2, "Probability of Run Smaller than Value")
box(col = "white")

# draw axes
axis(side = 1, at = v_grid, labels = v_grid/1000, las = 2)
axis(side = 2, at = h_grid, labels = paste0(h_grid * 100, "%"), las = 2)

# close the graphics device
dev.off()
