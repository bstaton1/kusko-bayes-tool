# query the list of all installed packages on this machine
pkg_have = rownames(installed.packages())

# the list of all packages used by the BayesTool
pkg_need = c(
  "shiny", "shinythemes", "shinyBS",
  "shinyjs", "shinyWidgets", "miniUI",
  "dplyr", "stringr", "reshape2",
  "mvtnorm", "scales", "coda",
  "gridExtra"
  )

# install the missing packages
if (all(pkg_need %in% pkg_have)) {
  cat("All necessary packages already installed")
} else {
  for (i in 1:length(pkg_need)) {
    if (!(pkg_need[i] %in% pkg_have)) {
      cat("Installing required package:", pkg_need[i], "\n")
      install.packages(pkg_need[i], quiet = T)
    }
  }
  
  cat("All necessary packages installed")
}
