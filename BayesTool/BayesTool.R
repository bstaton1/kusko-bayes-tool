## ------------------------------------------------------------- ##
## IN-SEASON CHINOOK RUN SIZE ASSESSMENT AND DECISION TOOL       ##
## ------------------------------------------------------------- ##

## -------------------------- ##
## CODE BY: BEN STATON        ##
## -------------------------- ##
## DATE V1 STARTED: 3/14/2018 ##
## -------------------------- ##
## DATE UPDATED:     8/2/2018 ##
## -------------------------- ##

## NECESSARY PACKAGES
library(shiny)           # for web apps
library(shinythemes)     # to change the theme
library(shinyBS)         # for pop up text
library(shinyjs)         # more shiny UI features
library(shinyWidgets)    # more shiny UI features
library(miniUI)          # for more UI options
library(dplyr)           # easy data manipulation
library(stringr)         # easy character manipulation
library(reshape2)        # easy data manipulation
library(mvtnorm)         # to draw random regression coefficients
library(scales)          # for transparent colors
library(coda)            # for mcmc diagnostics
library(gridExtra)       # for turning tables to plots

##### SESSION SETUP #####
rm(list = ls())

# Sys.setenv("R_ZIPCMD" = "C:/Rtools/bin/zip.exe")
# Sys.getenv("R_ZIPCMD")

# mcmc defaults
ni_def = 76000     # total samples (including thin and burnin)
nb_def = 1000      # number of burnin
nc_def = 4         # number of chains
prop_sig_def = 0.7 # sd of the lognormal proposal dist.

# forecast mean and cv defaults
fcst_mean_def = 133000
fcst_cv_def = 0.27

# btf ccpue and date defaults
cdate_def = "6/1"
ccpue_def = 0

# downstream harvest mean and cv defaults
ds_harv_mean_def = 0
ds_harv_cv_def = 0.15

# total harvest mean and cv defaults
harv_mean_def = 0
harv_cv_def = 0.15

# risk analysis defaults
H_cand_lim_def = c(0, 50000)
Sobj_def = 92500
p.star_def = 0.5

# number to sample from likelihood and prior
n.boot = 1e6

# maximum run size to calculate density for
dens.max = 5e6

# load in necessary functions to sample from likelihood and conduct MCMC
source("Updating_Functions.R")

# prepare the input data:
# historical btf: for timing
# data for N past btf (N.btf) vs. eos regression (fitted as lm(log(N.btf) ~ log(eos)))
source("Data_prep.R")

# view data
# N.eos.dat
# head(hist.btf)

# calculate times
sec_prior_like = calc_prior_like_time()
sec_per_mcmc = calc_mcmc_time()

# save(list = ls(), file = "StartUpWorkspace")
# load(file = "StartUpWorkspace")

##### UI() #####
ui = navbarPage(strong("In-season Chinook Salmon Bayesian Risk Assessment Tool"), 
                windowTitle = "Bayes' Tool",
                theme = shinytheme("cerulean"), 

                # ESTIMATION TAB
                navbarMenu(menuName = "ests",
                  "Estimation", icon = icon("calculator", "fa-1x"),
                  tabPanel("Single-Day",
                           useShinyjs(), 
                           sidebarLayout(
                             sidebarPanel(width = 5,
                                          
                                          h4(strong("Run Size Forecast")),
                                          splitLayout(
                                            numericInput(inputId = "prior.mean", label = "Mean:", min = 0, max = 4e6, step = 5000, value = fcst_mean_def),
                                            hidden(
                                              div(id = "prior_cv_id",
                                                  numericInput(inputId = "prior.cv", label = "CV:", min = 0, max = 3, step = 0.01, value = fcst_cv_def)
                                              )
                                            )
                                          ),
                                          a(id = "toggle_prior_cv", "Show/Hide Uncertainty"),  
                                          
                                          
                                          h4(strong("In-Season Index Information")),
                                          splitLayout(
                                            textInput(inputId = "cdate", label = "Today's Date:", value = cdate_def),
                                            numericInput(inputId = "ccpue", label = "BTF CCPUE:", min = 0, max = 5000, step = 1, value = ccpue_def)
                                          ),
                                          
                                          h4(strong("Total Harvest Downstream of BTF")),
                                          splitLayout(
                                            numericInput(inputId = "ds.harv", label = "Mean:", min = 0, max = 100000, step = 1000, value = ds_harv_mean_def),
                                            hidden(
                                              div(id = "ds_harv_cv_id",
                                                  numericInput(inputId = "ds.harv.cv", label = "CV:", min = 0, max = 1, step = 0.01, value = ds_harv_cv_def)
                                              )  
                                            )
                                          ),
                                          
                                          
                                          a(id = "toggle_timing", 
                                            p("Show/Hide Included Timing Scenarios", popify(icon("info-circle"), title = NULL,
                                                                                            content = "Only exclude timing scenarios if you have a really good reason, this setting can have a big influence on the results of the estimation.", placement = "right")
                                              
                                              
                                            )),
                                          
                                          hidden(
                                            div(id = "timing_id",
                                                checkboxGroupInput(inputId = "timing",label = NULL,
                                                                   choices = list("Early" = 1, "Average" = 2, "Late" = 3),
                                                                   selected = 1:3, inline = T)
                                            )
                                          ),
                                          
                                          
                                          # mcmc settings
                                          h4(strong("Information Update Settings")),
                                          splitLayout(cellWidths = c("35%", "65%"),
                                                      checkboxInput(inputId = "do.mcmc", label = "Perform Update", value = T),
                                                      hidden(
                                                        div(id = "toggle_mcmc_def",
                                                            checkboxInput(inputId = "show_mcmc_opts", label = "Display/Change Settings", value = F)
                                                        )
                                                      )
                                          ),
                                          
                                          
                                          hidden(
                                            div(id = "mcmc_opts",
                                                splitLayout(cellWidths = c("30%", "30%", "20%", "20%"),
                                                            numericInput(inputId = "ni", label = "Samples:", min = 0, max = 1e6, step = 5000, value = ni_def),
                                                            numericInput(inputId = "nb", label = "Burn-in:", min = 0, max = 1e6, step = 5000, value = nb_def),
                                                            numericInput(inputId = "nc", label = "Chains:", min = 1, max = 4, step = 1, value = nc_def),
                                                            numericInput(inputId = "prop.sig", label = "Tune:", min = 0.01, max = 1, step = 0.05, value = prop_sig_def)
                                                )
                                            )
                                          ),
                                          
                                          hidden(
                                            div(id = "reset_mcmc_opts_id",
                                                # actionButton(inputId = "reset_mcmc_opts", label = "Restore Defaults", icon = icon("undo", "fa-1x"))
                                                a(id = "reset_mcmc_opts", p(icon("undo", "fa-1x"), "Restore Default Settings"))
                                            )
                                          ),
                                          
                                          uiOutput("time_text"),
                                          
                                          
                                          actionButton(inputId = "EstUpdate", label = "Calculate", icon = icon("cogs", "fa-1x")),
                                          actionButton(inputId = "EstReset", label = "Clear", icon = icon("trash", "fa-1x")),
                                          hr(), 
                                          dropdownButton(inputId = "export_dropdown", circle = F, status = "primary",
                                                         icon = icon("download"), label = "Export Estimates", 
                                                         
                                                         textInput(inputId = "export_suffix", label = NULL, placeholder = "Optional File Suffix"),
                                                         radioButtons(inputId = "export_ext", label = NULL, 
                                                                      choices = list("CSV" = ".csv", "TXT" = ".txt"), selected = ".csv",
                                                                      inline = T),
                                                         uiOutput("export_FileName"),
                                                         downloadLink(outputId = "DL_Ests", label = "Export File")
                                          )
                                          
                             ),
                             
                             mainPanel(width = 7,
                                       tabsetPanel(type = "tabs",
                                                   tabPanel(p(icon("area-chart", "fa-1x"), 
                                                              strong("Distributions"),
                                                              style = "font-size:16px"),
                                                            plotOutput("EstPlot1"),
                                                            br(),
                                                            uiOutput("DL_Dist_Plot_Button")),
                                                   tabPanel(p(icon("table", "fa-1x"), 
                                                              strong("Table"),
                                                              style = "font-size:16px"),
                                                            tableOutput("EstTable"),
                                                            br(),
                                                            uiOutput("DL_EstTable_Button"),
                                                            br()),
                                                   tabPanel(p(icon("bar-chart-o", "fa-1x"), 
                                                              strong("Components"),
                                                              style = "font-size:16px"), 
                                                            plotOutput(outputId = "EstPlot2"),
                                                            br(),
                                                            plotOutput(outputId = "EstPlot3"),
                                                            br(),
                                                            plotOutput(outputId = "EstPlot4")),
                                                   tabPanel(p(icon("thermometer-half", "fa-1x"),
                                                              strong("Diagnostics"),
                                                              style = "font-size:16px"),
                                                            plotOutput(outputId = "MCMCPlot"),
                                                            br(),
                                                            tableOutput(outputId = "MCMCTable"),
                                                            br(),
                                                            a(id = "toggle_MCMC_help", "Show/Hide More Information"),
                                                            
                                                            
                                                            hidden(
                                                              div(id = "MCMC_Diag_Help",
                                                                  wellPanel(
                                                                    h4(strong("What Are These Plots?")),
                                                                    p("Density Plot", popify(icon("info-circle"), title = NULL,
                                                                                             content = "This is a plot showing the frequency of sampling of each value of N (run size). You want to see the different lines drawn to be very similar to each other.", placement = "right")),
                                                                    
                                                                    p("Trace Plot", popify(icon("info-circle"), title = NULL,
                                                                                           content = "This is a plot time series of sampling for each value of N (run size). You want to see lines look like random jumps around the same average.", placement = "right")),
                                                                    
                                                                    h4(strong("What Are These Numbers?")),
                                                                    
                                                                    p("Saved Samples", popify(icon("info-circle"), title = NULL,
                                                                                              content = "This is the total number of posterior samples saved.", placement = "right", trigger = "hover")),
                                                                    
                                                                    p("BGR", popify(icon("info-circle"), title = NULL,
                                                                                    content = "This quantity provides an indication if the algorithm has converged to sampling from its target distribution. BGR value < 1.1 are generally indicative of convergence (i.e., good).", placement = "right", trigger = "hover"),
                                                                      a("See the Technical Details", href = "http://www.stat.columbia.edu/~gelman/research/published/brooksgelman2.pdf")),
                                                                    
                                                                    p("Effective Samples", popify(icon("info-circle"), title = NULL,
                                                                                                  content = "This quantity provides an estimate of the independent number of samples saved. Approximately 4,000 effective samples are required to precisely estimate the location of various quantiles in the tails.", placement = "right", trigger = "hover"),
                                                                      a("See the Technical Details", href = "http://people.ee.duke.edu/~lcarin/raftery92how.pdf")),
                                                                    
                                                                    p("Acceptance Rate", popify(icon("info-circle"), title = NULL,
                                                                                                content = "This is the fraction of all proposed values that were accepted in the MCMC algorithm. It has to do with the efficiency of the algorithm. If rate << 0.2, decrease the tune parameter, if rate >> 0.4, increase the tune parameter.", placement = "right", trigger = "hover"),
                                                                      a("See the Technical Details", href = "http://probability.ca/jeff/ftpdir/mylene2.pdf"))
                                                                  )
                                                              )
                                                            )
                                                            
                                                   )
                                       )
                                       
                             )
                             
                           )
                           
                  ),
                  tabPanel(
                    "Multi-Day",
                    fluidRow(
                      column(width = 3,
                             wellPanel(
                               h4(icon("file", "fa-1x"), strong("Select the Input File")),
                               fileInput(inputId = "MEst_in_file", label = NULL, multiple = F),
                               actionButton(inputId = "MEst_read", label = "Import", icon = icon("file", "fa-1x")),
                               actionButton(inputId = "MEst_start_over", label = "Start Over", icon = icon("refresh", "fa-1x")),
                               br(),
                               a(id = "toggle_example_file", "Show/Hide Example Input File"),
                               uiOutput("MEst_read_message")
                             )
                      ),
                      
                      column(width = 3,
                             hidden(
                               div(id = "MEst_dates",
                                   wellPanel(
                                     h4(icon("calendar", "fa-1x"), strong("Select Dates to Estimate")),
                                     splitLayout(
                                       textInput(inputId = "MEst_fdate", label = "First Date:"),
                                       textInput(inputId = "MEst_ldate", label = "Last Date:")
                                     ),
                                     uiOutput("MEst_date_message"),
                                     selectizeInput(inputId = "MEst_freq", 
                                                    label = p("Interval to Produce Estimates",
                                                              popify(el = icon("info-circle", "fa-1x"), title = NULL,
                                                                     content = "1 is every day in range, 2 is every other day, etc. The last day provided in the range will always be estimated")), 
                                                    choices = 1:10, selected = 1, multiple = F),
                                     
                                     actionButton(inputId = "MEst_select_dates", label = "Select Dates", icon = icon("check", "fa-1x")),
                                     
                                     
                                     uiOutput("MEst_selected_dates")
                                   )
                               )
                             )
                      ),
                      
                      column(width = 6,
                             hidden(
                               div(id = "MEst_settings",
                                   wellPanel(
                                     h4(icon("gear", "fa-1x"), strong("Estimation Settings")),
                                     h5(strong("Run Size Forecast")),
                                     splitLayout(
                                       numericInput(inputId = "MEst_prior_mean", label = "Mean:", min = 0, max = 4e6, step = 5000, value = fcst_mean_def),
                                       hidden(
                                         div(id = "MEst_prior_cv_id",
                                             numericInput(inputId = "MEst_prior_cv", label = "CV:", min = 0, max = 3, step = 0.01, value = fcst_cv_def)
                                         )
                                       )
                                     ),
                                     a(id = "toggle_MEst_prior_cv", "Show/Hide Uncertainty"),
                                     
                                     a(id = "toggle_MEst_timing", 
                                       p("Show/Hide Timing Scenarios", popify(icon("info-circle"), title = NULL,
                                                                              content = "Only exclude timing scenarios if you have a really good reason, this setting can have a big influence on the results of the estimation.", placement = "right")
                                         
                                         
                                       )),
                                     
                                     hidden(
                                       div(id = "MEst_timing_id",
                                           checkboxGroupInput(inputId = "MEst_timing",label = NULL,
                                                              choices = list("Early" = 1, "Average" = 2, "Late" = 3),
                                                              selected = 1:3, inline = T)
                                       )
                                     ),
                                     
                                     splitLayout(cellWidths = c("35%", "65%"),
                                                 checkboxInput(inputId = "MEst_do.mcmc", label = "Perform Update", value = T),
                                                 hidden(
                                                   div(id = "toggle_MEst_mcmc_def",
                                                       checkboxInput(inputId = "show_MEst_mcmc_opts", label = "Display/Change Settings", value = F)
                                                   )
                                                 )
                                     ),
                                     
                                     
                                     hidden(
                                       div(id = "MEst_mcmc_opts",
                                           splitLayout(cellWidths = c("30%", "30%", "20%", "20%"),
                                                       numericInput(inputId = "MEst_ni", label = "Samples:", min = 0, max = 1e6, step = 5000, value = ni_def),
                                                       numericInput(inputId = "MEst_nb", label = "Burn-in:", min = 0, max = 1e6, step = 5000, value = nb_def),
                                                       numericInput(inputId = "MEst_nc", label = "Chains:", min = 1, max = 4, step = 1, value = nc_def),
                                                       numericInput(inputId = "MEst_prop.sig", label = "Tune:", min = 0.01, max = 1, step = 0.05, value = prop_sig_def)
                                           )
                                       )
                                     ),
                                     
                                     hidden(
                                       div(id = "reset_MEst_mcmc_opts_id",
                                           # actionButton(inputId = "reset_mcmc_opts", label = "Restore Defaults", icon = icon("undo", "fa-1x"))
                                           a(id = "reset_MEst_mcmc_opts", p(icon("undo", "fa-1x"), "Restore Default Settings"))
                                       )
                                     ),
                                     
                                     uiOutput("MEst_time_text"),
                                     
                                     hr(),
                                     
                                     h4(icon("gear", "fa-1x"), strong("Output File Settings")),
                                     textInput(width = "300px", inputId = "zip_suffix", label = NULL, placeholder = "Optional File Suffix"),
                                     uiOutput("MEst_filename"),
                                     radioButtons(inputId = "MEst_ext", label = "Format of Individual Files", choices = c("CSV" = ".csv", "TXT" = ".txt"), selected = ".csv", inline = T),
                                     
                                     actionButton(inputId = "MEst_update", label = "Calculate", icon = icon("cogs", "fa-1x")),
                                     downloadButton("DL_MEst_zip", "Download Zip File")
                                     
                                     
                                   )
                                   
                               )
                             )
                      )
                      
                    ),
                    
                    fluidRow(
                      column(width = 5, 
                             hidden(
                               div(id = "ex_file",
                                   wellPanel(
                                     h4(icon("file", "fa-1x"), strong("Example Input File")),
                                     helpText(em("Must be saved as either a .csv or .txt file and formatted exactly as shown below")),
                                     helpText(tags$ul(
                                       tags$li(strong("date:"), "the date corresponding to the data for that row"),
                                       tags$li(strong("ccpue:"), "the BTF cumulative CPUE on that date; does not need to be an integer (whole number) value. Zeros are allowed."),
                                       tags$li(strong("ds_charv:"), "cumulative harvest to date downstream of the BTF site"),
                                       tags$li(strong("ds_charv_cv:"), "the coefficient of variation for ds_charv")
                                     )),
                                     helpText("These data are from 2017, if you wish, you can download this file as a template.", em("The harvest values are approximate.")),
                                     downloadButton("DL_ex_file", "Download This File"),
                                     br(), br(),
                                     tableOutput("example_file")
                                   ))
                             )
                      )
                      
                    )
                    )
                  ),
                
                # RISK ANALYSIS TAB
                navbarMenu("Risk Analysis", icon = icon("filter", "fa-1x"),
                           useShinyjs(),
                           tabPanel("Choose Harvest Objective",
                                    sidebarLayout(
                                      sidebarPanel(width = 5,
                                        div(id = "RA1_message",
                                            h4("Please complete the", icon("calculator", "fa-1x"), strong("Single-Day Estimation"), "option before performing risk analysis")
                                            ),
                                        hidden(
                                        div(id = "RA1",
                                        h4(strong("Run Size Estimate(s)")),
                                        uiOutput("PlcyInfoChoice1"),
                                        
                                        h4(strong("Current Total Harvest")),
                                        splitLayout(
                                          numericInput(inputId = "harv", label = "Mean:", min = 0, max = 100000, step = 1000, value = harv_mean_def),
                                          hidden(
                                            div(id = "harv_cv_id",
                                                numericInput(inputId = "harv.cv", label = "CV:", min = 0, max = 1, step = 0.01, value = harv_cv_def)
                                            )  
                                          )
                                        ),
                                        
                                        a(id = "toggle_Hcand_lim", "Show/Hide Range of Considered Harvest Candidates"),
                                        hidden(
                                          div(id = "Hcand_lim_id", 
                                              sliderInput(inputId = "Hcand_lim", label = NULL, ticks = F,
                                                          min = 0, max = 300000, value = H_cand_lim_def, step = 5000)
                                          )
                                        ),
                                        
                                        h4(strong("Escapement Target")),
                                        sliderInput(inputId = "Sobj", label = NULL, ticks = F,
                                                    min = 65000, max = 120000, value = Sobj_def, step = 500),
                                        
                                        h4(strong("Risk Tolerance")),
                                        uiOutput("p_star_help"),
                                        uiOutput("p_star_slide"),
                                        
                                        h5(strong("Frame Risk As:")),
                                        radioButtons(inputId = "risk_direction", label = NULL,
                                                     choices = list("Pr(S < Target)" = 1,
                                                                    "Pr(S > Target)" = 2), selected = 1, inline = T),
                                        
                                        miniButtonBlock(
                                          actionButton(inputId = "Plcy1Update", label = "Update", icon = icon("cogs", "fa-1x")),
                                          actionButton(inputId = "Plcy1Reset", label = "Clear", icon = icon("trash", "fa-1x")),
                                          downloadButton(outputId = "DL_Plcy_Plot", label = "Plot"),
                                          downloadButton(outputId = "DL_Plcy_Table", "Table")
                                        )
                                        
                                      
                                      ))),
                                      
                                      mainPanel(width = 7,
                                        plotOutput("PolicyPlot"),
                                        br(),
                                        uiOutput("DL_Plcy_Plot_Button"),
                                        br(),
                                        tableOutput("PolicyTable"),
                                        uiOutput("DL_Plcy_Table_Button"),
                                        br()
                                      )
                                    )
                           ),
                           
                           tabPanel("Compare Harvest Objectives",
                                    sidebarLayout(
                                      sidebarPanel(width = 5,
                                       div(id = "RA2_message",
                                           h4("Please complete the", icon("calculator", "fa-1x"), strong("Single-Day Estimation"), "option before performing risk analysis")
                                       ),
                                        hidden(
                                          div(id = "RA2",
                                        h4(strong("Run Size Estimate")),
                                        uiOutput("PlcyInfoChoice2"),
                                        
                                        h4(strong("Current Total Harvest")),
                                        
                                        splitLayout(
                                          numericInput(inputId = "harv2", label = "Mean:", min = 0, max = 100000, step = 1000, value = harv_mean_def),
                                          hidden(
                                            div(id = "harv2_cv_id",
                                                numericInput(inputId = "harv2.cv", label = "CV:", min = 0, max = 1, step = 0.01, value = harv_cv_def)
                                            )  
                                          )
                                        ),
                                       
                                        h4(strong("Candidate Harvest Objectives")),
                                        helpText("In addition to harvest currently taken listed above"),
                                        splitLayout(
                                          numericInput("Hobj1", "H1:", value = 0, min = 0, max = 1.5e5, step = 1000),
                                          numericInput("Hobj2", "H2:", value = 0, min = 0, max = 1.5e5, step = 1000),
                                          numericInput("Hobj3", "H3:", value = 0, min = 0, max = 1.5e5, step = 1000)
                                        ),
                                        actionButton("Plcy2Update", "Update", icon = icon("cogs", "fa-1x")),
                                        actionButton("Plcy2Reset", "Clear", icon = icon("trash", "fa-1x")),
                                        downloadButton("DL_PlcyTable2", "Download Table")
                                      ))),
                                      mainPanel(width = 7,
                                        tableOutput("PlcyTable2")
                                      )
                                    )
                           )
                ),
                
                # TRACK LEARNING TAB
                tabPanel("Track Learning", icon = icon("line-chart", "fa-1x"),
                         sidebarLayout(
                           sidebarPanel(width = 5,
                                        h4(icon("file", "fa-1x"), strong("Select Files Containing Desired Estimates"), 
                                           popify(icon("question-circle", "fa-1x"),title = NULL,
                                                  content = "These are the individual output files from downloading estimates from the estimation tabs. If the individual files are in a .zip file, you must paste them into a new folder before importing")),
                                        # helpText(),
                                        fileInput(inputId = "learn_files", label = NULL, multiple = T, placeholder = "Select 2 or More Files"),
                                        actionButton("read_learn", "Import", icon = icon("file", "fa-1x")),
                                        actionButton("learn_start_over", "Start Over", icon = icon("refresh", "fa-1x")),
                                        uiOutput("learn_read_text"),
                                        uiOutput("learn_read_text2"),
                                        
                                        hidden(
                                          div(id = "learn_opts",
                                              h4(icon("calendar", "fa-1x"), strong("Dates to View")),
                                              splitLayout(
                                                textInput("learn_fdate", label = "First Date:"),
                                                textInput("learn_ldate", label = "Last Date:")
                                              ),
                                              
                                              h4(strong("Quantity to View")),
                                              selectizeInput(inputId = "learn_est", label = NULL,
                                                             choices = list("Forecast" = "prior","BTF Only" = "likelihood","Updated" = "posterior","Run Timing" = "p","EOS BTF CCPUE" = "eos"),
                                                             selected = 1),
                                              
                                              h4(strong("Y-Axis Range")),
                                              uiOutput("LearnSlide"),
                                              
                                              h4(strong("Legend")),
                                              selectizeInput(inputId = "legend_loc", label = NULL,
                                                             choices = list("None" = "none", "Top Left" = "topleft", "Top Right" = "topright", "Bottom Left" = "bottomleft", "Bottom Right" = "bottomright"),
                                                             selected = "none"),
                                              
                                              actionButton("learn_update", label = "Update", icon = icon("refresh", "fa-1x")),
                                              actionButton("learn_clear", label = "Clear", icon = icon("trash", "fa-1x")),
                                              downloadButton("DL_learn_plot", label = "Download Plot")
                                          )
                                        )
                                        
                           ),
                           
                           mainPanel(width = 7,
                                     plotOutput("LearnPlot")
                           )
                         )
                ),
                
                navbarMenu(
                  "About",icon = icon("info-circle"),
                  tabPanel(
                    p(icon("info-circle"), "Overview", style = "margin:0;"),
                    fluidRow(
                      column(
                        width = 8,
                        wellPanel(
                          h4(strong("Overview")),
                          p("This tool was developed to allow users to easily perform run size estimation and evaluate
                            the consequences of different harvest alternatives in a risk analysis framework.",
                            strong("As currently constructed, the tool performs these tasks for Chinook salmon in the Kuskokwim River only, 
                                   and does not assess risk of failing to meet tributary-level escapement goals.")),
                          p("The tool updates run size estimates from the pre-season forecast using in-season data from the Bethel Test Fishery
                            using a statistical framework called", em("Bayesian Inference."), "Bayesian Inference provides updated probabilities
                            of the different run size outcomes managers are considering after obtaining new information, which can then
                            be used to assess the probability of meeting or failing to meet various levels of drainage-wide escapement if different numbers of fish were harvested."),
                          
                          
                          
                          h5(strong("What is Bayesian Inference?")),
                          p("This tool uses a statistical method known as", em("Bayesian Inference"), "to perform the update.
                            Bayesian Inference combines two sources of information to obtain an updated understanding of the value of a quantity we are uncertain about (the current year's run size is the uncertain quantity in this case). 
                            The first peice of information is called the ", em("prior knowledge"), "and is expressed as the probability of each possible run size", em("before observing new information about that quantity."),
                            "In this case, the pre-season forecast is the prior knowledge about the run size, and it can be expressed probabilistically (e.g., there is an X% chance the run will be smaller than Y).", 
                            "The second peice of information is known as ", em("new data"), "and in our case it is the sum of BTF catches through the current day of the run.",
                            "Bayesian Inference uses the new data to update the prior understanding to obtain an updated understanding using the laws of probability. This new understanding is referred to as", em("posterior knowledge"), "but to avoid confusion, in this tool it is referred to as ", em("updated knowledge."),
                            "A key property of Bayesian Inference is that the information source that is most precise (i.e., least uncertain) will contribute the most to updated knowledge, thus Bayesian Inference is a variance-based approach to weighting two sources of information."),
                          
                          h5(strong("What is MCMC?")),
                          p("MCMC stands for Markov Chain Monte Carlo simulation, and it is a commonly-used family of numerical algorithms we use to implement Bayesian Inference. In most Bayesian problems, the calculations are too complex to be performed analytically (i.e., by solving equations with alegbra or calculus), which is why MCMC methods must be used.
                            It is because of MCMC that this tool takes some time to perform the update. You can find more details on MCMC in the technical documentation under 'Posterior Formulation'."
                          )
                          
                          
                          )
                          )
                          )
                  ),
                  tabPanel(
                    p(icon("question-circle"), "User Manual", style = "margin:0;"),
                    fluidRow(
                      column(
                        width = 6,
                        wellPanel(
                          h4(strong("User Manual")),
                          p("We strongly encourage users of the tool with any questions about how to use this tool to consult the", a("user manual webpage.", href = "https://bstaton.shinyapps.io/BayesTool_UserMan/"),
                            "The user manual was developed specifically to be a walkthrough tutorial for how to use this tool, and we hope users will find it informative."),
                          p("If your questions were not answered by the user manual, please see the", strong(icon("envelope"), "Contact Page"),
                            "to send us an email. After answering your question(s), we will consider updating the user manual to help future users."
                          )
                        )
                      )
                    )
                  ),
                  tabPanel(
                    p(icon("file"), "Technical Documentation", style = "margin:0;"),
                    fluidRow(
                      column(
                        width = 6,
                        wellPanel(
                          h4(strong("Technical Documentation")),
                          p("Should users be interested in the details of how this tool works, the motivations for its development, and other technical information, we encourage them to visit the", a("technical documentation webpage", href = "https://bstaton.shinyapps.io/BayesTool_TechDoc/"),
                            "devoted to this purpose. Though we caution, much of the webpage is devoted to the statistical details of the tool, and many users will not find it light reading material."),
                          p("If you have questions or concerns about the information found in the technical documentation, please visit the ", strong(icon("envelope"), "Contact Page"),
                            "to send us an email."
                          )
                        )
                      )
                    )
                  ),
                  tabPanel(
                    p(icon("calendar"), "Historical Data", style = "margin:0;"),
                    fluidRow(
                      column(
                        width = 8,
                        wellPanel(
                          h4(strong("Historical Data")),
                          p("This tab is provided so users can investigate the behavior of the tool across different types of years.
                             It provides all of the information needed to use the tool pretending it had been available going back to 2008.
                             Click the button below to download all of the information in a single document."),
                          downloadButton("DL_HistData", "Download the Historical Data")
                        )
                      )
                    )
                  ),
                           
                  tabPanel(p(icon("users"), "Developers", style = "margin:0;"),
                           fluidRow(
                             column(
                               width = 8,
                               wellPanel(
                                 h4(strong("Tool Developers")),
                                 p("The interface of this tool was developed and is maintained by Ben Staton. 
                                   Ben is a Graduate Student Researcher and USFWS Pathways Quantitative Ecologist,
                                   who has studied the statistical aspects of Chinook salmon assessment
                                   and management issues in the Kuskokwim River since 2014. 
                                   The statistical framework for the tool was developed jointly by Ben and his advisor, Dr. Matt Catalano,
                                   who is the Principal Investigator of the", a("Quantitative Fisheries Laboratory", href = "http://sfaas.auburn.edu/catalano/quantitative-fisheries-lab/"),
                                   "at Auburn University and has studied salmon assessment methods in Alaska since 2009."),
                                 
                                 h4(strong("Contributors")),
                                 p("Many people have given incredibly useful feedback on this tool and on its utility for aiding in management discussions (listed alphabetically):",
                                   tags$ul(
                                     tags$li(strong("Dr. Bill Bechtol,"), "statistical consultant for the Kuskokwim River Inter-tribal Fisheries Commission"),
                                     tags$li(strong("Dr. Lew Coggins,"), "biology program supervisor, U.S. Fish and Wildlife Service, Yukon Delta National Wildlife Refuge"),
                                     tags$li(strong("Zach Liller,"), "Research Coordinator for the Arctic-Yukon-Kuskokwim Region, Alaska Department of Fish and Game, Commercial Fisheries Division"),
                                     tags$li(strong("Nick Smith,"), "Kuskokwim Area Fisheries Research Biologist, Alaska Department of Fish and Game, Commercial Fisheries Division"),
                                     tags$li(strong("Kuskokwim area fisheries managers,"), "from the Alaska Department of Fish and Game, Kuskokwim River Inter-tribal Fisheries Commission, and the U.S. Fish and Wildlife Service")
                                   )
                                 ),
                                 
                                 p("The statistical framework used in the tool has been evaluated and is currently undergoing peer-review for publication in the primary scientific literature."),
                                 p("The statistical aspects of the tool have undergone informal review by Alaska Department of Fish and Game biometric staff."),
                                 p("While the tool and its philosophy is novel to the Kuskokwim River Chinook salmon stock, it builds on pioneering work by Carl Walters, Sandy Buckingham, Stephen Fried, and Ray Hilborn done in the 1970s and 1980s.
                                   To our knowledge, the tool is the first to allow users to easily perform these types of relatively complex calculations without having to interact with spreadsheets or code."),
                                 p("The interface of this tool was developed using", a("Program R", href = "https://www.r-project.org/"), "which handles the statistical computing and an package called", a("Shiny", href = "https://shiny.rstudio.com/"),
                                   "which allows for integration of R code into the web-interface."),
                                 p("Funding for the development of this tool was provided by the Arctic-Yukon-Kuskokwim Sustainable Salmon Initiative, administered through the Bering Sea Fishermen's Association, through a research grant to Dr. Matt Catalano.")
                                 )
                             )
                           )
                             ),
                  tabPanel(p(icon("envelope"), "Contact", style = "margin:0;"),
                           fluidRow(
                             column(
                               width = 6,
                               wellPanel(
                                 h4(strong("Contact")),
                                 p("Should users have any questions about the tool that were not answered by
                                   the user manual or technical documentation, please feel free to contact the lead tool developer, Ben Staton, 
                                   at bas0041@auburn.edu"),
                                 p("We are also ready and willing to hear about any feedback, suggestions for future functionality, or bugs you may find.")
                               )
                             )
                           )
                  )
                  )
)


##### SERVER() #####

server = function(input, output, session) {
  
  ### EVERYTHING FOR SINGLE DAY ESTIMATION TAB ###
  # toggle the prior cv option
  onclick("toggle_prior_cv", toggle(id = "prior_cv_id", anim = F))
  
  # restore mcmc defaults
  onclick(id = "reset_mcmc_opts", reset("mcmc_opts"))
  
  # toggle MCMC help
  onclick(id = "toggle_MCMC_help", toggle(id = "MCMC_Diag_Help", anim = T))
  
  # toggle timing scenario options
  onclick(id = "toggle_timing", toggle(id = "timing_id", anim = T))
  
  # toggle other options on this tab
  observe({
    toggle(id = "ds_harv_cv_id", condition = input$ds.harv != 0 & !is.na(input$ds.harv))  # toggle the harvest cv input
    toggle(id = "toggle_mcmc_def", condition = input$do.mcmc, anim = F)                   # toggle the change settings checkbox
    toggle(id = "mcmc_opts", condition = input$show_mcmc_opts & input$do.mcmc, anim = T)  # toggle the mcmc options
    toggle(id = "reset_mcmc_opts_id", condition = {                                       # toggle the reset link
      !(input$ni == ni_def & input$nb == nb_def &
          input$nc == nc_def & input$prop.sig == prop_sig_def) & input$do.mcmc
    })
    toggleState("EstReset", condition = !is.null(out$EstOut))                             # toggle the reset button
    toggleState("export_dropdown", condition = !is.null(out$EstOut))
  })
  
  # export file name
  output$export_FileName = renderUI({
    helpText(icon("file"), paste("File Name:", create_export_filename(date = input$cdate, suffix = input$export_suffix, ext = input$export_ext)))
  })
  
  output$time_text = renderUI({
    
    if (input$do.mcmc) { # if doing mcmc
        mcmc_time = input$ni * input$nc * sec_per_mcmc # time is a function of defaults
    } else {
      mcmc_time = 0 # if not doing mcmc don't calculate this time
    }
    
    total_time = ceiling(sec_prior_like + mcmc_time)
    time_amount = ifelse(total_time == 1, "second", "seconds")
    
    helpText(icon("clock-o", "fa-1x"), paste("Estimated Time:", total_time, time_amount))
    
  })
  
  out = reactiveValues()
  
  observeEvent(input$EstUpdate, {
    starttime = Sys.time()
    # fit the historical regression
    fit.eos = fit.reg.mod(dat = N.eos.dat)
    
    hist.p = filter(hist.btf, date == input$cdate & rt.type %in% input$timing) %>% select(p.ccpue) %>% unlist %>% unname
    # hist.p = filter(hist.btf, date == input$cdate) %>% select(p.ccpue) %>% unlist %>% unname
    rt.mu = estBetaParams(mean(hist.p), var(hist.p))
    withProgress(message = "Drawing random samples:", value = 0, {
      
      # sample the prior
      incProgress(1/2, detail = "Forecast")
      prior.samp = exp(rnorm(n.boot, log(input$prior.mean) - 0.5 * cv2sig(input$prior.cv)^2, cv2sig(input$prior.cv)))
      
      # generate N samples based on BTF
      incProgress(2/2, detail = "BTF Expansion")
      boot_out = sample.likelihood(n.boot = n.boot,
                                   obs.ccpue = input$ccpue,
                                   fit.eos = fit.eos,
                                   rt.method = "beta",
                                   rt.mu = rt.mu, 
                                   sub.harv = input$ds.harv, 
                                   sub.cv = input$ds.harv.cv,
                                   com.harv = 0
      )
    })
    
    if (input$do.mcmc) {
      # obtain pdf's
      prior.fun = approxfun(density(prior.samp, from = 0, to = dens.max))
      like.fun = approxfun(density(boot_out[,"N"], from = 0, to = dens.max))
      
      withProgress(message = "Running MCMC:", value = 0, {
        # run the mcmc
        mcmc.out = lapply(1:input$nc, FUN = function(i) {
          incProgress((i-1)/input$nc, detail = paste("Chain #", i))
          MH(init = runif(1, 1e5, 3e5),
             ni = input$ni, nb = input$nb, nt = 1, prop.sig = input$prop.sig, 
             like.fun = like.fun, prior.fun = prior.fun)
        }) # lapply
      }) # progress
    } # if do.mcmc
    
    elapsed = Sys.time() - starttime
    
    # make the dataframe to export
    prior_summ = summ(prior.samp)
    like_summ = summ(boot_out[,"N"])
    harv_summ = summ(boot_out[,"sub.harv"])
    p_summ = summ(boot_out[,"p"])
    eos_summ = summ(boot_out[,"eos"])
    
    if (input$do.mcmc) {
      post_summ = summ(unlist(lapply(mcmc.out, function(x) x$post.samp)))
    } else {
      post_summ = rep(NA, 9)
    }
    
    q = names(prior_summ)
    
    export = data.frame(date = input$cdate,
               source = rep(c("prior", "likelihood", "posterior", "p", "eos", "ds_harv"), each = length(q)),
               stat = rep(q, 6),
               value = c(prior_summ, like_summ, post_summ, p_summ, eos_summ, harv_summ))
    
    
    # make the output list
    if (input$do.mcmc) {
      EstOut = list(
        did_mcmc = T,
        boot_out = boot_out,
        prior_samp = prior.samp,
        post_samp = matrix(unlist(lapply(mcmc.out, function(x) x$post.samp)),  input$ni - input$nb, input$nc),
        accept = unlist(lapply(mcmc.out, function(x) x$accept.rate)),
        cdate = input$cdate,
        mcmc_settings = c(ni = input$ni, nb = input$nb, nc = input$nc, tune = input$prop.sig),
        elapsed = elapsed,
        export = export
      )
    } else {
      EstOut = list(
        did_mcmc = F,
        boot_out = boot_out,
        prior_samp = prior.samp,
        cdate = input$cdate,
        elapsed = elapsed,
        export = export
      )
    }
    out$EstOut = EstOut
  })
  
  observeEvent(input$EstReset, {
    out$EstOut = NULL
  })
  
  output$DL_Ests = downloadHandler(
    filename = function() {
      create_export_filename(date = input$cdate, suffix = input$export_suffix, ext = input$export_ext)
      },
    content = function(file) {
      if (input$export_ext == ".csv") {
        write.csv(out$EstOut$export, file, row.names = F)
      } else {
        write.table(out$EstOut$export, file, row.names = F)
      }
    }
  )
  
  ### PLOTS/TABLES FOR ESTIMATION TAB ###
  output$EstPlot1 = renderPlot({
    create_dist_plot(out)
  })
  
  output$DL_Dist_Plot_Button = renderUI({
    if (!is.null(out$EstOut$boot_out)) {
      downloadButton(outputId = "DL_Dist_Plot", label = "Download Plot")
    } else {
      NULL
    }
  })
  
  output$DL_Dist_Plot = downloadHandler(
    filename = function() {
      paste("DistPlot_", paste(unlist(str_split(input$cdate, pattern = "/")), collapse = "_"), ".png", sep = "")
    },
    content = function(file) {
      ppi = 600
      png(file, h = 5 * ppi, w = 7 * ppi, res = ppi)
      create_dist_plot(out)
      dev.off()
    }
  )
  
  output$EstPlot2 = renderPlot({
    if (!is.null(out$EstOut$boot_out)) {
      # if (F) {
      par(mar = c(4,5,2,1.5), xaxs = "i", yaxs = "i", cex.axis = 1.5, cex.main = 1.5)
      rand.i = sample(x = 1:n.boot, size = round(n.boot/3), replace = F)
      ymax = max(hist(out$EstOut$boot_out[rand.i,"p"], breaks = seq(0, 1, 0.01), plot = F)$counts)
      hist(out$EstOut$boot_out[rand.i,"p"], yaxt = "n", xaxt = "n",
           col = "skyblue2", xlim = c(0,1), ylim = c(0, ymax) * 1.05, breaks = seq(0,1, 0.01),
           main = "% of Run Complete", xlab = "", ylab = "")
      axis(side = 1, at = seq(0, 1, 0.2), labels = paste(seq(0, 1, 0.2) * 100, "%", sep = ""), lwd = 4)
      mtext(side = 2, line = 2, "Relative Probability", cex = 1.4)
      box(lwd = 4)
    }
  })
  
  output$EstPlot3 = renderPlot({
    if (!is.null(out$EstOut$boot_out)) {
      # if (F) {
      par(mar = c(4,5,2,1.5), xaxs = "i", yaxs = "i", cex.axis = 1.5, cex.main = 1.5)
      # rand.i = sample(x = 1:n.boot, size = round(n.boot/3), replace = F)
      hist_samps = out$EstOut$boot_out[,"eos"]
      rand.i = which(hist_samps <= 1500)
      ymax = max(hist(hist_samps[rand.i], breaks = seq(0, 1500, length = 100), plot = F)$counts)
      hist(hist_samps[rand.i], yaxt = "n", col = "skyblue2",
           xlim = c(0, max(1501)), ylim = c(0, ymax) * 1.05, breaks = seq(0, 1500, length = 100),
           main = "Predicted EOS BTF CCPUE", xlab = "", ylab = "", xaxt = "n")
      axis(side = 1, lwd = 4)
      mtext(side = 2, line = 2, "Relative Probability", cex = 1.4)
      box(lwd = 4)
    }
  })
  
  output$EstPlot4 = renderPlot({
    if (!is.null(out$EstOut$boot_out)) {
      # if (F) {
      par(mar = c(4,5,2,1.5), cex.main = 1.5, xaxs = "i", yaxs = "i", cex.axis = 1.5, cex.lab = 1.4)
      rand.i = sample(x = 1:n.boot, size = 500, replace = F)
      plot(out$EstOut$boot_out[rand.i,"N.btf"] ~ out$EstOut$boot_out[rand.i,"eos"],
           pch = 16, cex = 1.8, col = alpha("skyblue3", 0.5),
           xlab = "Predicted EOS BTF CCPUE", ylab = "",
           xlim = c(0, 1500), ylim = c(0, 400000),
           main = "EOS BTF CCPUE Expansion", yaxt = "n", xaxt = "n")
      axis(side = 2, at = seq(0, 400000, 50000), labels = seq(0, 400, 50), las = 2, lwd = 4)
      axis(side = 1, at = seq(0, 1500, 200), labels = seq(0, 1500, 200), lwd = 4)
      mtext(side = 2, line = 3.5, "Run Size Past BTF (1000s)", cex = 1.4)
      fit.eos = fit.reg.mod(dat = N.eos.dat)
      new_eos = seq(0, 1500, 10)
      new_N = exp(fit.eos$coef[1] + fit.eos$coef[2] * log(new_eos))
      lines(new_N ~ new_eos, lwd = 3)
      
      box(lwd = 2, lwd = 4)
    }
  })
  
  output$EstTable = renderTable(align = "c", striped = T, hover = T, bordered = T, {
    if (!is.null(out$EstOut$boot_out)) {
      create_est_table(out)
    } else {
      NULL
    }
  })
  
  output$DL_EstTable_Button = renderUI({
    if (!is.null(out$EstOut$boot_out)) {
      downloadButton("DL_EstTable", "Download Table")
    } else {
      NULL
    }
  })
  
  output$DL_EstTable = downloadHandler(
    filename = function() {
      paste(paste("Est_Summary", paste(unlist(str_split(input$cdate, pattern = "/")), collapse = "_"), sep = "_"), ".png", sep = "")
    },
    content = function(file) {
      ppi = 600
      png(file, h = 3 * ppi, w = 5.5 * ppi, res = ppi)
      grid.table(create_est_table(out), rows = NULL)
      dev.off()
    }
  )
  
  output$MCMCPlot = renderPlot({
    if (!is.null(out$EstOut$post_samp)) {
      diag.plots(post.samp = out$EstOut$post_samp)
    }
  })
  
  output$MCMCTable = renderTable(align = "c", {
    if (!is.null(out$EstOut$post_samp)) {
      
      n_saved = prettyNum(length(as.numeric(out$EstOut$post_samp)), big.mark = ",", scientific = F)
      ess = unname(effectiveSize(as.mcmc.list(lapply(as.data.frame(out$EstOut$post_samp), mcmc))))
      ess = prettyNum(round(ess), big.mark = ",", scientific = F)
      accept = round(mean(out$EstOut$accept),2)
      if (out$EstOut$mcmc_settings["nc"] > 1) {
        bgr = round(unname(gelman.diag(x = as.mcmc.list(lapply(as.data.frame(out$EstOut$post_samp), mcmc)), multivariate = F)[[1]][1,1]),2)
      } else {
        bgr = NA
      }
      
      df = data.frame(n_saved = n_saved, ESS = ess, bgr = bgr, accept = accept)
      colnames(df) = c("Saved Samples", "Effective Samples", "BGR", "Acceptance Rate")
      df
    }
  })
  
  ### EVERYTHING FOR MULTI-DAY ESTIMATION ###
  # the example file to view
  output$example_file = renderTable(bordered = T, striped = T, {
    
    tab = read.csv("Inputs/Example_MultiDay_Input.csv", stringsAsFactors = F)
    tab
    
  })
  
  # download example file
  output$DL_ex_file = downloadHandler(
    filename = function() {
      "Example_MultiDay_Input.csv"
    },
    content = function(file) {
      tab = read.csv("Inputs/Example_MultiDay_Input.csv", stringsAsFactors = F)
      write.csv(tab, file, row.names = F)
    }
  )
  
  # toggle the Import button when file is uploaded
  observe({
    toggleState("MEst_read", condition = !is.null(input$MEst_in_file))
    toggleState("MEst_start_over", condition = !is.null(out$MEstOut))
    toggleState("DL_MEst_zip", condition = !is.null(out$MEstSaved))
    toggle("MEst_dates", condition = !is.null(out$MEstOut$dat), anim = T)
    toggle("MEst_settings", condition = !is.null(out$MEstDat), anim = T)
    toggle(id = "toggle_MEst_mcmc_def", condition = input$MEst_do.mcmc, anim = F)                   # toggle the change settings checkbox
    toggle(id = "MEst_mcmc_opts", condition = input$show_MEst_mcmc_opts & input$MEst_do.mcmc, anim = T)  # toggle the mcmc options
    toggle(id = "reset_MEst_mcmc_opts_id", condition = {                                       # toggle the reset link
      !(input$MEst_ni == ni_def & input$MEst_nb == nb_def &
          input$MEst_nc == nc_def & input$MEst_prop.sig == prop_sig_def) & input$MEst_do.mcmc
    })
  })
  
  onclick(id = "toggle_example_file", toggle("ex_file", anim = T))
  onclick("toggle_MEst_prior_cv", toggle("MEst_prior_cv_id"))
  onclick("toggle_MEst_timing", toggle("MEst_timing_id", anim = T))
  onclick(id = "reset_MEst_mcmc_opts", reset("MEst_mcmc_opts"))
  
  # read in the file when import button clicked
  observeEvent(input$MEst_read, {
    dat = read.csv(input$MEst_in_file$datapath, stringsAsFactors = F)
    
    if (all(names(dat) == c("date", "ccpue", "ds_charv", "ds_charv_cv"))) {
      message = "File read in successfully and is in the right format"
      dat = merge(dat, dates, by = "date")
      dat = dat[order(dat$day, decreasing = F),]
      fdate = dat[dat$day == min(dat$day),"date"]
      ldate = dat[dat$day == max(dat$day),"date"]
      updateTextInput(session = session, inputId = "MEst_fdate", value = fdate)
      updateTextInput(session = session, inputId = "MEst_ldate", value = ldate)
    } else {
      message = "File read in successfully, but is in the wrong format. Must have columns 'date', 'ccpue', 'ds_charv', and 'ds_charv_cv'. Cannot Proceed."
      dat = NULL
      fdate = NULL
      ldate = NULL
    } 
    out$MEstOut = list(message = message, dat = dat, fdate = fdate, ldate = ldate)
  })
  
  observeEvent(input$MEst_start_over, {
    reset(id = "MEst_in_file")
    out$MEstOut = NULL
    out$MEstDat = NULL
  })
  
  # message giving import status
  output$MEst_read_message = renderUI({
    if (!is.null(out$MEstOut$message)) {
      div(
        hr(),
        helpText(out$MEstOut$message)
      )
    } else {
      NULL
    }
  })
  
  # message giving input date error
  output$MEst_date_message = renderUI({
    if (!(input$MEst_fdate %in% out$MEstOut$dat$date) | !(input$MEst_ldate %in% out$MEstOut$dat$date)) {
      div(
        helpText("One of the dates provided is not within the range in the input file.
                 The first date available is", out$MEstOut$fdate, "and the last date available is", out$MEstOut$ldate),
        helpText("Please select valid dates within this range, otherwise the app will crash")
        )
    } else {
      NULL
    }
  })
  
  # filter dates when "Select Dates" button clicked
  observeEvent(input$MEst_select_dates, {
    dat = out$MEstOut$dat
    dat = dat[dat$day >= dat$day[dat$date == input$MEst_fdate],]
    dat = dat[dat$day <= dat$day[dat$date == input$MEst_ldate],]
    
    fday = min(dat$day)
    lday = max(dat$day)
    
    keep_days = seq(fday, lday, by = as.numeric(input$MEst_freq))
    if (keep_days[length(keep_days)] != lday) keep_days = c(keep_days, lday)
    
    out$MEstDat = dat[dat$day %in% keep_days,]
  })
  
  # message showing the selected dates to estimate
  output$MEst_selected_dates = renderUI({
    if (!is.null(out$MEstDat)) {
      div(
        helpText(strong(paste(nrow(out$MEstDat), "Dates Selected:"))),
        helpText(paste(out$MEstDat$date, collapse = ", "))
      )
    } else {
      NULL
    }
  })
  
  # time for calculations
  output$MEst_time_text = renderUI({
    if (input$MEst_do.mcmc) { # if doing mcmc
      mcmc_time = input$MEst_ni * input$MEst_nc * sec_per_mcmc  # time is a function of defaults
    } else {
      mcmc_time = 0 # if not doing mcmc don't calculate this time
    }
    
    total_time_seconds = ceiling(sec_prior_like + mcmc_time) * nrow(out$MEstDat)
    
    if (total_time_seconds > 60) {
      total_time = round(total_time_seconds/60, 1)
      time_amount = "minutes"
    } else {
      total_time = total_time_seconds
      time_amount = ifelse(total_time == 1, "second", "seconds")
    }
    
    helpText(icon("clock-o", "fa-1x"), paste("Estimated Time:", total_time, time_amount))
  })
  
  # preview file name for output zip file
  output$MEst_filename = renderUI({
    helpText(icon("file", "fa-1x"), "File Name:", 
             create_export_filename(date = c(input$MEst_fdate, input$MEst_ldate), suffix = input$zip_suffix, ext = ".zip"))
  })
  
  # perform the estimation
  observeEvent(input$MEst_update, {
    
    # fit the historical regression
    fit.eos = fit.reg.mod(dat = N.eos.dat)
    
    # sample the prior
    prior_samps = exp(rnorm(n.boot, log(input$MEst_prior_mean) - 0.5 * cv2sig(input$MEst_prior_cv)^2, cv2sig(input$MEst_prior_cv)))
    prior_summ = summ(prior_samps)
    prior_fun = approxfun(density(prior_samps, from = 0, to = dens.max))
    
    dat = out$MEstDat
    est_dates = dat$date
    
    # placeholder if not doing MCMC
    if (!input$MEst_do.mcmc) post_summ = rep(NA, 9)
    
    post_name = rep("posterior", 9)
    like_name = rep("likelihood", 9)
    prior_name = rep("prior", 9)
    p_name = rep("p", 9)
    eos_name = rep("eos", 9)
    harv_name = rep("ds_harv", 9)
    
    accept_rate = list()
    bgr = list()
    ess = list()
    ests = list()
    
    boot_out = array(NA, dim = c(n.boot, 4, length(est_dates)))
    dimnames(boot_out) = list(1:n.boot, c("p", "eos", "N", "sub.harv"), est_dates)
    withProgress(message = "Estimating Day:", value = 0, min = 0, max = length(est_dates), {
      for (d in 1:length(est_dates)) {
        incProgress(amount = d - 1, detail = est_dates[d])
        
        # sample likelihood
        hist.p = filter(hist.btf, date == est_dates[d] & rt.type %in% input$MEst_timing) %>% select(p.ccpue) %>% unlist %>% unname
        rt.mu = estBetaParams(mean(hist.p), var(hist.p))
        
        # generate N samples based on BTF
        boot_out[,,d] = sample.likelihood(n.boot = n.boot,
                                          obs.ccpue = dat$ccpue[dat$date == est_dates[d]],
                                          fit.eos = fit.eos,
                                          rt.method = "beta",
                                          rt.mu = rt.mu, 
                                          sub.harv = dat$ds_charv[dat$date == est_dates[d]], 
                                          sub.cv = dat$ds_charv_cv[dat$date == est_dates[d]],
                                          com.harv = 0
        )[,c("p", "eos", "N", "sub.harv")]
        
        if (input$MEst_do.mcmc) {
          
          # obtain pdf
          like_fun = approxfun(density(boot_out[,"N",d], from = 0, to = dens.max))
          
          # withProgress(message = "Running MCMC:", value = 0, min = 0, max = input$MEst_nc, {
          # run the mcmc
          mcmc.out = lapply(1:input$MEst_nc, FUN = function(i) {
            
            # incProgress(i-1, detail = paste("Chain #", i))
            
            MH(init = runif(1, 1e5, 3e5),
               ni = input$MEst_ni, nb = input$MEst_nb, nt = 1, prop.sig = input$MEst_prop.sig, 
               like.fun = like_fun, prior.fun = prior_fun)
          }) # lapply
          
          post_samp = matrix(unlist(lapply(mcmc.out, function(x) x$post.samp)),  input$MEst_ni - input$MEst_nb, input$MEst_nc)
          
          # MCMC diagnostics
          ess[[d]] = unname(effectiveSize(as.mcmc.list(lapply(as.data.frame(post_samp), mcmc))))
          accept_rate[[d]] = unlist(lapply(mcmc.out, function(x) x$accept.rate))
          if (input$MEst_nc > 1) {
            bgr[[d]] = round(unname(gelman.diag(x = as.mcmc.list(lapply(as.data.frame(post_samp), mcmc)), multivariate = F)[[1]][1,1]),2)
          } else {
            bgr[[d]] = NA
          }
          # }) # progress over chains
          
          post_summ = summ(as.numeric(post_samp))
        } # if do.mcmc
        
        like_summ = summ(boot_out[,"N",d])
        p_summ = summ(boot_out[,"p",d])
        eos_summ = summ(boot_out[,"eos",d])
        harv_summ = summ(boot_out[,"sub.harv",d])
        
        ests[[d]] = data.frame(date = est_dates[d], 
                               source = c(prior_name, like_name, post_name, p_name, eos_name, harv_name),
                               stat = rep(names(like_summ), 6),
                               value = c(prior_summ, like_summ, post_summ, p_summ, eos_summ, harv_summ))
        
        
        
      } # for()
    }) # progress over days
    
    out$MEstSaved = list(
      est_dates = est_dates,
      bgr = unlist(bgr),
      accept_rate = unname(unlist(lapply(accept_rate, mean))),
      ess = unname(unlist(ess)),
      ests = ests
    )
  }) # observeEvent
  
  output$DL_MEst_zip = downloadHandler(
    filename = function() {
      name = create_export_filename(date = c(input$MEst_fdate, input$MEst_ldate), suffix = input$zip_suffix, ext = ".zip")
      return(name)
    },
    content = function(file) {
      files = NULL
      for (d in 1:length(out$MEstSaved$est_dates)) {
        fileName = create_export_filename(date = out$MEstSaved$est_dates[d], suffix = NULL, ext = input$MEst_ext)
        if (input$MEst_ext == ".csv") {
          write.csv(out$MEstSaved$ests[[d]], fileName, row.names = F)
        } else {
          write.table(out$MEstSaved$ests[[d]], fileName, row.names = F)
        }
        files = c(fileName, files)
      }
      zip(file, files)
    }
  )
  
  ### EVERYTHING FOR POLICY TAB1 ###
  output$PlcyInfoChoice1 = renderUI({
    if (input$do.mcmc) {
      checkboxGroupInput(inputId = "knowledge", label = NULL, choices = list("Forecast" = 1, "BTF Only" = 2, "Updated" = 3), selected = 1, inline = T)
    } else {
      checkboxGroupInput(inputId = "knowledge", label = NULL, choices = list("Forecast" = 1, "BTF Only" = 2), selected = 1, inline = T)
    }
  })
  
  # toggle the Hcand_lim slider when link is clicked
  onclick("toggle_Hcand_lim", toggle("Hcand_lim_id", anim = T))
  
  # other toggles for this tab
  observe({
    toggle(id = "harv_cv_id", condition = input$harv != 0 & !is.na(input$harv))  # toggle the harvest cv input
    toggleState("Plcy1Update", condition = !is.null(out$EstOut))                 # toggle the update button
    toggleState("Plcy1Reset", condition = !is.null(out$PolicyOut))               # toggle the reset button
    toggleState("DL_Plcy_Plot", condition = !is.null(out$PolicyOut))             # toggle the download plot button
    toggleState("DL_Plcy_Table", condition = !is.null(out$PolicyOut))            # toggle the download table button
    toggle(id = "RA1", condition = !is.null(out$EstOut))
    toggle(id = "RA1_message", condition = is.null(out$EstOut))
  })
  
  output$p_star_help = renderUI({
    if (input$risk_direction == 1) {
      helpText("Maximum acceptable Pr(S < Target)")
    } else {
      helpText("Minimum acceptable Pr(S > Target)")
    }
  })
  
  output$p_star_slide = renderUI({
    if (input$risk_direction == 1) {
      sliderInput(inputId = "p.star", label = NULL, ticks = F,
                  min = 0, max = 1, value = p.star_def)
      
    } else {
      sliderInput(inputId = "p.star", label = NULL, ticks = F,
                  min = 0, max = 1, value = 1 - p.star_def)
    }
  })
  
  observeEvent(input$Plcy1Update, {
    # bin dimensions
    bin_dims = get_bin_dim(0, 1e6, 1001)
    N = bin_dims$mp
    Hcand = seq(input$Hcand_lim[1], input$Hcand_lim[2], 1000)
    
    # determine which direction risk is framed
    direction = ifelse(input$risk_direction == 1, "less", "greater")
    
    # produce random current harvests taken
    harv_samp = gen_harv_samps(n.boot, mean = input$harv, cv = input$harv.cv)
    
    # critical escapement levels for table
    Scrit = c(65000, 120000, 95000, 110000); nms = paste("less", Scrit/1000, sep = "")
    if (1 %in% input$knowledge) {
      prior_samp = out$EstOut$prior_samp
      prior_cdf = ecdf(prior_samp - harv_samp)
      risk_p_prior = find.risk.p(Sobj = input$Sobj, Hcand = Hcand, cdf = prior_cdf, direction = direction)
      Htarget_prior = Hcand[closest.index(risk_p_prior, input$p.star)]
      ES_prior = mean(prior_samp - harv_samp - Htarget_prior); ES_prior = ifelse(ES_prior < 0, 0, ES_prior)
      Pcrit_prior = find.risk.p(Sobj = Scrit, Hcand = Htarget_prior, cdf = prior_cdf, direction = "less"); names(Pcrit_prior) = nms
    } else {
      risk_p_prior = NULL
      Htarget_prior = NULL
      ES_prior = NULL
      Pcrit_prior = NULL
    }
    
    if (2 %in% input$knowledge) {
      like_samp = out$EstOut$boot_out[,"N"]
      like_cdf = ecdf(like_samp - harv_samp)
      risk_p_like = find.risk.p(Sobj = input$Sobj, Hcand = Hcand, cdf = like_cdf, direction = direction)
      Htarget_like = Hcand[closest.index(risk_p_like, input$p.star)]
      ES_like = mean(like_samp - harv_samp - Htarget_like); ES_like = ifelse(ES_like < 0, 0, ES_like)
      Pcrit_like = find.risk.p(Sobj = Scrit, Hcand = Htarget_like, cdf = like_cdf, direction = "less"); names(Pcrit_like) = nms
    } else {
      risk_p_like = NULL
      Htarget_like = NULL
      ES_like = NULL
      Pcrit_like = NULL
    }
    
    if (3 %in% input$knowledge) {
      post_samp = as.numeric(out$EstOut$post_samp)
      post_samp = post_samp - sample(x = harv_samp, size = length(post_samp), replace = T)
      post_cdf = ecdf(post_samp)
      risk_p_post = find.risk.p(Sobj = input$Sobj, Hcand = Hcand, cdf = post_cdf, direction = direction)
      Htarget_post = Hcand[closest.index(risk_p_post, input$p.star)]
      ES_post = mean(post_samp - Htarget_post); ES_post = ifelse(ES_post < 0, 0, ES_post)
      Pcrit_post = find.risk.p(Sobj = Scrit, Hcand = Htarget_post, cdf = post_cdf, direction = "less"); names(Pcrit_post) = nms
    } else {
      risk_p_post = NULL
      Htarget_post = NULL
      ES_post = NULL
      Pcrit_post = NULL
    }
    
    if (input$risk_direction == 1) {
      ylab = paste("Probability of Escapement Less than", prettyNum(input$Sobj, big.mark = ",", scientific = F))
    } else {
      ylab = paste("Probability of Escapement Greater than", prettyNum(input$Sobj, big.mark = ",", scientific = F))
    }
    
    out$PolicyOut = list(
      Hcand = Hcand,
      N = N,
      risk_p_prior = risk_p_prior,
      risk_p_like = risk_p_like,
      risk_p_post = risk_p_post,
      Htarget_prior = Htarget_prior,
      Htarget_like = Htarget_like,
      Htarget_post = Htarget_post,
      ylab = ylab,
      p.star = input$p.star,
      know.ind = as.numeric(input$knowledge),
      ES_prior = ES_prior,
      ES_like = ES_like,
      ES_post = ES_post,
      Pcrit_prior = Pcrit_prior,
      Pcrit_like = Pcrit_like,
      Pcrit_post = Pcrit_post,
      direction = direction
    )
  })
  
  observeEvent(input$Plcy1Reset, {
    out$PolicyOut = NULL
  })
  
  output$PolicyPlot = renderPlot({
    create_plcy_plot(out, input)
  })
  
  output$DL_Plcy_Plot = downloadHandler(
    filename = function() {
      paste("PlcyPlot_", paste(unlist(str_split(input$cdate, pattern = "/")), collapse = "_"), ".png", sep = "")
    },
    content = function(file) {
      ppi = 600
      png(file, h = 5 * ppi, w = 7 * ppi, res = ppi)
      create_plcy_plot(out, input)
      dev.off()
    }
  )
  
  output$PolicyTable = renderTable(striped = T, bordered = T, hover = T, {
    if (!is.null(out$PolicyOut$know.ind)) {
      create_plcy_table(out)
    } else {
      NULL
    }
  })
  
  output$DL_Plcy_Table = downloadHandler(
    filename = function() {
      paste(paste("PlcyTable", input$p.star, input$Sobj, paste(unlist(str_split(input$cdate, pattern = "/")), collapse = "_"), sep = "_"), ".png", sep = "")
    },
    content = function(file) {
      ppi = 600
      png(file, h = 2.5 * ppi, w = 5 * ppi, res = ppi)
      grid.table(create_plcy_table(out), rows = NULL)
      dev.off()
    }
  )
  
  ### EVERYTHING FOR POLICY TAB2 ###
  # toggles for this tab
  observe({
    toggle(id = "harv2_cv_id", condition = input$harv2 !=0 & !is.na(input$harv2))  # toggle the harvest cv input
    toggleState("Plcy2Update", condition = !is.null(out$EstOut))                   # toggle the update button
    toggleState("Plcy2Reset", condition = !is.null(out$PlcyOut2))                  # toggle the clear button
    toggleState("DL_PlcyTable2", condition = !is.null(out$PlcyOut2))               # toggle the download button
    toggle(id = "RA2", condition = !is.null(out$EstOut))
    toggle(id = "RA2_message", condition = is.null(out$EstOut))
  })
  
  output$PlcyInfoChoice2 = renderUI({
    if (input$do.mcmc) {
      radioButtons(inputId = "knowledge2", label = NULL, choices = list("Forecast" = 1, "BTF Only" = 2, "Updated" = 3), selected = 1, inline = T)
    } else {
      radioButtons(inputId = "knowledge2", label = NULL, choices = list("Forecast" = 1, "BTF Only" = 2), selected = 1, inline = T)
    }
  })
  
  observeEvent(input$Plcy2Update, {
    
    if (input$knowledge2 == "1") {
      samps = out$EstOut$prior_samp
      n.h.boot = n.boot
    }
    if (input$knowledge2 == "2") {
      samps = out$EstOut$boot_out[,"N"]
      n.h.boot = n.boot
    }
    if (input$knowledge2 == "3") {
      samps = as.numeric(out$EstOut$post_samp)
      n.h.boot = length(samps)
    }
    
    # samples of harvest already taken
    harv_samps = gen_harv_samps(n.h.boot, input$harv2, input$harv2.cv)
    
    # critical escapement levels for table
    Scrit = c(65000, 120000, 95000, 110000); nms = paste("less", Scrit/1000, sep = "")
    
    # cdf of escapement if no more fish are harvested
    cdf = ecdf(samps - harv_samps)
    
    # calculate expected escapement under the three candidates
    ES1 = mean(samps - harv_samps - input$Hobj1); ES1 = ifelse(ES1 < 0, 0, ES1)
    ES2 = mean(samps - harv_samps - input$Hobj2); ES2 = ifelse(ES2 < 0, 0, ES2)
    ES3 = mean(samps - harv_samps - input$Hobj3); ES3 = ifelse(ES3 < 0, 0, ES3)
    
    Pcrit_1 = round(find.risk.p(Sobj = Scrit, Hcand = input$Hobj1, cdf = cdf, direction = "less"), 2); names(Pcrit_1) = nms
    Pcrit_2 = round(find.risk.p(Sobj = Scrit, Hcand = input$Hobj2, cdf = cdf, direction = "less"), 2); names(Pcrit_2) = nms
    Pcrit_3 = round(find.risk.p(Sobj = Scrit, Hcand = input$Hobj3, cdf = cdf, direction = "less"), 2); names(Pcrit_3) = nms
    
    p_in_EG1 = unname(Pcrit_1["less120"] - Pcrit_1["less65"])
    p_in_EG2 = unname(Pcrit_2["less120"] - Pcrit_2["less65"])
    p_in_EG3 = unname(Pcrit_3["less120"] - Pcrit_3["less65"])
    
    p_above120_1 = 1 - Pcrit_1["less120"]
    p_above120_2 = 1 - Pcrit_2["less120"]
    p_above120_3 = 1 - Pcrit_3["less120"]
    
    p1 = unname(c(Pcrit_1["less65"], p_in_EG1, p_above120_1, Pcrit_1["less95"], Pcrit_1["less110"]))
    p2 = unname(c(Pcrit_2["less65"], p_in_EG2, p_above120_2, Pcrit_2["less95"], Pcrit_2["less110"]))
    p3 = unname(c(Pcrit_3["less65"], p_in_EG3, p_above120_3, Pcrit_3["less95"], Pcrit_3["less110"]))

    out$PlcyOut2 = list(
      ES1 = ES1,
      ES2 = ES2,
      ES3 = ES3,
      p1 = p1,
      p2 = p2,
      p3 = p3,
      Hobj1 = input$Hobj1,
      Hobj2 = input$Hobj2,
      Hobj3 = input$Hobj3
    )
  })
  
  observeEvent(input$Plcy2Reset, {
    out$PlcyOut2 = NULL
  })
  
  output$PlcyTable2 = renderTable(striped = T, bordered = T, hover = T, {
    if (!is.null(out$PlcyOut2$p1)) {
      create_plcy2_table(out)
    } else {
      NULL
    }
  })
  
  output$DL_PlcyTable2 = downloadHandler(
    filename = function() {
      src = ifelse(input$knowledge2 == "1", "fcst", ifelse(input$knowledge2 == "2", "btf", "update"))
      d = paste(unlist(str_split(input$cdate, pattern = "/")), collapse = "_")
      paste(paste("Compare_Hobj", src, d, sep = "_"), ".png", sep = "")
    },
    content = function(file) {
      ppi = 600
      png(file, h = 2 * ppi, w = 5.5 * ppi, res = ppi)
      grid.table(create_plcy2_table(out), rows = NULL)
      dev.off()
    }
  )
  
  ### EVERYTHING FOR TRACK LEARNING ###
  observe({
    toggleState("read_learn", condition = {
      if (!is.null(input$learn_files)) {
        ifelse(nrow(input$learn_files) < 2, F, T)
      } else {
        F
      }
    })
    toggleState("learn_start_over", condition = !is.null(out$LearnOut))
    toggle("learn_opts", condition = !is.null(out$LearnOut), anim = T)
    toggleState("learn_clear", condition = !is.null(out$LearnOut$plot_dat))
    toggleState("DL_learn_plot", condition = !is.null(out$LearnOut$plot_dat))
  })
  
  output$learn_read_text = renderUI({
    if (is.null(input$learn_files)) {
      # helpText("No files selected")
      NULL
    } else {
      if (nrow(input$learn_files) < 2) {
        helpText("Please select 2 or more files")
      } else {
        NULL
      }
    }
  })
  
  observeEvent(input$read_learn, {
    learn_dates = NULL
    for (i in 1:nrow(input$learn_files)) {
      tmp_file = as.character(input$learn_files$datapath[i])
      tmp_ext = substr(tmp_file, nchar(tmp_file) - 3, nchar(tmp_file))
      if (tmp_ext == ".csv") {
        tmp = read.csv(tmp_file, stringsAsFactors = F, na.strings = "NA")
      } else {
        tmp = read.table(tmp_file, stringsAsFactors = F, header = T, na.strings = "NA")
      }
      tmp_date = unique(tmp$date)
      learn_dates = c(learn_dates, tmp_date)
      if (i == 1) dat = tmp else dat = rbind(dat, tmp)
    }
    
    dat = merge(dat, dates, by = "date")
    fdate = unique(dat$date[dat$day == min(dat$day)])
    ldate = unique(dat$date[dat$day == max(dat$day)])
    
    updateTextInput(session = session, inputId = "learn_fdate", label = "First Date:", value = fdate)
    updateTextInput(session = session, inputId = "learn_ldate", label = "Last Date:", value = ldate)
    
    out$LearnOut = list(
      learn_dat = dat,
      fdate = fdate,
      ldate = ldate,
      learn_dates = learn_dates
    )
  })
  
  output$learn_read_text2 = renderUI({
    if (!is.null(out$LearnOut)) {
      if (any(duplicated(out$LearnOut$learn_dates))) {
        helpText(strong("WARNING:"),  "Two or more files from the same date were uploaded. It is advisable to re-select the files so only one per date is used.")
      } else {
        NULL
      }
    } else {
      NULL
    }
  })
  
  observeEvent(input$learn_start_over, {
    out$LearnOut = NULL
    reset(id = "learn_files")
  })
  
  output$LearnSlide = renderUI({
    if (input$learn_est %in% c("prior", "likelihood", "posterior")) {
      sliderInput(inputId = "learn_ylim", label = NULL,
                  min = 0, max = 800000, step = 10000, value = c(0, 250000), ticks = F)
    } else {
      if (input$learn_est == "eos") {
        sliderInput(inputId = "learn_ylim", label = NULL,
                    min = 0, max = 2000, step = 50, value = c(0, 700), ticks = F)
      } else {
        sliderInput(inputId = "learn_ylim", label = NULL,
                    min = 0, max = 1, step = 0.01, value = c(0, 1), ticks = F)
      }
    }
  })
  
  observeEvent(input$learn_update, {
    dat = out$LearnOut$learn_dat
    dat = dat[dat$day >= dates$day[dates$date == input$learn_fdate],]
    dat = dat[dat$day <= dates$day[dates$date == input$learn_ldate],]
    dat = dat[dat$source == input$learn_est,]
    
    out$LearnOut$plot_dat = list(
      dat = dat,
      ylim = input$learn_ylim,
      xlim = c(dates$day[dates$date == input$learn_fdate], dates$day[dates$date == input$learn_ldate]),
      ylab = ifelse(input$learn_est %in% c("prior", "likelihood", "posterior"), "Run Size (1000s)",
                    ifelse(input$learn_est == "eos", "EOS BTF CCPUE", "Proportion of Run Complete")),
      legend_loc = input$legend_loc
    )
  })
  
  observeEvent(input$learn_clear, {
    out$LearnOut$plot_dat = NULL
  })
  
  output$LearnPlot = renderPlot({
    if (!is.null(out$LearnOut$plot_dat)) {
      create_learn_plot(out)
    } else {
      NULL
    }
  })
  
  output$DL_learn_plot = downloadHandler(
    file = function() {
      create_export_filename(date = c(input$learn_fdate, input$learn_ldate),
                             suffix = NULL, ext = ".png", base = paste(input$learn_est, "_", sep = ""))
    },
    content = function(file) {
      ppi = 600
      png(file, h = 5 * ppi, w = 7 * ppi, res = ppi)
      create_learn_plot(out)
      dev.off()
    }
  )
  
  # download handler for the historical data tab
  output$DL_HistData = downloadHandler(
    filename = function() {
      "BayesTool_HistData.html"
    },
    
    content = function(file) {
      file.copy(paste(data_dir, "BayesTool_HistData.html", sep = "/"), file)
    }
  )
  
}

##### SHINY() #####
shinyApp(ui = ui, server = server)

