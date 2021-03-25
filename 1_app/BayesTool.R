## ------------------------------------------------------------- ##
## IN-SEASON CHINOOK RUN SIZE ASSESSMENT AND DECISION TOOL       ##
## ------------------------------------------------------------- ##

## -------------------------- ##
## CODE BY: BEN STATON        ##
## -------------------------- ##
## DATE V1 STARTED: 3/14/2018 ##
## -------------------------- ##
## DATE UPDATED:    4/16/2020 ##
## -------------------------- ##

##### NECESSARY PACKAGES #####
# source("0a_install_packages.R")
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
library(gridExtra)       # for turning tables to images

##### SESSION SETUP #####
# CLEAR THE WORKSPACE
rm(list = ls(all = T))

# LOAD ALL NECESSARY FUNCTIONS 
source("0b_util_funcs.R")
source("0c_update_funcs.R")
source("0d_output_funcs.R")

data_dir = "Inputs"

# PREPARE DATA FILES
# a date lookup key
dates = dates.prep(year = 2020)

# historical BTF data
btf_data = hist.index.prep(
  f.date = "6/1", l.date = "8/24",
  chinook = read.csv(file.path(data_dir, "Chinook Daily BTF.csv"), stringsAsFactors = F)
)

# historical total run abundance data
N_data = read.csv(file.path(data_dir, "run_size_ests.csv"), stringsAsFactors = F)

##### SET DEFAULT INPUT VALUES #####
# forecast mean and cv defaults
fcst_mean_def = round(N_data$N[nrow(N_data)], -3)
true_runs = N_data[2:nrow(N_data),"N"]; fcst_runs = N_data[1:(nrow(N_data)-1),"N"]
fcst_cv_def = round(sig2cv(sd(log(fcst_runs) - log(true_runs))),2)

# btf ccpue and date defaults
cdate_def = "6/1"
ccpue_def = 0

# mcmc defaults
ni_def = 76000     # total samples (including thin and burnin)
nb_def = 1000      # number of burnin
nc_def = 4         # number of chains
prop_sig_def = 0.7 # sd of the lognormal proposal dist.

# total harvest mean and cv defaults
harv_mean_def = 0
harv_cv_def = 0.1

# risk analysis defaults
H_cand_lim_def = c(0, 100000)
Sobj_def = 92500
p.star_def = 0.5

# number to sample from likelihood and prior
n_mc = 1e6

# maximum run size to calculate density for
dens.max = 5e6

# calculate times
sec_prior_like = 1.05
sec_per_mcmc = 1.919999e-05

# copy the documentation to the Inputs location to ensure they are the most current versions
# file.copy(from = "../2_docs/BayesTool_TechDoc/BayesTool_TechDoc.pdf", to = "Inputs/BayesTool_TechDoc.pdf", overwrite = T)
# file.copy(from = "../2_docs/BayesTool_UserMan/BayesTool_UserMan.pdf", to = "Inputs/BayesTool_UserMan.pdf", overwrite = T)
# file.copy(from = "../2_docs/BayesTool_HistData/BayesTool_HistData.pdf", to = "Inputs/BayesTool_HistData.pdf", overwrite = T)

##### CREATE USER INTERFACE #####
ui = navbarPage(
  title = strong("Chinook Salmon In-season Bayesian Risk Assessment Tool"), 
  windowTitle = "Bayes' Tool",
  theme = shinytheme("cerulean"),
  footer = "Version 1.4.0 (for use in 2021)",
  
  # ESTIMATION TAB
  tabPanel(
    title = "Estimation", 
    icon = icon("calculator", "fa-1x"),
    useShinyjs(),
    
    sidebarLayout(
      sidebarPanel(
        width = 5,
        
        # run size forecast input
        h4(strong("Run Size Forecast")),
        splitLayout(
          numericInput(inputId = "prior.mean", label = "Mean:", min = 0, max = 4e6, step = 5000, value = fcst_mean_def),
          hidden(
            div(
              id = "prior_cv_id",
              numericInput(inputId = "prior.cv", label = "CV:", min = 0, max = 3, step = 0.01, value = fcst_cv_def)
            )
          )
        ),
        a(id = "toggle_prior_cv", "Show/Hide Uncertainty"),
        
        # BTF information input 
        h4(strong("In-Season Index Information")),
        splitLayout(
          textInput(inputId = "cdate", label = "Today's Date:", value = cdate_def),
          numericInput(inputId = "ccpue", label = "BTF CCPUE:", min = 0, max = 5000, step = 1, value = ccpue_def)
        ),
        
        ### MCMC SETTINGS ###
        # whether to do update and whether to display settings
        h4(strong("Information Update Settings")),
        splitLayout(
          cellWidths = c("35%", "65%"),
          checkboxInput(inputId = "do.mcmc", label = "Perform Update", value = T),
          hidden(
            div(
              id = "toggle_mcmc_def",
              checkboxInput(inputId = "show_mcmc_opts", label = "Display/Change Settings", value = F)
            )
          )
        ),
        
        # number of mcmc samples of each type
        hidden(
          div(
            id = "mcmc_opts",
            splitLayout(
              cellWidths = c("30%", "30%", "20%", "20%"),
              numericInput(inputId = "ni", label = "Samples:", min = 0, max = 1e6, step = 5000, value = ni_def),
              numericInput(inputId = "nb", label = "Burn-in:", min = 0, max = 1e6, step = 5000, value = nb_def),
              numericInput(inputId = "nc", label = "Chains:", min = 1, max = 4, step = 1, value = nc_def),
              numericInput(inputId = "prop.sig", label = "Tune:", min = 0.01, max = 1, step = 0.05, value = prop_sig_def)
            )
          )
        ),
        
        # option to reset mcmc settings if they were changed
        hidden(
          div(
            id = "reset_mcmc_opts_id",
            a(id = "reset_mcmc_opts", p(icon("undo", "fa-1x"), "Restore Default Settings"))
          )
        ),
        
        # insert the time estimate
        uiOutput("time_text"),
        
        # buttons for update and clear results
        actionButton(inputId = "EstUpdate", label = "Calculate", icon = icon("cogs", "fa-1x")),
        actionButton(inputId = "EstReset", label = "Clear", icon = icon("trash", "fa-1x")),
        
        # options for downloading summarized output
        dropdownButton(
          inputId = "export_dropdown", circle = F, status = "primary",
          icon = icon("download"), label = "Export Estimates", 
          
          textInput(inputId = "export_suffix", label = NULL, placeholder = "Optional File Suffix"),
          radioButtons(inputId = "export_ext", label = NULL, 
                       choices = list("CSV" = ".csv", "TXT" = ".txt"), selected = ".csv",
                       inline = T),
          uiOutput("export_FileName"),
          downloadLink(outputId = "DL_Ests", label = "Export File")
        )
      ),
      
      mainPanel(
        width = 7,
        tabsetPanel(
          type = "tabs",
          
          # panel for distribution plot
          tabPanel(
            title = p(icon("area-chart", "fa-1x"), 
                      strong("Distributions"),
                      style = "font-size:16px"),
            # message to perform estimation
            div(
              id = "E1_message",
              br(),
              wellPanel(
                h4(
                  "No estimates generated. Please enter the current 
                  In-Season Index Information and click the",
                  icon("cogs", "fa-1x"), strong("Calculate"), "button.")
              )
              ),
            plotOutput("EstPlot1"),
            br(),
            uiOutput("DL_Dist_Plot_Button")
            ),
          
          # tab for estimates table
          tabPanel(
            title = p(icon("table", "fa-1x"), strong("Table"), style = "font-size:16px"),
            # message to perform estimation
            div(
              id = "E2_message",
              br(),
              wellPanel(
                h4(
                  "No estimates generated. Please enter the current 
                  In-Season Index Information and click the",
                  icon("cogs", "fa-1x"), strong("Calculate"), "button.")
              )
              ),
            br(),
            tableOutput("EstTable"),
            br(),
            uiOutput("DL_EstTable_Button"),
            br()
            ),
          
          # tab for relationship plot
          tabPanel(
            title = p(icon("bar-chart-o", "fa-1x"), strong("Relationship"), style = "font-size:16px"),
            # dropdown menu for plot options
            hidden(
              div(
                id = "EstPlot2_dropdown",
                dropdownButton(
                  circle = T, icon = icon("gear"), status = "primary",
                  sliderInput(inputId = "EstPlot2_xlim", min = 0, max = 1200, value = c(0, 700), step = 10, label = "CCPUE Range", ticks = F, dragRange = T),
                  sliderInput(inputId = "EstPlot2_ylim", min = 0, max = 450000, value = c(0, 250000), step = 5000, label = "Run Size Range", ticks = F, dragRange = T),
                  checkboxInput(inputId = "EstPlot2_predata", label = "Display Pre-2008 Data", value = F)
                )
              )
            ),      
            
            # message to perform the estimation
            div(
              id = "E3_message",
              br(),
              wellPanel(
                h4(
                  "No estimates generated. Please enter the current 
                  In-Season Index Information and click the",
                  icon("cogs", "fa-1x"), strong("Calculate"), "button."
                )
              )
              ),
            plotOutput(outputId = "EstPlot2"),
            uiOutput("DL_EstPlot2_Button")
            ),
          
          # tab for MCMC diagnostics
          tabPanel(
            title = p(icon("thermometer-half", "fa-1x"), strong("Diagnostics"), style = "font-size:16px"),
            
            # message to perform the estimation       
            div(
              id = "E4_message",
              br(),
              wellPanel(
                h4(
                  "No estimates generated. Please enter the current 
                  In-Season Index Information and click the",
                  icon("cogs", "fa-1x"), strong("Calculate"), "button."
                )
              )
              ),
            plotOutput(outputId = "MCMCPlot"),
            br(),
            tableOutput(outputId = "MCMCTable"),
            br(),
            a(id = "toggle_MCMC_help", "Show/Hide More Information"),
            
            # dialog box with info about MCMC diags
            hidden(
              div(
                id = "MCMC_Diag_Help",
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
                ) # well panel
              ) # div
            ) # hidden
            ) # mcmc diagnostics tab
          ) # tabset panel
          ) # main panel
      ) # sidebarpanel
  ), # estimation tab
  
  # RISK ANALYSIS TAB
  navbarMenu(
    title = "Risk Analysis", icon = icon("filter", "fa-1x"),
    useShinyjs(),
    tabPanel(
      title = "Choose Harvest Target",
      sidebarLayout(
        sidebarPanel(
          width = 5,
          div(id = "RA1_message",
              h4("Please complete the", icon("calculator", "fa-1x"), strong("Estimation"), "option before performing risk analysis")
          ),
          hidden(
            div(
              id = "RA1",
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
              
              h4(strong("Escapement Limit")),
              sliderInput(inputId = "Sobj", label = NULL, ticks = F,
                          min = 65000, max = 120000, value = Sobj_def, step = 500),
              
              h4(strong("Risk Tolerance")),
              uiOutput("p_star_help"),
              uiOutput("p_star_slide"),
              
              h5(strong("Frame Risk As:")),
              radioButtons(inputId = "risk_direction", label = NULL,
                           choices = list("Pr(S < Limit)" = 1,
                                          "Pr(S > Limit)" = 2), selected = 1, inline = T),
              
              miniButtonBlock(
                actionButton(inputId = "Plcy1Update", label = "Update", icon = icon("cogs", "fa-1x")),
                actionButton(inputId = "Plcy1Reset", label = "Clear", icon = icon("trash", "fa-1x")),
                downloadButton(outputId = "DL_Plcy_Plot", label = "Plot"),
                downloadButton(outputId = "DL_Plcy_Table", "Table")
              )
            )
          )
        ),
                        
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
    
    tabPanel(
      title = "Compare Harvest Targets",
      sidebarLayout(
        sidebarPanel(
          width = 5,
          div(id = "RA2_message",
              h4("Please complete the", icon("calculator", "fa-1x"), strong("Estimation"), "option before performing risk analysis")
          ),
          hidden(
            div(
              id = "RA2",
              h4(strong("Run Size Estimate")),
              uiOutput("PlcyInfoChoice2"),
              
              h4(strong("Current Total Harvest")),
              
              splitLayout(
                numericInput(inputId = "harv2", label = "Mean:", min = 0, max = 100000, step = 1000, value = harv_mean_def),
                hidden(
                  div(
                    id = "harv2_cv_id",
                    numericInput(inputId = "harv2.cv", label = "CV:", min = 0, max = 1, step = 0.01, value = harv_cv_def)
                  )  
                )
              ),
              
              h4(strong("Candidate Harvest Targets")),
              helpText("In addition to harvest currently taken listed above"),
              splitLayout(
                numericInput("Hobj1", "H1:", value = 0, min = 0, max = 1.5e5, step = 1000),
                numericInput("Hobj2", "H2:", value = 0, min = 0, max = 1.5e5, step = 1000),
                numericInput("Hobj3", "H3:", value = 0, min = 0, max = 1.5e5, step = 1000)
              ),
              actionButton("Plcy2Update", "Update", icon = icon("cogs", "fa-1x")),
              actionButton("Plcy2Reset", "Clear", icon = icon("trash", "fa-1x")),
              downloadButton("DL_PlcyTable2", "Download Table")
            )
          )
        ),
        mainPanel(
          width = 7,
          tableOutput("PlcyTable2")
        )
      )
    )
  ),
  
  ### ABOUT TAB ###
  navbarMenu(
    title = "About", icon = icon("info-circle"),
    tabPanel(
      p(icon("info-circle"), "Overview", style = "margin:0;"),
      fluidRow(
        column(
          width = 8,
          wellPanel(
            h4(strong("Overview")),
            p("This Tool was developed to facilitate the probabilistic treatment of multiple sources of run size information when 
              considering allowable harvest targets for in-season management of Chinook salmon in the Kuskokwim River.
              Treating run size probabilistically has several advantages, for example, it allows managers to consider the 
              probability that escapement will be below some critical threshold if a given number of fish are harvested.",
              strong("As currently constructed, the tool performs these tasks for Chinook salmon in the Kuskokwim River only, 
                     and does not assess risk of failing to meet tributary-level escapement goals.")),
            p("The Tool allows users to easily perform run size predictions by updating the pre-season forecast with new information
              from the Bethel Test Fishery as it accumulates. This update is performed using a statistical framework called", em("Bayesian Inference."),
              "Bayesian Inference provides updated probabilities of different run size outcomes managers are considering after obtaining new information, which can then
              be used to assess the probability of having favorable or unfavorable escapment levels at various levels of harvest."),
            
            p("The statistical framework used in the tool has been evaluated and an", a("article summarizing this work", href = "https://www.nrcresearchpress.com/doi/10.1139/cjfas-2018-0176"), "has been published in a peer-reviewed scientific journal."),
            
            p("For more details on how to use the Tool, please see the", strong(icon("question-circle"), "User Manual.")),
            p("For details on how the Tool works (including a worked example of Bayesian calculations for a 
              simplified problem in this context), see the", strong(icon("file"), "Technical Documentation."))
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
            p(
              "We strongly encourage users of the tool with any questions about how to use this tool
              to consult the user manual. The manual was developed specifically to be a walkthrough
              tutorial for how to use this tool, and we hope users will find it informative."
            ),
            p(
              "If your questions were not answered by the user manual, please see the",
              strong(icon("envelope"), "Contact Page"), "to send us an email.
              After answering your question(s), we will consider updating the user manual to help future users."
            ),
            downloadLink(outputId = "DL_UserMan", label = p(icon("download"), "Download the User Manual (PDF)"))
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
            p("Should users be interested in the details of how this tool works, 
              the motivations for its development, and other technical information,
              we encourage them to download and review the technical documentation."),
            
            p("If you have questions or concerns about the information found in the technical documentation,
              please visit the ", strong(icon("envelope"), "Contact Page"),
              "to send us an email."),
            
            downloadLink(outputId = "DL_TechDoc", label = p(icon("download"), "Download the Technical Documentation (PDF)"))
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
            downloadLink("DL_HistData", label = p(icon("download"), "Download the Historical Data (PDF)"))
            )
          )
        )
        ),
    
    tabPanel(
      p(icon("users"), "Developers", style = "margin:0;"),
      fluidRow(
        column(
          width = 8,
          wellPanel(
            h4(strong("Tool Developers")),
            p("The interface of this tool was developed and is maintained by Dr. Ben Staton. 
              Ben was a Graduate Student Researcher and USFWS Pathways Quantitative Ecologist at the time of
              initial development and continues to update the tool for use by interested parties on a voluntary basis.
              Ben's graduate research focused on the statistical aspects of Chinook salmon assessment
              and management issues in the Kuskokwim River beginning in 2014 and ending in 2019. 
              The statistical framework for the tool was developed jointly by Ben and his advisor, Dr. Matt Catalano,
              who is the Principal Investigator of the", a("Quantitative Fisheries Laboratory", href = "http://sfaas.auburn.edu/catalano/quantitative-fisheries-lab/"),
              "at Auburn University and has studied salmon assessment methods in Alaska since 2009."),
            
            h4(strong("Contributors")),
            p("Many people have given incredibly useful feedback on this tool and on its utility for aiding in management discussions (listed alphabetically):",
              tags$ul(
                tags$li("Dr. Bill Bechtol"),
                tags$li("Dr. Lew Coggins"),
                tags$li("Zach Liller"),
                tags$li("Nick Smith"),
                tags$li("Kuskokwim area fisheries managers")
              )
            ),
            
            p("While the tool and its philosophy are novel to the Kuskokwim River Chinook salmon stock, it builds on pioneering work by Carl Walters, Sandy Buckingham, Stephen Fried, and Ray Hilborn done in the 1970s and 1980s.
              To our knowledge, the tool is the first to allow users to easily perform these types of relatively complex calculations without having to interact with spreadsheets or code."),
            p("The interface of this tool was developed using", a("Program R", href = "https://www.r-project.org/"), "which handles the statistical computing and a package called", a("Shiny", href = "https://shiny.rstudio.com/"),
              "which allows for integration of R code into the web-interface."),
            p("Funding for the development of this tool was provided by the Arctic-Yukon-Kuskokwim Sustainable Salmon Initiative, administered through the Bering Sea Fishermen's Association, through a research grant to Dr. Matt Catalano.")
          )
        )
      )
    ),
    
    tabPanel(
      p(icon("envelope"), "Contact", style = "margin:0;"),
      fluidRow(
        column(
          width = 6,
          wellPanel(
            h4(strong("Contact")),
            p("Should users have any questions about the tool that were not answered by
              the user manual or technical documentation, please feel free to contact the lead tool developer, Ben Staton, 
              at bas0041@auburn.edu"),
            p("We are also ready and willing to hear about any feedback, suggestions for future functionality, or bugs you may find.")
          ) # well
        ) # column
      ) # row
    ), # tabpanel
    
    tabPanel(
      p(icon("code"), "Source Code", style = "margin:0;"),
      fluidRow(
        column(
          width = 6,
          wellPanel(
            h4(strong("Source Code")),
            p("Users of the Tool with experience writing and reading R code may find it beneficial to explore the
              code for the Tool. All files needed to execute the Tool from any users local computer can be found 
              in the", strong(icon("github"), "Github"), a("repository.", href = "https://github.com/bstaton1/kusko-bayes-tool")),
            p("Advanced R users are invited to fork the repository and suggest changes via pull request.")
          ) # well
        ) # column
      ) # row
    ) # tabpanel
    
    
  ) # navbarmenu
) # end of UI


##### CREATE SERVER TASKS #####

server = function(input, output, session) {
  
  ### EVERYTHING FOR SINGLE DAY ESTIMATION TAB ###
  # toggle the prior cv option
  onclick("toggle_prior_cv", toggle(id = "prior_cv_id", anim = F))
  
  # restore mcmc defaults
  onclick(id = "reset_mcmc_opts", reset("mcmc_opts"))
  
  # toggle MCMC help
  onclick(id = "toggle_MCMC_help", toggle(id = "MCMC_Diag_Help", anim = T))
  
  # toggle other options on this tab
  observe({
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
  
  observe({
    toggle(id = "E1_message", condition = is.null(out$EstOut))
    toggle(id = "E2_message", condition = is.null(out$EstOut))
    toggle(id = "E3_message", condition = is.null(out$EstOut))
    toggle(id = "E4_message", condition = is.null(out$EstOut))
  })
  
  out = reactiveValues()
  
  observeEvent(input$EstUpdate, {
    starttime = Sys.time()
    
    # subset historical data
    fit_data = prepare_fit_data(dt = input$cdate, index = btf_data, N = N_data)
    
    # fit the historical regression
    fit = lm(log(N) ~ q * ccpue, data = fit_data)

    # randomization
    withProgress(message = "Drawing random samples:", value = 0, {
      
      # sample the prior
      incProgress(1/2, detail = "Forecast")
      prior.samp = exp(rnorm(n_mc, log(input$prior.mean) - 0.5 * cv2sig(input$prior.cv)^2, cv2sig(input$prior.cv)))
      
      # generate N samples based on index
      incProgress(2/2, detail = "BTF Expansion")
      boot_out = sample_likelihood(fit, pred_ccpue = input$ccpue, pred_q = 2, n_mc = n_mc)

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
    
    if (input$do.mcmc) {
      post_summ = summ(unlist(lapply(mcmc.out, function(x) x$post.samp)))
    } else {
      post_summ = rep(NA, 9)
    }
    
    q = names(prior_summ)
    
    export = data.frame(date = input$cdate,
               source = rep(c("prior", "likelihood", "posterior"), each = length(q)),
               stat = rep(q, 6),
               value = c(prior_summ, like_summ, post_summ))
    
    
    # make the output list
    if (input$do.mcmc) {
      EstOut = list(
        did_mcmc = T,
        boot_out = boot_out,
        prior_samp = prior.samp,
        post_samp = matrix(unlist(lapply(mcmc.out, function(x) x$post.samp)),  input$ni - input$nb, input$nc),
        accept = unlist(lapply(mcmc.out, function(x) x$accept.rate)),
        cdate = input$cdate,
        ccpue = input$ccpue,
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
        ccpue = input$ccpue,
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
  
  output$DL_EstPlot2_Button = renderUI({
    if (!is.null(out$EstOut$boot_out)) {
      downloadButton(outputId = "DL_EstPlot2", label = "Download Plot")
    } else {
      NULL
    }
  })
  
  observe(toggle(id = "EstPlot2_dropdown", condition = !is.null(out$EstOut)))
  output$EstPlot2 = renderPlot({
    create_EstPlot2(input, out, btf_data, N_data)
  })
  
  output$DL_EstPlot2 = downloadHandler(
    filename = function() {
      paste("RelationshipPlot_", paste(unlist(str_split(input$cdate, pattern = "/")), collapse = "_"), ".png", sep = "")
    },
    content = function(file) {
      ppi = 600
      png(file, h = 5 * ppi, w = 7 * ppi, res = ppi)
      create_EstPlot2(input, out, btf_data, N_data)
      dev.off()
    }
  )
  
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
      helpText("Maximum acceptable Pr(S < Limit)")
    } else {
      helpText("Minimum acceptable Pr(S > Limit)")
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
    harv_samp = gen_harv_samps(n_mc, mean = input$harv, cv = input$harv.cv)
    
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
      n.h.boot = n_mc
    }
    if (input$knowledge2 == "2") {
      samps = out$EstOut$boot_out[,"N"]
      n.h.boot = n_mc
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
  
  # download handler for the User Manual
  output$DL_UserMan = downloadHandler(
    filename = function() {
      "BayesTool_UserMan.PDF"
    },
    content = function(file) {
      file.copy(file.path(data_dir, "BayesTool_UserMan.pdf"), file)
    }
  )
  
  # download handler for the User Manual
  output$DL_TechDoc = downloadHandler(
    filename = function() {
      "BayesTool_TechDoc.pdf"
    },
    content = function(file) {
      file.copy(file.path(data_dir, "BayesTool_TechDoc.pdf"), file)
    }
  )
  
  # download handler for the historical data tab
  output$DL_HistData = downloadHandler(
    filename = function() {
      "BayesTool_HistData.pdf"
    },
    content = function(file) {
      file.copy(file.path(data_dir, "BayesTool_HistData.pdf"), file)
    }
  )
}

##### RUN THE APPLICATION #####
shinyApp(ui = ui, server = server)

