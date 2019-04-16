##### APP FOR HOSTING THE USER MANUAL #####

library(shiny)
library(shinythemes)

ui = fluidPage(
  
  wellPanel(
    width = 6,
    h4(strong("Please Help Us Out")),
    
    p("If you plan to view this webpage for an extended period of time,
        or if you plan to revisit it frequently, we request that you download the documentation to your computer.
        The downloaded HTML file can be opened in your web-brower, will display just as below, 
        and can be viewed without an internet connection."),
    p("We have a limited amount of time per month that the server that hosts the Bayes' Tool and the associated documentation can be running,
        and downloading the documentation as opposed to viewing it on this webpage will help
        ensure this time is used for actual use of the tool. Please take your time exploring the tool itself."),
    
    p(strong("Thank You!"), style = "margin: 0;"),
    
    p(em("-The Bayes' Tool Developers"), style = "margin: 0;"),
    br(),
    downloadLink(outputId = "DL", label = p(icon("download"), "Download this Documentation")),
    p("After downloading, please close this tab in your browser.")
  ),
  
  # includeHTML(path = "BayesTool_TechDoc.html")
  includeHTML(path = "Dummy.html")
  
)

server = function(input, output) {
  output$DL = downloadHandler(
    filename = function() {
      # "BayesTool_TechDoc.html"
      "Dummy.html"
    },
    content = function(file) {
      # file.copy("BayesTool_TechDoc.html", file)
      file.copy("Dummy.html", file)
    }
  )
}

shinyApp(ui = ui, server = server)
