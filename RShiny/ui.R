#**************************************************************
#-                     TLD SHINY UI
#-                Oelkuct, M., Reed, J. 
#-
#**************************************************************

#**************************************************************
#- Libraries
#**************************************************************

library(shiny)
library(dplyr)
library(shinydashboard)
library(ggplot2)
library(plyr)
library(plotly)
library(optparse)
library(stats)
library(parallel)
library(foreach)
library(doParallel)
library(shinyjs)
for (pkg in c('shiny', 'dplyr', 'shinydashboard', 'plyr', 'plotly', 'optparse',
              'stats', 'parallel', 'foreach', 'doParallel', 'shinyjs')){
  suppressPackageStartupMessages()
}


#**************************************************************
#- Functions
#**************************************************************

#- Minimal file input button
fileInputNoExtra<-function(inputId, label, multiple = FALSE, accept = NULL, width = NULL, buttonLabel = "Browse...", placeholder = "No file selected"){
  
  restoredValue <- restoreInput(id = inputId, default = NULL)
  if (!is.null(restoredValue) && !is.data.frame(restoredValue)) {
    warning("Restored value for ", inputId, " has incorrect format.")
    restoredValue <- NULL
  }
  if (!is.null(restoredValue)) {
    restoredValue <- toJSON(restoredValue, strict_atomic = FALSE)
  }
  inputTag <- tags$input(id = inputId, name = inputId, type = "file", 
                         style = "display: none;", 
                         `data-restore` = restoredValue)
  if (multiple) 
    inputTag$attribs$multiple <- "multiple"
  if (length(accept) > 0)
    inputTag$attribs$accept <- paste(accept, collapse = ",")
  
  tags$label(class = "input-group-btn", type="button", style=if (!is.null(width)) paste0("width: ", validateCssUnit(width),";","padding-right: 5px; padding-bottom: 5px; display: inline-block;"),
             span(class = "btn btn-default btn-file",type="button", buttonLabel, inputTag, style=if (!is.null(width)) paste0("width: ", validateCssUnit(width),";","border-radius: 5px; padding-bottom: 5px;"))
  )
}


#**************************************************************
#- Layout
#**************************************************************

body <- dashboardBody(
  fluidRow(
    tabBox(
      id = "tabset1", height = "2400px",  width = 12, selected = 'Upload Data File' ,
      
      tabPanel('Length Determination and Significance Overview', '', width = 12, height = 12,
               tags$head(tags$style(
                 type="text/css",
                 "#image img {max-width: 100%; width: 180%; margin-left:-30px;%;margin-top:2%;}"
               )),
               tags$head(tags$style(
                 type="text/css",
                 "#image2 img {max-width: 300%; width:100%;margin-top:45%;margin-left:0px;}"
               )),
               tags$head(tags$style(
                 type="text/css",
                 "#image3 img {max-width: 100%; width: 220%; margin-left:0px;margin-top:80%;}"
               )),
               
               
               
               tags$head(tags$style(
                 type="text/css",
                 "#image4 img {max-width: 100%; width:100%;margin-top:-220%;margin-left:105%;}"
               )),
               
               
               
               
               
               box(
                 imageOutput("image"), 
                 
                 imageOutput("image2"),
                 imageOutput("image3"),
                 #
                 
                 tags$head(tags$style(
                   type="text/css",
                   
                   "
            #image8{
            max-width: 100vh;
            margin: auto;
            width:100%;margin-top:-15%;margin-left:99%;position:relative;
        }
            #image7{
            max-width: 100vh;
            margin: auto;
            width:100%;margin-top:-4%;margin-left:99%;position:relative;
        }
            #image6 {
            max-width: 50vh;
            margin: auto;
            width:100%;margin-top:-65%;margin-left:99%;position:relative;
        }
            #image5 {
            max-width: 100vh;
            margin: auto;
            width:100%;margin-top:-115%;margin-left:95%;position:relative;
        }
            #image4 {
            max-width: 95vh;
            max-height: 125vh;
            margin: auto;
            width:100%;margin-top:7.7%;margin-left:-8.5%;position:relative;
        }")),
                 
                 (div(style='display:grid;height:100vh;margin:0;padding:0;',
                      
                      imageOutput('image4'),
                      
                      
                      
                      imageOutput('image5'),
                      imageOutput('image6'),
                      imageOutput('image7'),
                      imageOutput('image8')
                 ))
               )
               
               
      ),
      
      tabPanel("References", "", width = 12, height = 12,
               
               tags$head(tags$style(
                 type="text/css",
                 "#image6 img {max-width: 130%; width: 250%; margin-left:-5px;margin-top:5%;}"
               )),
               # box(imageOutput("image8"))
               
               
      ),
      
      tabPanel('Upload Data File', '', width = 12, height = 12,
               (div(style='height:720px;overflow-y:scroll;',
                    DT::dataTableOutput("table")               )
               ),
               (div(style='height:100px;position:absolute;overflow-y:hidden;left:700px;top:50px;',
                    
                    fileInputNoExtra("FileInput",label="",accept=".csv",buttonLabel=list(icon("folder"),"Upload Results Data FIle"),width="300px"))
                
                
               )
      ),
      tabPanel('Generate Plots', '',
               tabBox( width = 12,
                       selected='Number of Reads by End',
                       
                       
                       tabPanel('Number of Reads by End', '',
                                (div(style='margin-left:5%;margin-top:5%;width:1600px;margin-top:2%;',
                                     fluidRow(
                                        plotlyOutput('myPlot', height = 800, width = 1800)
                                     )
                                     
                                )
                                )
                       ),
                       tabPanel('Telomere Length Distribution', '', width = 12, height = 12,
                                
                                
                                plotOutput('distplot', height = 870, width = 1820)
                                
                                
                                
                                
                       ),
                       tabPanel('Length by Chromosome End (Violin Plot)', '',
                                (div(style='height:768px;overflow-y:scroll;margin-left:5%;',
                                     plotlyOutput('violinPlot', height = 530, width = 740),
                                     br(),
                                     br(),
                                     plotlyOutput('violinPlot2', height = 530, width = 740)
                                     
                                )
                                )
                                
                                
                       )
                       
               )
      ),
      
      
      
    )
  )
  
)


sidebar <- dashboardSidebar(
  width = "0px"
)

ui = dashboardPage(
  dashboardHeader(title = "TeCAT"),
  sidebar,
  body
)