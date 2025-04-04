#**************************************************************
#-                     TLD SHINY Server
#-                Oelkuct, M., Reed, J. 
#-
#**************************************************************

#**************************************************************
#- Libraries + SRC
#**************************************************************

source("./scripts/def_plotting.R", local = F)
source("./scripts/modes.R", local = F)
source("./scripts/shiny_plots.R", local = F)

for (pkg in c('shiny', 'dplyr', 'shinydashboard', 'ggplot2', 'plyr', 'plotly', 'optparse',
              'stats', 'parallel', 'foreach', 'doParallel', 'shinyjs')){
  suppressPackageStartupMessages(pkg)
}

options(shiny.autoreload = TRUE)
registerDoParallel(cores = 4) 

server <- function(input, output, session) {
  
#-- custom css mod, remove sidebar-toggle
removeUI("a.sidebar-toggle", multiple = T)
addClass(selector = "body", class = "sidebar-collapse")
  
#**************************************************************
#- DataSet & Plot Functions
#**************************************************************

datasetInput <- reactive({ 
  if (is.null(input$FileInput)){
    result.df <- read.csv("./data/example_data/output_results.csv")
    result.df
  } else {
    infile <- input$FileInput
    read.csv(infile$datapath, header = TRUE)
  }
  })
  
  output$table <- DT::renderDataTable(datasetInput())
  
  output$myPlot <- renderPlotly({
    plot_num_reads(datasetInput())
  })
  
  # Telomere Read Count Plot
  output$myPlot2 <- renderPlotly({
    barpldf <- datasetInput() %>%
      group_by(seqnames, s.name, chr.end) %>%
      summarise(`Telomere Read Count` = n())
    
    fig2 <- filter(barpldf, chr.end == "5'") %>%
      plot_ly(width = 900) %>%
      layout(yaxis = list(title = list(text = "Telomere Read Count", font = list(size = 15))), xaxis = list(title = list(text = "Chromosome", font = list(size = 15))))
    fig2 <- fig2 %>%
      add_trace(x = filter(barpldf, s.name == "s288c")$seqnames, y = filter(barpldf, s.name == "s288c")$`Telomere Read Count`, type = "bar", name = "s288c") %>%
      layout(yaxis = list(title = list(text = "Telomere Read Count", font = list(size = 15))), xaxis = list(title = list(text = "Chromosome", font = list(size = 15))), title = list(text = "5' End Number of Reads by Sample", font = list(size = 20), y = 1))
    foo <- filter(barpldf, !s.name == "s288c")

    for (item in input$show_vars) {
      if (item != "s288c") {
        fig2 <- fig2 %>%
          add_trace(x = filter(foo, s.name == item)$seqnames, y = filter(foo, s.name == item)$`Telomere Read Count`, type = "bar", name = item) %>%
          layout(yaxis = list(title = list(text = "Telomere Read Count", font = list(size = 17))), xaxis = list(title = list(text = "Chromosome", font = list(size = 17))))
      }
    }
    fig2
  })

  # Distribution Plot
  output$distplot <- renderPlot({
    get_dens_plot(datasetInput())
  })
 # Bar Plot
  output$plot2 <- renderPlot({
    bar_pl
  })
  #----- Violin Plot Select User Input -----
  output$violinPlot <- renderPlotly({
    viola_plot(datasetInput(), "5'")
  })


  output$violinPlot2 <- renderPlotly({
    viola_plot(datasetInput(), "3'")

  })

  #-endbam DataTable
  output$endbam <- DT::renderDT({
    DT::datatable(head(select(result.df, -1)))
  })

  #**************************************************************
  #- Image Output
  #**************************************************************
  
  output$png <- renderUI({
    tags$a(img(src = "www/paper.png", width = "200px", height = "200px"), href = "https://www.ebi.ac.uk/gxa/experiments/E-MTAB-2770/Results")
  })
  output$paperpng <- renderImage(
    {
      filename <- "www/githubclean.png"
      list(
        src = filename
      )
    },
    deleteFile = FALSE
  )
  images <- c("www/tel1.png", "www/tel2.png", "www/tel3.png", "www/tel4.png", "www/tel5.png", "www/tel6.png", "www/tel7.png", "www/vign6.png")
  output$image <- renderImage(
    {
      filename <- normalizePath(file.path(paste0(images[1])))
      list(
        src = filename
      )
    },
    deleteFile = FALSE
  )
  output$image2 <- renderImage(
    {
      filename <- normalizePath(file.path(paste0(images[2])))
      list(
        src = filename
      )
    },
    deleteFile = FALSE
  )
  output$image3 <- renderImage(
    {
      filename <- normalizePath(file.path(paste0(images[3])))
      list(
        src = filename
      )
    },
    deleteFile = FALSE
  )
  output$image4 <- renderImage(
    {
      filename <- normalizePath(file.path(paste0(images[4])))
      list(
        src = filename
      )
    },
    deleteFile = FALSE
  )
  output$image5 <- renderImage(
    {
      filename <- normalizePath(file.path(paste0(images[5])))
      list(
        src = filename
      )
    },
    deleteFile = FALSE
  )
  output$image6 <- renderImage(
    {
      filename <- normalizePath(file.path(paste0(images[6])))
      list(
        src = filename
      )
    },
    deleteFile = FALSE
  )
  output$image7 <- renderImage(
    {
      filename <- normalizePath(file.path(paste0(images[7])))
      list(
        src = filename
      )
    },
    deleteFile = FALSE
  )
  output$image8 <- renderImage(
    {
      filename <- normalizePath(file.path(paste0(images[8])))
      list(
        src = filename
      )
    },
    deleteFile = FALSE
  )
  
  # Current Selected Tab
  output$tabset1Selected <- renderText({
    input$tabset1
  })
  output$tabset2Selected <- renderText({
    input$tabset2
  })
}
