library(shiny)
library(shinyAce)
data_path <- "//DATA4/Mass Spec DataRrepository/imgs/"
data_path_rpt <- "//DATA4/Mass Spec DataRrepository/reportfiles/current_report_files/"
SO_options <- list.dirs(data_path, full.names = FALSE )
QCB_options <- list.files(data_path_rpt, pattern = "*.rpt", full.names = FALSE )

shinyUI(fluidPage(
        titlePanel("QC data lookup", "DataLookup"),
        
        sidebarLayout(
                sidebarPanel(
                        selectInput(
                                inputId = "SO", 
                                label = "Choose SO Number",
                                c(Choose="", SO_options),
                                
                                selectize=TRUE),
                        actionButton("refresh","Refresh SO Choices"),
                        actionButton("newimages","Find New Images"),
                        selectInput(
                                inputId = "QCB", 
                                label = "Input QC Barcode",
                                c(Choose="", QCB_options),
                                selectize = TRUE,
                                multiple = TRUE),
                        actionButton("goButton", "Parse Data"),
                        downloadButton('download', 'Download')
                ),
                
                mainPanel(
                         uiOutput("plots"),
                         textOutput("message"),
                         tableOutput("reports")
                        )
                        
                )
        )
)


