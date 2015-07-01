library(shiny)
data_path <- "//DATA4/Mass Spec DataRrepository/reportfiles/current_report_files/"
QCB_options <- list.files(data_path, pattern = "*.rpt", full.names = FALSE )

shinyUI(fluidPage(
        titlePanel("QC data lookup"),
        
        sidebarLayout(
                sidebarPanel(
                        selectInput(
                                inputId = "QCB", 
                                label = "Input QC Barcode",
                                c(Choose="", QCB_options),
                                selectize = TRUE,
                                multiple = TRUE),
                        br(),
                        actionButton("goButton", "Parse Data"),
                        br(),
                        downloadButton('download', 'Download')
                ),
                
                mainPanel(
                        textOutput("message"),
                        tableOutput("reports")
                )
        )
))
