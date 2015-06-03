library(shiny)
library(shinyAce)
data_path <- "//DATA4/Mass Spec DataRrepository/imgs/"
data_path_rpt <- "//DATA4/Mass Spec DataRrepository/reportfiles/current_report_files/"
SO_options <- list.dirs(data_path, full.names = FALSE )
QCB_options <- list.files(data_path_rpt, pattern = "*.rpt", full.names = FALSE )

shinyUI(navbarPage("QC data lookup",
                   tabPanel("Select Images",
                            sidebarLayout(
                                    sidebarPanel(
                                            includeText("dir.txt"),
                                            br(),
                                            br(),
                                            selectInput(
                                                    inputId = "SO", 
                                                    label = "Choose SO Number",
                                                    c(Choose="", SO_options),
                                                    selectize=TRUE),
                                            uiOutput("ui"),
                                            actionButton("refresh","Refresh SO Choices"),
                                            br(),
                                            br(),
                                            includeText("emaildir.txt"),
                                            br(),
                                            br(),
                                            textInput("to", "To:", "example@anything.com"),
                                            textInput("subject", "Subject:", value="Your Traces"),
                                            aceEditor("message", value="Enter your message here", height="200px"),
                                            actionButton("send", "Send Email")
                                    ),
                                    
                                    mainPanel(
                                            tabsetPanel(
                                                    tabPanel("MS Trace",plotOutput("mztrace2")),
                                                    tabPanel("LC Trace",plotOutput("lctrace2")),
                                                    tabPanel("Both Traces", plotOutput("mztrace"),
                                                             plotOutput("lctrace")),
                                                    tabPanel("Flr Trace", plotOutput("fltrace"))
                                                    
                                            )
                                            
                                    ))
                           ),
                   tabPanel("All Images",
                            sidebarLayout(
                                    sidebarPanel(
                                            selectInput(
                                                    inputId = "SO2", 
                                                    label = "Choose SO Number",
                                                    c(Choose="", SO_options),
                                                    
                                                    selectize=TRUE),
                                            actionButton("refresh2","Refresh SO Choices"),
                                            tags$br(),
                                            textInput("dlpath", "Download Path"),
                                            actionButton("browse", "Set Download Path"),
                                            actionButton("downloadplots", "Download Plots")
                                    ),
                                    mainPanel(uiOutput("plots"))
                                    
                            )
                            
                            
                   ),
                   tabPanel("Generate Images",
                            sidebarLayout(
                                    sidebarPanel(
                                            selectInput(
                                                    inputId = "QCB", 
                                                    label = "Input QC Barcode",
                                                    c(Choose="", QCB_options),
                                                    selectize = TRUE,
                                                    multiple = TRUE),
                                            actionButton("goButton", "Parse Data"),
                                            downloadButton('parseddata', 'Download')
                                    ),
                                    mainPanel(
                                            textOutput("message"),
                                            tableOutput("reports")
                                    )
                            )
                   )
                   
                   
                   
                   
)
)


