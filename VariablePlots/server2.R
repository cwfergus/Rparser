library(stringr)
library(sendmailR)
data_path <- "//DATA4/Mass Spec DataRrepository/imgs/"
lc_end <- "_lc.png"
mz_end <- "_mz.png"
mailControl=list(smtpServer="EXCHANGE2.bti.biosearchtech.com")

shinyServer(function(input, output, session){
        
        output$ui <- renderUI({
                fpath <- paste(data_path, input$SO, "/", sep="")
                flist <- list.files(fpath, pattern="*.png")
                sslist <- unique(str_sub(flist, end=-8))
                selectInput(
                        inputId = "SS", 
                        label = "Choose SS Number",
                        c(Choose ="", sslist),
                        selectize=TRUE)
        })
        
        output$plots <- renderUI({
                plot_output_list <- lapply(1:length(flist), function(i){
                        plotname <- flist[i]
                        plotOutput(plotname)
                })
                do.call(taglist, plot_output_list)
        })
        
        for (i in 1:length(flist)) {
                local({
                        my_i <- i
                        plotname <- flist[my_i]
                        output[[plotname]] <- renderImage({
                                fn <- paste(data_path, input$SO, "/", plotname, sep="")
                                list(src = fn,
                                     contentType = 'image/png',
                                     width = 680,
                                     height = 440)
                        }, deletefile=FALSE)
                        })
                })
        }
        
        
        
        output$mztrace <- renderImage({
                fn <- paste(data_path, input$SO,"/", input$SS, mz_end, sep="")
                list(src = fn,
                     contentType = 'image/png',
                     width = 680,
                     height = 440,
                     alt = "Choose an SO and SS to begin")
        }, deleteFile=FALSE)
        output$mztrace2 <- renderImage({
                fn <- paste(data_path, input$SO,"/", input$SS, mz_end, sep="")
                list(src = fn,
                     contentType = 'image/png',
                     width = 680,
                     height = 440,
                     alt = "Choose an SO and SS to begin")
        }, deleteFile=FALSE)
        output$lctrace <- renderImage({ 
                fn <- paste(data_path, input$SO, "/", input$SS, lc_end, sep="")
                list(src = fn,
                     contentType = 'image/png',
                     width = 680,
                     height = 440,
                     alt = "")
        }, deleteFile=FALSE)
        output$lctrace2 <- renderImage({ 
                fn <- paste(data_path, input$SO, "/", input$SS, lc_end, sep="")
                list(src = fn,
                     contentType = 'image/png',
                     width = 680,
                     height = 440,
                     alt = "")
        }, deleteFile=FALSE)
        
        observe({
                if(is.null(input$send) || input$send==0) return(NULL)
                from <- "qcdatacheck@biosearchtech.com"
                to <- isolate(input$to)
                subject <- isolate(input$subject)
                msg <- isolate(input$message)
                MZattachment_Path <- isolate(paste(data_path, input$SO,"/", input$SS, mz_end, sep=""))
                MZattachment_Name <- isolate(paste("SO: ",input$SO," ", input$SS, mz_end, sep=""))
                MZattachment <- mime_part(x=MZattachment_Path, name=MZattachment_Name)
                LCattachment_Path <- isolate(paste(data_path, input$SO, "/", input$SS, lc_end, sep=""))
                LCattachment_Name <- isolate(paste("SO: ",input$SO, "/", input$SS, lc_end, sep=""))
                LCattachment <- mime_part(x=LCattachment_Path, name=LCattachment_Name)
                if (MZattachment_Path == "//DATA4/Mass Spec DataRrepository/imgs//_mz.png" |
                            LCattachment_Path == "//DATA4/Mass Spec DataRrepository/imgs//_lc.png"){
                        body <- msg
                } else {body <- list(msg, MZattachment, LCattachment)}
                
                sendmail(from, to, subject, body, control=mailControl)
        })
        observe({
                if(is.null(input$refresh) || input$refresh==0) return(NULL)
                SO_choices2 <- list.dirs(data_path, full.names = FALSE )
                updateSelectizeInput(session, "SO", "Choose SO Number", 
                                     choices = c(Refreshed="", SO_choices2))
                
        })
        
        session$onSessionEnded(function(){
                stopApp()
        })
})
