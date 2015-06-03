library(shiny)

shinyServer(function(input, output, session) {
        
        observe({
                
                if (input$browse == 0) return()
                
                updateTextInput(session, "path",  value = file.choose())
        })
        
        contentInput <- reactive({ 
                
                if(input$upload == 0) return()
                
                isolate({
                        data <- read.table(input$path,sep="\t", col.names = c("SO", "SS"), stringsAsFactors=FALSE)
                        splitdata <- split(data, data$SO)
                        
                        SO_list <- unique(data$SO)
                        message <- data.frame(SO="SO", suc_move="Images Processed?", stringsAsFactors=FALSE)
                        for (i in 1:length(SO_list)){
                                SO <- SO_list[i]
                                img_loc <- paste(input$path, "imgs", SO, sep="/")
                                ship_loc <- paste(input$path, "Shipping", SO, sep="/")
                                SOchr <- as.character(SO)
                                tempdata <- splitdata[[SOchr]]
                                SS_list_lc <- paste(unique(tempdata$SS), "_lc.png", sep="")
                                SS_list_mz <- paste(unique(tempdata$SS), "_mz.png", sep="")
                                SS_list <- c(SS_list_lc, SS_list_mz)
                                list.files(img_loc) -> all_imgs
                                all_imgs[(all_imgs%in%SS_list)] -> print_these
                                from_print_paths <- paste(img_loc, print_these, sep="/")
                                to_print_paths <- paste(ship_loc, print_these, sep="/")
                                if (!file.exists(ship_loc)){
                                        dir.create(ship_loc)
                                }
                                file.copy(from = from_print_paths, to = to_print_paths, overwrite=FALSE) -> temp_msg
                                printed_loc <- paste(input$path,"/imgs/", "printed_", SO, sep="")
                                file.rename(img_loc, printed_loc)
                                suc_move <- as.character(unique(temp_msg))
                                if (length(suc_move)==2) {
                                        suc_move <- "Partially"
                                } 
                                row <- cbind(SO, suc_move)
                                message <- rbind(message, temp_msg)
                        }
                        message
                        })
        })
        
        output$content <- renderPrint({
                contentInput()
        })
        
})