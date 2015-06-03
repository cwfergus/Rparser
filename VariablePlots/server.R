library(stringr)
library(sendmailR)
data_path <- "//DATA4/Mass Spec DataRrepository/imgs/"
data_path_rpt <- "//DATA4/Mass Spec DataRrepository/reportfiles/current_report_files/"
lc_end <- "_lc.png"
mz_end <- "_mz.png"
fl_end <- "_fl.png"
mailControl=list(smtpServer="EXCHANGE2.bti.biosearchtech.com")
data <- data.frame()

shinyServer(function(input, output, session){
        ##### All Images Tabs #####
        observe({
                
                if (input$browse == 0) return()

                        updateTextInput(session, "dlpath",  value = choose.dir(default = "C:/Users/", caption="Select Save Location"))
                
        })
        observe({
                if(input$SO2=="") return(NULL)
                isolate({fpath <- paste(data_path, input$SO2, "/", sep="")
                         flist <- list.files(fpath, pattern="*.png")
                         sslist <- unique(str_sub(flist, end=-8))
                         output$plots <- renderUI({
                                 plot_output_list <- lapply(1:length(flist), function(i){
                                         plotname <- flist[i]
                                         plotOutput(plotname)
                                 })
                                 do.call(tagList, plot_output_list)
                         })
                         
                         for (i in 1:length(flist)) {
                                 local({
                                         my_i <- i
                                         plotname <- flist[my_i]
                                         output[[plotname]] <- renderImage({
                                                 fn <- paste(data_path, input$SO2, "/", plotname, sep="")
                                                 list(src = fn,
                                                      contentType = 'image/png',
                                                      width = 680,
                                                      height = 440)
                                         }, deleteFile=FALSE)
                                 })
                         }
                         
                         })})
        observe({
                if(is.null(input$refresh2) || input$refresh2==0) return(NULL)
                SO_choices2 <- list.dirs(data_path, full.names = FALSE )
                updateSelectizeInput(session, "SO2", "Choose SO Number", 
                                     choices = c(Refreshed="", SO_choices2))
                
        }, priority=1)
        
        observe({
                if(is.null(input$downloadplots) || input$downloadplots==0) return(NULL)
                isolate({fpath <- fpath <- paste(data_path, input$SO2, "/", sep="")
                         flist <- list.files(fpath, pattern="*.png", full.names=TRUE)
                         file.copy(from = flist, to = input$dlpath)
                         })
        })
        ##### Image Generater ####
        observe({
                filename_reactive <- eventReactive(input$goButton, {
                        input$QCB
                })
                filename <- filename_reactive()
                filenamefull <- paste(data_path_rpt, filename, sep="")
                for (i in 1:length(filenamefull)) {
                        tempdata <- scan(filenamefull[i],
                                         what= list(Type="", Value = ""), #Reads into two columns as a list
                                         sep = "\t", #determines the two columns by the seperater tab
                                         multi.line = TRUE, #does not break 
                                         fill=TRUE,
                                         comment.char = "{",
                                         allowEscapes = FALSE)
                        tempdata <- data.frame(do.call('cbind', tempdata),
                                               stringsAsFactors = FALSE)
                        data <- rbind(data, tempdata) 
                }
                data$Value[grep("Well", data$Type)] -> f2 
                data$Value[grep("UserName", data$Type)] -> UserName
                data$Value[grep("^Instrument$", data$Type)] -> Instrument
                data$Value[grep("^SampleID$", data$Type)] -> ProdNum
                data$Type[grep("Formula", data$Type)+1] -> ExpectMass
                data$Value[grep("SampleDescription", data$Type)] -> SO
                data$Value[(grep("JobCode", data$Type))]-> ESBarCode
                #Takes f2, which is the location of the sample like 16:4,F, and extracts the
                #different parts of this.
                str_extract(f2, "[[:digit:]]+") -> PlateLoc
                str_extract(str_extract(f2, "[,]{1}[[:alpha:]]"), "[[:alpha:]]+") -> Row
                str_extract(str_extract(f2, "[[:digit:]]+[,]"), "[[:digit:]]+") -> Col
                #Put the row and col together like Filemaker wants it.
                paste(Row, Col, sep = ":") -> Well
                
                #Remove the -01 etc from the ProdNum to get the SSID
                str_extract(ProdNum, "[[:digit:]]+") -> SSID
                
                #Finds the date, converts to date class, then converts back to character string
                # in the format that filemaker wants
                Date <- data$Value[grep("^Date$", data$Type)] %>%
                        as.Date(format="%d-%b-%Y") %>%
                        format(format = "%m/%d/%Y")
                
                #Find Time values with a :, then converts to Date/Time format, then back to 
                # character string in Time format FM wants.
                Time <- grep(":", data$Value[grep("^Time$", data$Type)], value=TRUE) %>%
                        strptime(format="%X") %>%
                        format(format="%r")
                
                # Splits data by SAMPLE, enabling looking at each sample for each value.
                samplesplit <- split(data, cumsum(data[,"Type"] == "[SAMPLE]"))
                
                #The loop looks within each Sample (now sperated in a list), and finds the BPM row
                #IT then collects the value next to it. If there is no value found, it makes it a 0
                #Then it generates a vector of Found masses
                # If it can't find a mass, then it puts 0 into its place. This solves dead wells.
                Results <- vector()
                for (i in 1:length(samplesplit)) {
                        workingdata <- as.data.frame(samplesplit[i])
                        colnames(workingdata) <- c("Type", "Value")
                        Foundmass <- workingdata$Value[grep("BPM", workingdata$Type)]
                        if (length(Foundmass) == 0) {
                                Foundmass <- 0
                        } 
                        Results <- append(Results, Foundmass)
                }
                
                #The loop does a similar thing but for Area. It also looks for FLR data being
                #available. If FLR data exists, it splits up the sample data into chromatogram
                #split data. Then it finds the Area values for both FLR and UV. 
                #If the length of the split is 1 then no area values exists, and it appends 0 for the
                #sample.
                Area <- vector()
                FlrArea <- vector()
                for (i in 1:length(samplesplit)) {
                        workingdata <- as.data.frame(samplesplit[i])
                        colnames(workingdata) <- c("Type", "Value")
                        chromosplit <- split(workingdata, cumsum(workingdata[,"Type"]== "[CHROMATOGRAM]"))
                        if (length(chromosplit) > 2) {
                                FLRlocations <- grep("ACQUITY FLR", chromosplit)
                                for (x in 1:length(FLRlocations)) {
                                        flr <- FLRlocations[x]
                                        workingdata2 <- as.data.frame(chromosplit[flr])
                                        colnames(workingdata2) <- c("Type", "Value")
                                        percentareas2 <- as.numeric(workingdata2$Value[grep("Area %Total", workingdata2$Type)])
                                        FlrArea <- append(FlrArea, max(percentareas2))
                                }
                                if (length(chromosplit) == 1) {
                                        Area <- append(Area, 0)
                                }
                                workingdata3 <- as.data.frame(chromosplit[2])
                                colnames(workingdata3) <- c("Type", "Value")
                                percentareas3 <- as.numeric(workingdata3$Value[grep("Area %Total", workingdata3$Type)])
                                Area <- append(Area, max(percentareas3))
                        } else {
                                if (length(chromosplit) == 1) {
                                        Area <- append(Area, 0)
                                } else{
                                        percentareas <- as.numeric(workingdata$Value[grep("Area %Total", workingdata$Type)])
                                        Area <- append(Area, max(percentareas)) }}
                }
                
                #Generates a Test Type vector, which is just the TestType repeated over and over
                rep_len(x = "MS", length.out = length(f2)) -> TestType_MS
                rep_len(x = "HPLC", length.out = length(f2)) -> TestType_HPLC
                
                
                
                #Now with all the data generated we have to put it together in the format that FM wants
                
                #Makes vectors with empty fields for each of the two TestType splits
                ExpectMass_null <- vector(mode = "character", length = length(ExpectMass))
                Results_null <- vector(mode = "character", length = length(Results))
                Area_null <- vector(mode = "character", length = length(Area))
                
                #Make a vector of column names just to make it easy.
                columnnames <- c("Well",
                                 "ProdNum",
                                 "Date",
                                 "Time",
                                 "UserName",
                                 "Instrument",
                                 "Results",
                                 "PlateLoc",
                                 "Row",
                                 "Col",
                                 "SSID",
                                 "ExpectMass",
                                 "TestType",
                                 "SO",
                                 "ESBarCode",
                                 "Area")
                #Bind together the different vectors into a data frame.
                msESMS <- data.frame(cbind(Well,
                                           ProdNum,
                                           Date,
                                           Time,
                                           UserName,
                                           Instrument,
                                           Results,
                                           PlateLoc,
                                           Row,
                                           Col,
                                           SSID,
                                           ExpectMass,
                                           TestType_MS,
                                           SO,
                                           ESBarCode,
                                           Area_null))
                
                areaESMS <- data.frame(cbind(Well,
                                             ProdNum,
                                             Date,
                                             Time,
                                             UserName,
                                             Instrument,
                                             Results_null,
                                             PlateLoc,
                                             Row,
                                             Col,
                                             SSID,
                                             ExpectMass_null,
                                             TestType_HPLC,
                                             SO,
                                             ESBarCode,
                                             Area))
                
                #Append the new column names, removing the Area_null name to Area
                colnames(msESMS) <- columnnames
                colnames(areaESMS) <- columnnames
                #puts the two TestType data frames together
                rbind(msESMS, areaESMS) -> likeESMS
                output$message <- renderPrint("You have successfully parsed the following reports:")
                filenametable <- as.data.frame(filename)
                colnames(filenametable) <- "Parsed Reports"
                output$reports <- renderTable(filenametable)
                output$download <- downloadHandler(
                        filename = "dataforimport.txt",
                        content = function(file) {
                                write.table(likeESMS, file, sep = "\t", col.names=FALSE, quote=FALSE, row.names=FALSE, na="")
                        })
                
        })         
        
  #####Select Images####              

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
        
        
        
        output$mztrace <- renderImage({
                fn <- paste(data_path, input$SO,"/", input$SS, mz_end, sep="")
                list(src = fn,
                     contentType = 'image/png',
                     width = 680,
                     height = 440,
                     alt = "MZ Trace not available")
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
                     alt = "LC trace not available")
        }, deleteFile=FALSE)
        output$lctrace2 <- renderImage({ 
                fn <- paste(data_path, input$SO, "/", input$SS, lc_end, sep="")
                list(src = fn,
                     contentType = 'image/png',
                     width = 680,
                     height = 440,
                     alt = "LC trace not available")
        }, deleteFile=FALSE)
        output$fltrace <- renderImage({
                fn <- paste(data_path, input$SO, "/", input$SS, fl_end, sep="")
                list(src = fn,
                     contentType = 'image/png',
                     width = 680,
                     height = 440,
                     alt = "FLR trace not available")
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
