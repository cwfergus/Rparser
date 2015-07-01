library(stringr)
library(xlsx)
data_path_rpt <- "//DATA4/Mass Spec DataRrepository/reportfiles/current_report_files/"
data <- data.frame()

shinyServer(function(input, output, session){
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
                data$Value[grep("^SampleID$", data$Type)] -> ProdNum
                data$Type[grep("Formula", data$Type)+1] -> ExpectMass
                
                samplesplit <- split(data, cumsum(data[,"Type"] == "[SAMPLE]"))
                
                #The loop looks within each Sample (now sperated in a list), and finds the BPM row
                #IT then collects the value next to it. If there is no value found, it makes it a 0
                #Then it generates a vector of Found masses
                # If it can't find a mass, then it puts 0 into its place. This solves dead wells.
                #Generates empty vectors for the results to go into
                Results <- vector()
                Area <- vector()
                #The following loop finds all the Test data. It looks at each sample individually and
                #generates data frames based on the type of data. Then its searches for the test result
                # If a specific test is missing, or is shorter than 20 lines (also means missing data)
                for (ss in 1:length(samplesplit)) {
                        workingdata <- as.data.frame(samplesplit[ss]) #Just 1 samples data
                        colnames(workingdata) <- c("Type", "Value") 
                        functionsplit <- split(workingdata, cumsum(workingdata[,"Type"]=="[FUNCTION]")) #split data by function
                        for (f in 2:length(functionsplit)) { #skip the first section, this is an intro
                                if (nrow(as.data.frame(functionsplit[f])) <= 20) {
                                        next #If the function contains less that 20 rows, ignore it
                                }
                                if (any(grepl("MS", functionsplit[f]))){ #Is this MS data?
                                        ms_data <- as.data.frame(functionsplit[f]) #put in data frame
                                        colnames(ms_data) <- c("Type", "Value")
                                }
                                if (any(grepl("DAD", functionsplit[f]))){ #Is it LC data?
                                        dad_data <- as.data.frame(functionsplit[f])
                                        colnames(dad_data) <- c("Type", "Value")
                                }
                        }
                        
                        ## Find Mass value:
                        if (exists("ms_data")) { #If Ms data exists, find the mass
                                Foundmass <- ms_data$Value[grep("BPM", ms_data$Type)]
                        } else { #Otherwise, put 0 for the mass and report missing data
                                Foundmass <- 0
                        } #append the sample result to all results
                        Results <- append(Results, Foundmass) 
                        
                        ## Find LC Area
                        if (exists("dad_data")) { #If LC data exists, find Max Area %
                                percentareas <- as.numeric(dad_data$Value[grep("Area %Total", dad_data$Type)])
                                Foundarea <- max(percentareas)
                        } else { #Otherwise put 0 for Area and report missing data
                                Foundarea <- 0
                        } #Append result to all results
                        Area <- append(Area, Foundarea)
                        
                }
                
                finaldata <- data.frame(cbind(ProdNum, Area, Results, ExpectMass))
                finaldata <- finaldata[c(2:nrow(finaldata)),]
                filename <- paste(str_sub(ProdNum[2], end=-4), "purity.xlsx", sep=" ")
                output$reports <- renderTable(finaldata)
                output$download <- downloadHandler(
                        filename =  filename,
                        content = function(file) {
                                write.xlsx(finaldata, file, row.names=FALSE)
                        })
        })
        
        session$onSessionEnded(function(){
                stopApp()
        })
})
        
