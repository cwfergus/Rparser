library(stringr)
library(calibrate)
library(dplyr)

localMaxima <- function(x) {
        # Use -Inf instead if x is numeric (non-integer)
        y <- diff(c(-.Machine$integer.max, x)) > 0L
        rle(y)$lengths
        y <- cumsum(rle(y)$lengths)
        y <- y[seq.int(1L, length(y), 2L)]
        if (x[[1]] == x[[2]]) {
                y <- y[-1]
        }
        y
}

coa_plot <- function(x,y,x_txt,y_txt,l_txt,pt_b,pt_c,fb,gr_t,labels,l_color="black", shd=NULL, so)
{
        # x,y are vectors; pt_b TRUE/FALSe to show point txt; pt_min: min value to show pt txt
        # fb = file base name; gr_t = mz or lc plot type; labels includes titles, axis labels, and other specs
                        # filename
        # Open a device: bmp & pdf are other options here
        plot(x, y, 
             main=labels[1], sub=labels[2],xlab=labels[6], ylab=labels[7], type="n", axes=TRUE,
             col.lab="blue", col.axis="blue"
        )
        if(pt_b) textxy(x_txt, y_txt, l_txt, cex = 1, offset=-0.8)
        if(pt_c) {
                textxy(x_txt, 0.96*y_txt, l_txt, cex=1, offset=-0.8)
                polygon(shd, density=NA, col="green", border="black")
                text (min(x), max(y), labels[9], pos=4)
        }
        lines(x,y, lty=1, type=labels[8], pch='', col=l_color)
        
        text(max(x),max(y), labels[3], pos=2)
        text(max(x),0.96*max(y), labels[4], pos=2)
        text(max(x),0.92*max(y), labels[5], pos=2)
        
        # Close the image device that was opened above; the actions in the meantime have been recorded
        return(fn)
}

curvecolor <- function(ploc, tdf, df, range) {
        
        
        shd_bl <- grep(df[ploc+4, 2], tdf[,1], fixed=TRUE)
        shd_el <- grep(df[ploc+5,1], tdf[,1], fixed=TRUE)
        if (length(shd_bl) == 0) {
                bt_m1 <- type.convert(df[ploc+4,2])-0.0001
                shd_bl_m1 <- grep(bt_m1, tdf[,1], fixed=TRUE)
                bt_p1 <- type.convert(df[ploc+4,2])+0.0001
                shd_bl_p1 <- grep(bt_p1, tdf[,1], fixed=TRUE)
                if (length(shd_bl_m1) > 0) shd_bl <- shd_bl_m1
                if (length(shd_bl_p1) > 0) shd_bl <- shd_bl_p1
        }
        if (length(shd_el) == 0) {
                et_m1 <- type.convert(df[ploc+5,1])-0.0001
                shd_el_m1 <- grep(et_m1, tdf[,1], fixed=TRUE)
                et_p1 <- type.convert(df[ploc+5,1])+0.0001
                shd_el_p1 <- grep(et_p1, tdf[,1], fixed=TRUE)
                if (length(shd_el_m1) > 0) shd_el <- shd_el_m1
                if (length(shd_el_p1) > 0) shd_el <- shd_el_p1
                
        }
        
        shd_b <- min(shd_bl)
        shd_e <- max(shd_el)
        shd <- type.convert(tdf[shd_b:shd_e,])
        shd[1,2] <- 1
        shd[length(shd[,1]),2] <- 1
        shd[,2] <- (shd[,2]/100)*range
        return (shd)
}

data_path <- "//DATA4/Mass Spec DataRrepository/reportfiles/current_report_files/"
data_pathrpt <- paste(data_path, ".rpt", sep="")
data <- data.frame()

shinyServer(function(input, output){
        observe( {
                filename_reactive <- eventReactive(input$goButton, {
                        input$QCB
                })
                filename <- filename_reactive()
                filenamefull <- paste(data_path, filename, sep="")
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
                
                
                
                ##### Output file generation ####
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
                
               
#         
#         
        #Takes f2, which is the location of the sample like 16:4,F, and extracts the
        #different parts of this.
        
        #Remove the -01 etc from the ProdNum to get the SSID
        
        #Finds the date, converts to date class, then converts back to character string
        # in the format that filemaker wants
        
        # Splits data by SAMPLE, enabling looking at each sample for each value
        })
        
})
