##### Script Description ####
# The following script reads a MassLynx LC/MS report file
# It then parses out sample information, It grabs all the information shown in the ESMS import
# script final screen, except that which is not in the rpt file.
# It then compiles a data frame with each sample as a row
# It then writes out a tab delimited text file

#####Data Read in ####
library(stringr)
library(MALDIquant)

data <- data.frame()

filename <- "y"
fn_list <- "You have successfully parsed the following reports:"
fn_full_list <- vector()

while (filename != "n") {
        filename <- readline("Please enter the QC set Barcode. 
                             If you are done adding reports enter n\t")
        filenamerpt <- paste(filename, ".rpt", sep="")
        filenamereprocess <- paste(filename, "_reprocess.rpt", sep="")
        
        data_path <- "//DATA4/Mass Spec DataRrepository/reportfiles/current_report_files/"
        
        filenamefull <- paste(data_path, filenamerpt, sep="")
        
        reprocess_filename <- paste(data_path, filename, "_reprocess.rpt", sep="")
        
        if (file.exists(reprocess_filename) == TRUE) {
                print("A more up-to-date version of this data may be availabe at:")
                print(filenamereprocess)
                decision <- readline("Do you want to import this instead? Enter y for yes\t")
                if (decision == "y") {
                        filenamefull <- reprocess_filename
                        filenamerpt <- filenamereprocess
                }
        }
        breakname <- paste(data_path, "n.rpt", sep="")
        if (filenamefull == breakname) {
                break
        }
        
        if (file.exists(filenamefull) == TRUE) {
                
                if (filenamefull %in% fn_full_list == TRUE) {
                        print("You have already parsed this report")
                } else {
                        tempdata <- scan(filenamefull, 
                                         what= list(Type="", Value = ""), #Reads into two columns as a list
                                         sep = "\t", #determines the two columns by the seperater tab
                                         multi.line = TRUE, #does not break 
                                         fill=TRUE,
                                         comment.char = "{",
                                         allowEscapes = FALSE)
                        tempdata <- data.frame(do.call('cbind', tempdata),
                                               stringsAsFactors = FALSE)
                        data <- rbind(data, tempdata) 
                        fn_full_list <- append(fn_full_list, filenamefull)
                        fn_list <- append(fn_list, filenamerpt)
                }
                
        } else {
                print("No file exists with the following name!")
                print(filenamefull)
        }
}

##### Parser ####

#Searches in column Type for the word Well, and returns the value in the column next to it
#Each of these commands makes a vector with all of the values for each sample.
data$Value[grep("Well", data$Type)] -> f2 #name chosen due to likeness with ESMS import
# Repeat the word well, for the number of samples, as determined by the pevious find
rep_len(x = "Well", length.out = length(f2)) -> f1

#Find the UserName row, take the info in the Value Column
data$Value[grep("UserName", data$Type)] -> UserName
#Find the Instrument row. The search uses a Regular expression for find the word Instrument,
#with nothing before it, and nothing after it.
data$Value[grep("^Instrument$", data$Type)] -> Instrument
#Str_extract extracts the first digit in the well location, which is actuall the Plate location
str_extract(f2, "[[:digit:]]+") -> PlateLoc
#Extract just the row, which is the first letter after a comma. Then remove the comma
str_extract(str_extract(f2, "[,]{1}[[:alpha:]]"), "[[:alpha:]]+") -> Row
#Extract just the column, which is the first number before the comma, then remove the comma
str_extract(str_extract(f2, "[[:digit:]]+[,]"), "[[:digit:]]+") -> Col
#Put the row and column together seperated by a : to represent the Well, w/o the plate location
paste(Row, Col, sep = ":") -> Well


#Synthesis Batch: found within FM, not possible in R

# SSID found within FM? I can probably do this here if necessary

#Searches for the SampleID row, finds the Production number, makes a vector with all of them
data$Value[grep("^SampleID$", data$Type)] -> ProdNum

str_extract(ProdNum, "[[:digit:]]+") -> SSID

#Searches for the Date row, makes vector with all the Dates
data$Value[grep("^Date$", data$Type)] -> Date
#Converts Dates to the actual R class for dates
as.Date(Date, format="%d-%b-%Y") -> Date
#Converts it out of this R format to the Filemaker format
format(Date, format = "%m/%d/%Y") -> Date

#Searches for the Time row, makes vector of times
data$Value[grep("^Time$", data$Type)] -> Time
#Removes Times that arn't actuall the time stamp, by selecting only those with :
grep(":", Time, value=TRUE) -> Time
#Converts to R type format then back into Filemaker format
format(strptime(Time, format = "%X"), format = "%r") -> Time

#Searches for the Expected mass, by finding the Word Formulat, then moving down 1 row
data$Type[grep("Formula", data$Type)+1] -> ExpectMass

#Splits the Data into sample groups, using the string [SAMPLE]. This creates a 
#list of data frames, each item being a seperate sample
samplesplit <- split(data, cumsum(data[,"Type"] == "[SAMPLE]"))
#Initiate an empty Results Vector
Results <- vector()
#The loop looks within each Sample (now sperated in a list), and finds the BPM row
#IT then collects the value next to it. If there is no value found, it makes it a 0
#Then it generates a vector of Found masses
for (i in 1:length(samplesplit)) {
        workingdata <- as.data.frame(samplesplit[i])
        colnames(workingdata) <- c("Type", "Value")
        Foundmass <- workingdata$Value[grep("BPM", workingdata$Type)]
        if (length(Foundmass) == 0) {
                Foundmass <- 0
        } 
        Results <- append(Results, Foundmass)
}

#The loop does a similar thing but for Area
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


#Finds the SO numbers
data$Value[grep("SampleDescription", data$Type)] -> SO

#Customer: in filemaker not in input

#SeqID: in filemaker not in input



#Finds the Barcode numbers
data$Value[(grep("JobCode", data$Type))]-> ESBarCode

##### Output file generation ####

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
#Bind together the different vectors into a data frame. Here we make the MS Testtype DF
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

#Append the new column names, removing the Area_null name to Area
colnames(msESMS) <- columnnames
#Bind together the HPLC test type stuff
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
#Remove the ExpectMass_null, and Results_null
colnames(areaESMS) <- columnnames
#puts the two TestType data frames together
rbind(msESMS, areaESMS) -> likeESMS
#Writes out the now ESMS replica data to the outputnume
write.table(likeESMS, "dataforimport.txt", sep = "\t", col.names=FALSE, quote=FALSE, row.names=FALSE, na="")
print(fn_list)
print("Don't forget to now import into Filemaker")

#####Graphing####
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

coa_plot <- function(x,y,x_txt,y_txt,l_txt,pt_b,fb,gr_t,labels,l_color="black")
{
        # x,y are vectors; pt_b TRUE/FALSe to show point txt; pt_min: min value to show pt txt
        # fb = file base name; gr_t = mz or lc plot type; labels includes titles, axis labels, and other specs
        fn <- paste("imgs/",fb,"_",gr_t,".png",sep="")                 # filename
        png(filename=fn, width=1400, height=500)	# Open a device: bmp & pdf are other options here
        if(pt_b) {
                peaks <- createMassPeaks(x, y)
                plot(peaks, col="red")
                labelPeaks(peaks, index=intensity(peaks)>25, adj=2, arrowLength=0.25, arrowLwd=1, arrowCol="blue")
                
        } else {
                plot(x, y, sub=labels[2], type="n", axes=TRUE) 
                textxy(x_txt, y_txt, l_txt, cex=1, offset=-0.8)
                lines(x,y, lty=1, type=labels[8], pch='', col=l_color)
        }
               
        
        
        text(max(x),max(y), labels[3], pos=2)
        text(max(x),0.96*max(y), labels[4], pos=2)
        text(max(x),0.92*max(y), labels[5], pos=2)
        dev.off()
        # Close the image device that was opened above; the actions in the meantime have been recorded
        return(fn)
}

ts_0 <- Sys.time()
fn_list <- ""        # brings in entire rpt file as list type; not the most efficient thing to do: seek alt method; also, is there a limit to amount of txt that can be processed this way?
grf <- grep("[SAMPLE]", data$Type, fixed=TRUE)        # gives vector with line numbers where each SAMPLE data starts
imax <- length(grf)

# for each SAMPLE in the file, do this whole section:
for (i in 1:imax)
{
        # find only the data belonging to the current sample
        io = grf[i] #Line number at beginning of sample
        i_n <- if(i==imax) length(data$Value) else grf[i+1]-1 #Line at end of sample
        txti <- cbind(data$Type[io:i_n], data$Value[io:i_n]) #All the lines in a sample.
        
        # find specific data sets within the sample
        g_mzi <- grep("[MS]", txti[,1], fixed=TRUE)	# start of mass spec data
        g_lci <- grep("[TRACE]", txti[,1], fixed=TRUE)# start of hplc data
        g_spi <- grep('}', txti[,1], fixed=TRUE)		# vector of all end-bracket locations
        g_pks <- grep("[PEAK]", txti[,1], fixed=TRUE)		# vector of all PEAK locations
        if (length(g_lci) == 0 | length(g_mzi) ==0)  { 
                err_msg <- paste (ProdNum[i], "is missing LC or MS data")
                print(err_msg) }
        else {
                # analyze mass spec data
                m_mzi <- g_spi[min(which(g_spi>g_mzi))]-1	# get the lowest position of } that's after [MS]
                t_mzi <- txti[(g_mzi+2):m_mzi,]			# get all data between [MS]{ and }
                x_mzi <- type.convert(t_mzi[,1])   		# every other point is x; make it numeric
                y_mzi <- type.convert(t_mzi[,2])     	# every other point is y; make it numeric
                # calculate text labels for mass spec; only the points that will be printed
                y_min <- 0.25*max(y_mzi)     # must be over minimum 5% threshold
                i_min <- which(y_mzi>y_min)
                X_min <- x_mzi[i_min]
                y_min <- y_mzi[i_min]
                i_txt <- localMaxima(y_min)  # only peak maxima instead of clusters of points, doesn't need change.
                x_txt <- X_min[i_txt]
                y_txt <- y_min[i_txt]
                l_txt <- as.character(x_txt)
                
                l# analyze hplc data
                m_lci <- g_spi[min(which(g_spi>g_lci))]-1	# get the lowest position of } that's after [TRACE]
                t_lci <- txti[(g_lci+1):m_lci,]			# get all data between [TRACE]{ and }
                x_lci <- type.convert(t_lci[,1])     	# every other point is x; make it numeric
                y_lci <- type.convert(t_lci[,2])     	# every other point is y; make it numeric
                # calculate text labels for LC
                peak_loc <- grep("[PEAK]", txti[,1], fixed=TRUE)
                areatable <- data.frame()
                for (l in 1:length(peak_loc)) {
                        rt <- type.convert(txti[peak_loc[l]+3,2])
                        area <- type.convert(txti[peak_loc[l]+11,2])
                        table <- cbind(rt, area)
                        areatable <- rbind(areatable,table)
                }
                y_lci_min <- 5
                minvalue <- areatable[which(areatable$area==max(areatable$area)),]
                x_lci_min <- minvalue$rt
                y_lci_min <- minvalue$area
                y_lci_min_txt <- paste("Area % Purity:", minvalue$area, "%")
                
                # analyze sample header data
                t_hdr <- txti[1:g_mzi,]				# everything before mass spec data is header text
                dfb   <- ProdNum[i]	# Finds prodNum
                mzt   <- type.convert(ExpectMass[i])
                machine <- t_hdr[grep("^Type", t_hdr[,1]),2]
                ionmode <- t_hdr[grep("IonMode", t_hdr[,1]),2]
                bpi <- type.convert(t_hdr[grep("^BPI", t_hdr[,1]),2])
                mzl_1 <- paste("Target Mass", mzt)	# set Expected Mass to label 1
                mzl_2 <- paste(machine, ionmode,bpi)  # set MS mode to label 2
                ptime <- t_hdr[max(grep("^Time", t_hdr[,1])),2]
                ProcDesc <- t_hdr[grep("ProcDesc", t_hdr[,1]),2]
                MaxEnt <- str_extract(ProcDesc, "[[:alpha:]]{6} [[:digit:]]{1,2}")
                combine <- unlist(str_split(ProcDesc, "; "))
                mzl_3 <- paste("Time", ptime, MaxEnt,combine[2])  # set peak info to label 3
                # line type: h=histogram-like spikes for mass, o=connect-dots for hplc
                mzlabs  <- c(paste("Mass Spectrum:", dfb),"",mzl_1,mzl_2,mzl_3,"m/z","%","h")
                
                lcl_1 <- "UV Detector: 259.494-260.505" 	# hplc conditions
                buf_col <- t_hdr[grep("^Conditions", t_hdr[,1]),2]
                lcl_2 <- buf_col  # hplc conditions
                lcl_3 <- ""  # hplc conditions
                lclabs  <- c(paste("HPLC:",dfb),"",lcl_1,lcl_2,lcl_3,"Time","Absorbance %","o")
                
                fn_mz <- coa_plot(x_mzi,y_mzi,x_txt,y_txt,l_txt,TRUE,dfb,"mz",mzlabs, l_col="red")        # mass spec
                fn_lc <- coa_plot(x_lci,y_lci,x_lci_min,y_lci_min,y_lci_min_txt,TRUE,dfb,"lc",lclabs)	# hplc
                
                fn_list <- paste(fn_list,fn_mz,fn_lc,sep='|')
        }
}


#rm(list=ls())

