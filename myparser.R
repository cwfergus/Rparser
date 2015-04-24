##### Script Description ####
# The following script reads a MassLynx LC/MS report file
# It then parses out sample information, It grabs all the information shown in the ESMS import
# script final screen, except that which is not in the rpt file.
# It then compiles a data frame with each sample as a row
# It then writes out a tab delimited text file

#####Data Read in ####
library(stringr)
#Asks user to specify the QC set barcode, ie the name of the .rpt file
filename <- readline("Enter the QC set barcode, without the .rpt\t")
#Ads the .rpt to the end 
filenamerpt <- paste0(filename, ".rpt", sep="")
#Makes a write out name that is based off the .rpt name
outputname <- paste0(filename, ".txt", sep ="")
#reads in the data from the filenamerpt
scan(filenamerpt, 
     what= list(Type="", Value = ""), #Reads into two columns as a list
     sep = "\t", #determines the two columns by the seperater tab
     multi.line = TRUE, #does not break 
     fill=TRUE,
     comment.char = "{",
     allowEscapes = FALSE) -> data
#Binds the list into a dataframe, and doesn't change characters to factors
data <- data.frame(do.call('cbind', data),
           stringsAsFactors = FALSE)

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
                workingdata3 <- as.data.frame(chromosplit[2])
                colnames(workingdata3) <- c("Type", "Value")
                percentareas3 <- as.numeric(workingdata3$Value[grep("Area %Total", workingdata3$Type)])
                Area <- append(Area, max(percentareas3))
        } else {
                percentareas <- as.numeric(workingdata$Value[grep("Area %Total", workingdata$Type)])
                Area <- append(Area, max(percentareas)) }
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

#Makes vectors with empty fields for each of the two TestType splits
ExpectMass_null <- vector(mode = "character", length = length(ExpectMass))
Results_null <- vector(mode = "character", length = length(Results))
Area_null <- vector(mode = "character", length = length(Area))

#Make a vector of column names just to make it easy.
columnnames <- c("f1", "f2", "UserName", "Instrument", "PlateLoc", "Row", "Well", "Col", "ProdNum", "Date", "Time", "ExpectMass", "Results", "Area", "TestType", "SO", "ESBarCode")
#Bind together the different vectors into a data frame. Here we make the MS Testtype DF
msESMS <- data.frame(cbind(f1,
                           f2,
                           UserName,
                           Instrument,
                           PlateLoc,
                           Row,
                           Well,
                           Col,
                           ProdNum,
                           Date,
                           Time,
                           ExpectMass,
                           Results,
                           Area_null,
                           TestType_MS,
                           SO,
                           ESBarCode))

#Append the new column names, removing the Area_null name to Area
colnames(msESMS) <- columnnames
#Bind together the HPLC test type stuff
areaESMS <- data.frame(cbind(f1,
                             f2,
                             UserName,
                             Instrument,
                             PlateLoc,
                             Row,
                             Well,
                             Col,
                             ProdNum,
                             Date,
                             Time,
                             ExpectMass_null,
                             Results_null,
                             Area,
                             TestType_HPLC,
                             SO,
                             ESBarCode))
#Remove the ExpectMass_null, and Results_null
colnames(areaESMS) <- columnnames
#puts the two TestType data frames together
rbind(msESMS, areaESMS) -> likeESMS
#Writes out the now ESMS replica data to the outputnume
write.table(likeESMS, "Importdata.txt", sep = "\t", col.names=FALSE, quote=FALSE, row.names=FALSE)
