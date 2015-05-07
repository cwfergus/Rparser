##### Script Description ####
# The following script reads a MassLynx LC/MS report file
# It then parses out sample information, It grabs all the information shown in the ESMS import
# script final screen, except that which is not in the rpt file.
# It then compiles a data frame with each sample as a row
# It then writes out a tab delimited text file

#####Data Read in ####
library(stringr)
#Asks user to specify the QC set barcode, ie the name of the .rpt file
data <- data.frame()

filename <- "y"

while (filename != "n") {
        filename <- readline("Please enter the QC set Barcode w/ .rpt. If you are done adding reports enter n \t")
        if (filename == "n") {
                break
        }
        tempdata <- scan(filename, 
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
#Searches in column Type for the word Well, and returns the value in the column next to it
#Each of these commands makes a vector with all of the values for each sample.
data$Value[grep("Well", data$Type)] -> Well #name chosen due to likeness with ESMS import
# Repeat the word well, for the number of samples, as determined by the pevious find

#Find the UserName row, take the info in the Value Column
data$Value[grep("UserName", data$Type)] -> UserName
#Find the Instrument row. The search uses a Regular expression for find the word Instrument,
#with nothing before it, and nothing after it.
data$Value[grep("^Instrument$", data$Type)] -> Instrument



#Synthesis Batch: found within FM, not possible in R

# SSID found within FM? I can probably do this here if necessary
#Searches for the SampleID row, finds the Production number, makes a vector with all of them
data$Value[grep("^SampleID$", data$Type)] -> Sample_Name

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

#Make a vector of column names just to make it easy.
columnnames <- c("Well", "Sample Name", "Expected Mass", "Results", "Area","Date", "Time", "User Name", "Instrument")
#Bind together the different vectors into a data frame. Here we make the MS Testtype DF
Outputdata <- data.frame(cbind(Well,
                           Sample_Name,
                           ExpectMass,
                           Results,
                           Area,
                           Date,
                           Time,
                           UserName,
                           Instrument))

#Append the new column names, removing the Area_null name to Area
colnames(Outputdata) <- columnnames

#Writes out the now ESMS replica data to the outputnume
write.table(Outputdata, "Parseddata.txt", sep = "\t", col.names=TRUE, quote=FALSE)
