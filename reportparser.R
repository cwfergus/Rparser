##### Script Description ####
# The following script reads a MassLynx LC/MS report file
# It then parses out sample information, It grabs all the information shown in the ESMS import
# script final screen, except that which is not in the rpt file.
# It then compiles a data frame with each sample as a row
# It then writes out a tab delimited text file

#####Data Read in ####

# Install/load required R libraries
if (!require(stringr)){
        install.packages("stringr")
        library(stringr)
} else library(stringr)
if (!require(dplyr)){
        install.packages("dplyr")
        library(dplyr)
} else library(dplyr)

# Initialize some vectors and data frames for use later
data <- data.frame()

filename <- "y"
fn_list <- "You have successfully parsed the following reports:"
fn_full_list <- vector()

#This loop allows read in of multiple report files. It assumes the computer has 
#a connection to DATA4///current_report_files
#It also checks if an _reprocess exists, or if you have entered a barcode twice

while (filename != "n") {
        #Ask for Barcode, append .rpt and _reprocess.rpt
        filename <- readline("Please enter the QC set Barcode. 
                             If you are done adding reports enter n\t")
        filenamerpt <- paste(filename, ".rpt", sep="")
        filenamereprocess <- paste(filename, "_reprocess.rpt", sep="")
        #Set Data Path
        data_path <- "//DATA4/Mass Spec DataRrepository/reportfiles/current_report_files/"
        #Paste Data Path to filenames
        filenamefull <- paste(data_path, filenamerpt, sep="")
        reprocess_filename <- paste(data_path, filename, "_reprocess.rpt", sep="")
        
        #Check if a reprocess exists and prompt user if they want to parse that instead
        if (file.exists(reprocess_filename) == TRUE) {
                print("A more up-to-date version of this data may be availabe at:")
                print(filenamereprocess)
                decision <- readline("Do you want to import this instead? Enter y for yes\t")
                if (decision == "y") {
                        filenamefull <- reprocess_filename
                        filenamerpt <- filenamereprocess
                }
        }
        #Check if the user entered n, and end the read in loop
        breakname <- paste(data_path, "n.rpt", sep="")
        if (filenamefull == breakname) {
                break
        }
        #check if the data exists
        if (file.exists(filenamefull) == TRUE) {
                # Check if the user has already imported the data, if not read it in.
                if (filenamefull %in% fn_full_list == TRUE) {
                        print("You have already parsed this report")
                } else {
                        #This read in is sub-optimal, but the best I can find
                        #reads in data as a two unit list
                        tempdata <- scan(filenamefull, 
                                         what= list(Type="", Value = ""), #Reads into two columns as a list
                                         sep = "\t", #determines the two columns by the seperater tab
                                         multi.line = TRUE, #does not break 
                                         fill=TRUE,
                                         comment.char = "{",
                                         allowEscapes = FALSE)
                        #Converts the list to a data.frame
                        tempdata <- data.frame(do.call('cbind', tempdata),
                                               stringsAsFactors = FALSE)
                        #Append the tempdata to the growing data chain
                        data <- rbind(data, tempdata) 
                        #Add the filenames to check for repeat list
                        fn_full_list <- append(fn_full_list, filenamefull)
                        #Add Barcode.rpt to the print out list
                        fn_list <- append(fn_list, filenamerpt)
                }
                
        } else { #if the data doesn't exist, print out the full path. If path is wrong
                #User can now tell.
                print("No file exists with the following name!")
                print(filenamefull)
        }
}

##### Parser ####
# The Parser uses grep commands to look for certain words/Reg Exp. Then it takes a value
# in either the row over or somewhere else related to the word.



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
#Writes out the now ESMS replica data to the outputnume
write.table(likeESMS, "dataforimport.txt", sep = "\t", col.names=FALSE, quote=FALSE, row.names=FALSE, na="")
# print out the list of barcodes parsed and imported
print(fn_list)
print("Don't forget to now import into Filemaker")
#remove variables.
rm(list=ls())