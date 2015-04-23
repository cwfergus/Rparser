filename <- readline("Enter the QC set barcode, without the .rpt\t")
filenamerpt <- paste0(filename, ".rpt", sep="")

outputname <- paste0(filename, ".txt", sep ="")

scan(filenamerpt, 
     what= list(Type="", Value = ""), 
     sep = "\t",
     multi.line = TRUE,
     fill=TRUE,
     comment.char = "{",
     allowEscapes = FALSE) -> test

df <- data.frame(do.call('cbind', test),
           stringsAsFactors = FALSE)

df$Value[grep("Well", df$Type)] -> well

df$Value[grep("UserName", df$Type)] -> UserName

df$Value[grep("^Instrument$", df$Type)] -> Instrument

#Plate is always equal to one?

#Convert Well to Row, and well w/o SO shelf number, and format: A:1, B:1 etc

#Column is the second number in reformated well

#Synthesis Batch: found within FM, not possible in R

# SSID found within FM? I can probably do this here if necessary

df$Value[grep("^SampleID$", df$Type)] -> ProdNum


df$Value[grep("^Date$", df$Type)] -> Date
as.Date(Date, format="%d-%b-%Y") -> Date
format(Date, format = "%m/%d/%Y") -> Date


df$Value[grep("^Time$", df$Type)] -> Time
grep(":", Time, value=TRUE) -> Time
format(strptime(Time, format = "%X"), format = "%r") -> Time

df$Type[grep("Formula", df$Type)+1] -> ExpectMass

ExpectMass[3] <- 13067

df$Value[(grep("BPM", df$Type))] -> Results

samplesplit <- split(df, cumsum(df[,"Type"] == "[SAMPLE]"))
Area <- vector()
for (i in 1:length(samplesplit)) {
        workingdata <- as.data.frame(samplesplit[i])
        colnames(workingdata) <- c("Type", "Value")
        percentareas <- as.numeric(workingdata$Value[grep("Area %Total", workingdata$Type)])
        Area <- append(Area, max(percentareas))
}

#TestType, is this necessary? I think I may be able to come up with something...

df$Value[grep("SampleDescription", df$Type)] -> SO

#Customer: in filemaker not in input

#SeqID: in filemaker not in input

df$Value[(grep("JobCode", df$Type))]-> ESBarCode

likeESMS <- data.frame(cbind(well, UserName, Instrument, ProdNum, Date, Time, ExpectMass, Results, Area, SO, ESBarCode))

write.table(likeESMS, outputname, sep = "\t", col.names=FALSE, quote=FALSE)
