#all report read in

filenames <- list.files("Reportfiles", pattern="*.rpt", full.names=TRUE)
# fix this again...
filestoreadin <- filenames[!filenames %in% readinfiles]

if (length(filestoreadin) == 0){
        print("No new report files to read in")
}

data <- ldply(lapply(filestoreadin, scan, 
       what= list(Type="", Value = ""),
       sep = "\t",
       multi.line = TRUE,
       fill=TRUE,
       comment.char = "{",
       allowEscapes = FALSE), as.data.frame, stringsAsFactors=FALSE) 


readinfiles <- filenames
