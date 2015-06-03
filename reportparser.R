##### Script Description ####
# The following script reads a MassLynx LC/MS report file
# It then parses out sample information, It grabs all the information shown in the ESMS import
# script final screen, except that which is not in the rpt file.
# It then compiles a data frame with each sample as a row
# It then writes out a tab delimited text file

#####Data Read in ####
if (!require(stringr)){
        install.packages("stringr")
        library(stringr)
} else library(stringr)
if (!require(dplyr)){
        install.packages("dplyr")
        library(dplyr)
} else library(dplyr)
if (!require(calibrate)){
        install.packages("calibrate")
        library(calibrate)
} else library(calibrate)
if (!require(plotrix)){
        install.packages("plotrix")
        library(plotrix)
} else library(plotrix)


# Initialize some vectors and data frames for use later
data <- data.frame()

filename <- "y"
fn_list <- "You have successfully parsed the following reports:"
fn_full_list <- vector()

#This loop allows read in of multiple report files. It assumes the computer has 
#a connection to DATA4///current_report_files
#It also checuks if an _reprocess exists, or if you have entered a barcode twice

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

#####Graphing functions####
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
        
          "//Data/it/DBMS/Integrations/ERPDev/imgs/lcms/" -> md
          dev_fn <- paste(md, fb, "_", gr_t, ".png", sep="")
          png(filename=dev_fn, width=1400, height=500) # Open a device: bmp & pdf are other options here
        if (pt_b) {ylim <- c(0,105)} else {ylim <- range(y)}
          
          plot(x, y, 
             main=labels[1], sub=labels[2],xlab=labels[6], ylab=labels[7], type="n", axes=TRUE,
             col.lab="blue", col.axis="blue", ylim = ylim
        )
        if (pt_b) {
                if (length(l_txt)>=10) {
                        min <- 0.5*max(y_txt)
                        y_min <- which(y_txt>min)
                        x_txt <- x_txt[y_min]
                        y_txt <- y_txt[y_min]
                        l_txt <- l_txt[y_min]
                        points(x_txt, y_txt, type="p", pch=c(15,16,17), col=rainbow(n=10, start=2/6), cex=2)
                        legend("right", l_txt, pch=c(15,16,17), col=rainbow(n=10, start=2/6))
                } else {
                        if (length(l_txt)==2){
                                x_txt <- append(0, x_txt)
                                y_txt <- append(0, y_txt)
                                l_txt <- append("", l_txt)
                        }
                        thigmophobe.labels(x_txt, y_txt, l_txt, cex = 1)} 
        }
        if(pt_c) {
                textxy(x_txt, 0.96*y_txt, l_txt, cex=1, offset=-0.8)
                polygon(shd, density=NA, col="green", border="black")
                text (min(x), max(y), labels[9], pos=4)
        }
        lines(x,y, lty=1, type=labels[8], pch='', col=l_color)
        
        text(max(x),max(y), labels[3], pos=2)
        text(max(x),0.96*max(y), labels[4], pos=2)
        text(max(x),0.92*max(y), labels[5], pos=2)
          dev.off()
        
        
        path <- paste("//DATA4/Mass Spec DataRrepository/imgs", so, sep="/")
        if (!file.exists(path)) {
                dir.create(path)
        }
        fn <- paste(path, "/", fb,"_",gr_t,".png",sep="")                 # filename
        png(filename=fn, width=680, height=440)	# Open a device: bmp & pdf are other options here
        plot(x, y, 
             main=labels[1], sub=labels[2],xlab=labels[6], ylab=labels[7], type="n", axes=TRUE,
             col.lab="blue", col.axis="blue"
        )
        lines(x,y, lty=1, type=labels[8], pch='', col=l_color)
        if (pt_b) {
                if (length(l_txt)>=10) {
                        min <- 0.5*max(y_txt)
                        y_min <- which(y_txt>min)
                        x_txt <- x_txt[y_min]
                        y_txt <- y_txt[y_min]
                        l_txt <- l_txt[y_min]
                        points(x_txt, y_txt, type="p", pch=c(15,16,17), col=rainbow(n=10, start=2/6), cex=2)
                        legend("right", l_txt, pch=c(15,16,17), col=rainbow(n=10, start=2/6))
                } else {
                        if (length(l_txt)==2){
                                x_txt <- append(0, x_txt)
                                y_txt <- append(0, y_txt)
                                l_txt <- append("", l_txt)
                        }
                        thigmophobe.labels(x_txt, y_txt, l_txt, cex = 1)} 
        }
        if(pt_b) {text(max(x),max(y), labels[3], pos=2)
                  mtext(labels[4], side=3, line=1, adj=1)
                  mtext(labels[5], side=3, line=0, adj=1)
        }
        if(pt_c) {
                textxy(x_txt, 0.96*y_txt, l_txt, cex=1, offset=-0.8)
                polygon(shd, density=NA, col="green", border="black")
                mtext (labels[9], side=3, adj=0)
                mtext(labels[3], side=3, line=1, adj=1)
                mtext(labels[4], side=3, line=0, adj=1)
        }
        
        
       
        dev.off()
        # Close the image device that was opened above; the actions in the meantime have been recorded
        return(fn)
}

curvecolor <- function(p_table, tdf, df, range) {
        
        selc_peak <- p_table[which(p_table$area==max(p_table$area)),]
        
        
        
        shd_bl <- grep(min(selc_peak$peak_s_time), tdf[,1], fixed=TRUE)
        shd_el <- grep(max(selc_peak$peak_e_time), tdf[,1], fixed=TRUE)
        
        if (length(shd_bl) == 0) {
                bt_m1 <- min(selc_peak$peak_s_time)-0.0001
                shd_bl_m1 <- grep(bt_m1, tdf[,1], fixed=TRUE)
                bt_p1 <- min(selc_peak$peak_s_time)+0.0001
                shd_bl_p1 <- grep(bt_p1, tdf[,1], fixed=TRUE)
                if (length(shd_bl_m1) > 0) shd_bl <- shd_bl_m1
                if (length(shd_bl_p1) > 0) shd_bl <- shd_bl_p1
        }
        if (length(shd_el) == 0) {
                et_m1 <- max(selc_peak$peak_e_time)-0.0001
                shd_el_m1 <- grep(et_m1, tdf[,1], fixed=TRUE)
                et_p1 <- max(selc_peak$peak_e_time)+0.0001
                shd_el_p1 <- grep(et_p1, tdf[,1], fixed=TRUE)
                if (length(shd_el_m1) > 0) shd_el <- shd_el_m1
                if (length(shd_el_p1) > 0) shd_el <- shd_el_p1
        }
        
        if (length(shd_bl)==0|length(shd_el)==0) {
                shd <- paste("Unable to generate shader for:", ProdNum[i])
                return (shd)
        } else {        shd_b <- min(shd_bl)
                        shd_e <- max(shd_el)
                        shd <- type.convert(tdf[shd_b:shd_e,])
                        shd[,2] <- (shd[,2]/100)*range
                        shd_b_au <- selc_peak$s_int[1]
                        shd_e_au <- selc_peak$e_int[length(selc_peak$pl)]
                        shd[1,2] <- min(shd_b_au, shd[1,2])
                        shd[length(shd[,1]),2] <- min(shd_e_au, shd[length(shd[,1]),2])
                        
                        return (shd)}       

                
}

peak_table <- function(df) {
        peak_loc <- grep("[PEAK]", df[,1], fixed=TRUE)
        areatable <- data.frame()
        for (l in 1:length(peak_loc)) {
                pl <- peak_loc[l]
                rt <- type.convert(df[peak_loc[l]+3,2])
                area <- type.convert(df[peak_loc[l]+11,2])
                peak_id <- type.convert(df[peak_loc[l]+1,2])
                peak_s_time <- type.convert(df[peak_loc[l]+4,2])
                peak_e_time <- type.convert(df[peak_loc[l]+5,1])
                peak_width <- type.convert(df[peak_loc[l]+12,2])
                s_int <- type.convert(df[peak_loc[l]+6,2])/1000000
                e_int <- type.convert(df[peak_loc[l]+7,1])/1000000
                table <- cbind(pl, peak_id, rt, area, peak_s_time, peak_e_time, peak_width, s_int, e_int)
                areatable <- rbind(areatable, table)
                
        }
        return(areatable)
}


#####Graphing Parser####

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
        g_spi <- grep('}', txti[,1], fixed=TRUE)	# vector of all end-bracket locations
        g_pks <- grep("[PEAK]", txti[,1], fixed=TRUE)	# vector of all PEAK locations
        if (length(g_lci) == 0 | length(g_mzi) ==0)  { 
                err_msg <- paste (ProdNum[i], "is missing LC or MS data")
                print(err_msg) }
        else { 
                # analyze mass spec data
                m_mzi <- g_spi[min(which(g_spi>g_mzi))]-1# get the lowest position of } that's after [MS]
                t_mzi <- txti[(g_mzi+2):m_mzi,]		# get all data between [MS]{ and }
                x_mzi <- type.convert(t_mzi[,1])   	# Type is x; make numeric
                y_mzi <- type.convert(t_mzi[,2])     	# Value is y; make it numeric
                # calculate text labels for mass spec; only the points that will be printed
                y_min <- 0.25*max(y_mzi)     # must be over minimum 25% threshold
                i_min <- which(y_mzi>y_min)
                X_min <- x_mzi[i_min]
                y_min <- y_mzi[i_min]
                if (length(y_min)== 1) {
                        centroid <- "yes"
                        x_txt <- x_mzi
                        y_txt <- y_mzi
                } else {i_txt <- localMaxima(y_min) 
                        x_txt <- X_min[i_txt]
                        y_txt <- y_min[i_txt]}
                 # only peak maxima instead of clusters of points, doesn't need change.
                
                #IF an FLR trace exists do the following:
                if (length(g_lci) > 1) {
                        
                        lcflr_range <- type.convert(txti[grep("MaxIntensity", txti[,1]),2])
                        lc_range <- lcflr_range[1]/1000000
                        flr_range <- lcflr_range[2]
                        lc_dis_range <- paste("Range:", round(lc_range, digits=3))
                        flr_dis_range <- paste("Range:", flr_range)
                        #Analyze the HPLC data
                        m_lci_LC <- g_spi[min(which(g_spi>g_lci[1]))]-1
                        t_lci_LC <- txti[(g_lci[1]+1):m_lci_LC,]
                        x_lci_LC <- type.convert(t_lci_LC[,1])
                        y_lci_LC <- type.convert(t_lci_LC[,2])
                        
                        y_lci_LC_au <- (y_lci_LC/100)*lc_range
                        #Generate labels for HPLC data
                        lc_data <- txti[g_lci[1]:(g_lci[2]-1),]
                        lc_areatable <- peak_table(lc_data)
                        peak_lc <- lc_areatable[which(lc_areatable$area==max(lc_areatable$area)),]
                        x_peak_lc <- peak_lc$rt
                        y_peak_lc <- peak_lc$area
                        lc_peak <- peak_lc$pl
                        y_peak_txt_lc <- paste("Area % Purity:", y_peak_lc, "%")
                        
                        lc_shd <- curvecolor(lc_areatable, t_lci_LC, lc_data, lc_range)
                        
                        
                        #Analyze FLR data
                        m_lci_flr <- g_spi[min(which(g_spi>g_lci[2]))]-1
                        t_lci_flr <- txti[(g_lci[2]+1):m_lci_flr,]
                        x_lci_flr <- type.convert(t_lci_flr[,1])
                        y_lci_flr <- type.convert(t_lci_flr[,2])
                        
                        y_lci_flr_au <- (y_lci_flr/100)*flr_range
                        #Generate labels for FLR data
                        flr_data <- txti[g_lci[2]:max(g_spi),]
                        flr_areatable <- peak_table(flr_data)
                        peak_flr <- flr_areatable[which(flr_areatable$area==max(flr_areatable$area)),]
                        x_peak_flr <- peak_flr$rt
                        y_peak_flr <- peak_flr$area
                        flr_peak <- peak_flr$pl
                        y_peak_txt_flr <- paste("Fluoresent % Purity:", y_peak_flr, "%")
                        
                        flr_shd <- curvecolor(flr_areatable, t_lci_flr, flr_data, flr_range)
                        
                        #Find Range Values
                
                
                
                } else {
                # analyze hplc data
                        lc_range <- (type.convert(txti[grep("MaxIntensity", txti[,1]),2]))/1000000
                        lc_dis_range <- paste("Range:", round(lc_range, digits=3))
                        
                m_lci <- g_spi[min(which(g_spi>g_lci))]-1	# get the lowest position of } that's after [TRACE]
                t_lci <- txti[(g_lci+1):m_lci,]			# get all data between [TRACE]{ and }
                x_lci <- type.convert(t_lci[,1])     	# every other point is x; make it numeric
                y_lci <- type.convert(t_lci[,2]) 
                y_lci_au <- (y_lci/100)*lc_range
                # every other point is y; make it numeric
                # calculate text labels for LC
                areatable <- peak_table(txti)
                peak <- areatable[which(areatable$area==max(areatable$area)),]
                x_peak <- peak$rt
                y_peak <- peak$area
                peakloc <- peak$pl
                y_peak_txt <- paste("Area % Purity:", y_peak, "%")
                
                #Generate Peak integration graphing data
                lc_shd <- curvecolor(areatable, t_lci, txti, lc_range)
                #Find LC range value
                
                }
                
                
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
                mzl_3 <- paste("Time", ptime, MaxEnt, grep("Combine", combine, value=TRUE))  # set peak info to label 3
                # line type: h=histogram-like spikes for mass, o=connect-dots for hplc
                if (exists("centroid")){
                        mzlabs  <- c(paste("Mass Spectrum:", dfb),"",mzl_1,mzl_2,mzl_3,"m/z","%","h")
                        m_y_txt <- y_txt[which(y_txt>0.25*max(y_txt))]
                        x_txt <- x_txt[which(y_txt==m_y_txt)]
                        y_txt <- m_y_txt
                        rm(centroid)
                } else {
                        mzlabs  <- c(paste("Mass Spectrum:", dfb),"",mzl_1,mzl_2,mzl_3,"m/z","%","l")
                }
                
                
                lcl_1 <- grep("DAD:", txti[,2], value=TRUE) 	# hplc conditions
                buf_col <- t_hdr[grep("^Conditions", t_hdr[,1]),2]
                lcl_2 <- buf_col  # hplc conditions
                lcl_3 <- ""  # hplc conditions
                lclabs  <- c(paste("HPLC:",dfb),"",lcl_1,lcl_2,lcl_3,"Time","AU","o", lc_dis_range)
                
                
                fn_mz <- coa_plot(x_mzi,y_mzi,x_txt,y_txt,x_txt,TRUE,FALSE,dfb,"mz",mzlabs, l_col="red", so=SO[i])        # mass spec
                if (class(lc_shd)=="character") {
                        print(lc_shd)
                } else {
                if (length(g_lci) == 1) {
                fn_lc <- coa_plot(x_lci,y_lci_au,x_peak,lc_range,y_peak_txt,FALSE,TRUE,dfb,"lc",lclabs, shd=lc_shd, so=SO[i])	# hplc
                } else {
                        fn_lc <- coa_plot(x_lci_LC, y_lci_LC_au, x_peak_lc, lc_range, y_peak_txt_lc,FALSE, TRUE, dfb, "lc", lclabs, shd=lc_shd, so=SO[i])
                        
                        flrl_1 <- grep("ACQUITY FLR", txti[,2], value=TRUE)
                        flrl_2 <- buf_col
                        flrl_3 <- ""
                        flrlabs <- c(paste("FLR", dfb),"",flrl_1, flrl_2, flrl_3, "Time", "Emission %", "o", flr_dis_range)
                        
                        fn_flr <- coa_plot(x_lci_flr, y_lci_flr_au, x_peak_flr, flr_range, y_peak_txt_flr,FALSE,TRUE,dfb,"fl",flrlabs, shd=flr_shd, so=SO[i])
                }
                }
                
        }
}


#rm(list=ls())

##### testing####

# infl_find <- function(values, loc) {
#         infl <- c(FALSE, diff(diff(values)>0)!=0)
#         point <- which(infl==TRUE)
#         if (length(point)==0 & loc =="bl") {
#                 point <- 1
#         }
#         if (length(point)==0 & loc =="el") {
#                 point <- length(values)
#         } 
#         if (length(point)>1) {
#                 point <- point[1]
#         }
#         return(point)
# }
# 
# curvecolor <- function(ploc, tdf, df, range) {
#         
#         shd_bl <- grep(df[ploc+4, 2], tdf[,1], fixed=TRUE)
#         #         shd_bl <- ((head(shd_bl, n=1)-40): (tail(shd_bl, n=1)))
#         shd_bl_v <- type.convert(tdf[shd_bl,2])
#         
#         shd_el <- grep(df[ploc+5,1], tdf[,1], fixed=TRUE)
#         #         shd_el <- ((head(shd_el, n=1)): (tail(shd_el, n=1)+10))
#         shd_el_v <- type.convert(tdf[shd_el,2])
#         
#         shd_b <- shd_bl[infl_find(shd_bl_v, "bl")]
#         shd_e <- shd_el[infl_find(shd_el_v, "el")]
#         
#         
#         
#         #         if (length(shd_bl) == 0) {
#         #                 bt_m1 <- type.convert(df[ploc+4,2])-0.0001
#         #                 shd_bl_m1 <- grep(bt_m1, tdf[,1], fixed=TRUE)
#         #                 bt_p1 <- type.convert(df[ploc+4,2])+0.0001
#         #                 shd_bl_p1 <- grep(bt_p1, tdf[,1], fixed=TRUE)
#         #                 if (length(shd_bl_m1) > 0) shd_bl <- shd_bl_m1
#         #                 if (length(shd_bl_p1) > 0) shd_bl <- shd_bl_p1
#         #         }
#         #         if (length(shd_el) == 0) {
#         #                 et_m1 <- type.convert(df[ploc+5,1])-0.0001
#         #                 shd_el_m1 <- grep(et_m1, tdf[,1], fixed=TRUE)
#         #                 et_p1 <- type.convert(df[ploc+5,1])+0.0001
#         #                 shd_el_p1 <- grep(et_p1, tdf[,1], fixed=TRUE)
#         #                 if (length(shd_el_m1) > 0) shd_el <- shd_el_m1
#         #                 if (length(shd_el_p1) > 0) shd_el <- shd_el_p1
#         #         }
#         
#         if (length(shd_bl)==0|length(shd_el)==0) {
#                 shd <- paste("Unable to generate shader for:", ProdNum[i])
#                 return (shd)
#         } else {
#                 shd <- type.convert(tdf[shd_b:shd_e,])
#                 (shd[,2]*range)/100 -> shd[,2]
#                 shd_b_au <- (type.convert(df[ploc+6,2])/1000000)
#                 shd_e_au <- (type.convert(df[ploc+7,1])/1000000)
#                 
#                 if(shd_b_au > shd[1,2]) {}
#                 
#                 shd[1,2] <- min((type.convert(df[ploc+6,2])/1000000), shd[1,2])
#                 shd[length(shd[,1]),2] <- min((type.convert(df[ploc+7,1])/1000000), shd[length(shd[,1]),2])
#                 
#                 return (shd)}
#         
# }
