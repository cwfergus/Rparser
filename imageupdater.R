data_path <- "//DATA4/Mass Spec DataRrepository/reportfiles/current_report_files/"
data <- data.frame()
file.list_reprocess <- list.files(data_path, pattern="[[:digit:]]{10}_reprocess.rpt$", full.names = TRUE)
file.list <- file.list_reprocess[456:length(file.list_reprocess)]
#file.list_reprocess <- file.list_reprocess[which(file.info(file.list_reprocess)$ctime >= "2015-06-01 00:00:01")]

file.list <- list.files(data_path, pattern="[[:digit:]]{10}.rpt$", full.names = TRUE)

file.list <- file.list[3692:length(file.list)]

file.list <- append(file.list, file.list_reprocess)

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

custom.labels <- function (og_x, og_y, og_labels = NULL, x_offsets = NA, ofs_amt=15) 
{
        
        
        rounder <- function(x) {round(x+10^-9)}
        
        posfinder <- function(x) {
                tens <- seq(0, 100, 10)
                
                for (i in 1:length(x)) {
                        value <- x[i]
                        x[i]<- tens[which.min(abs(tens-value))]
                        tens <- tens[which(tens!=x[i])]
                }
                x
        }
        
        limit_finder <- function(og_y) {x <- seq(0, 1, .05)
                                        
                                        for (i in 1:length(x)) {
                                                threshold <- x[i]
                                                min <- threshold*max(og_y)
                                                y_min <- which(og_y>min)
                                                if (length(y_min)<=10) {
                                                        threshold
                                                        break
                                                }
                                        }
                                        threshold
        }
        
        ny <- length(og_y)
        if (ny == 1) {
                text(og_x, og_y, og_labels, pos=3)
        } else {
                
                limit <- limit_finder(og_y)
                min <- limit*max(og_y)
                y_min <- which(og_y>min)
                og_x <- og_x[y_min]
                og_y <- og_y[y_min]
                og_labels <- og_labels[y_min]
                ny <- length(og_y)
                divide <- rounder(ny/2)
                
                if (is.na(x_offsets)) {
                        x_offset <- diff(par("usr")[1:2])/ofs_amt
                        x_offsets <- rep(c(rep(-x_offset, divide), rep(x_offset, divide)), ny/2 + 1)[1:ny]
                }
                
                new_x <- og_x + x_offsets
                
                x_left_pts <- new_x[1:divide]
                left_labels <- og_labels[1:divide]
                
                x_right_pts <- new_x[(divide+1):length(new_x)]
                right_labels <- og_labels[(divide+1):length(new_x)]
                
                
                
                y_left_pts <- og_y[1:divide]
                y_right_pts <- og_y[(divide+1):length(og_y)]
                
                y_left_pts <- posfinder(y_left_pts)
                y_right_pts <- posfinder(y_right_pts)
                
                new_y <- c(y_left_pts, y_right_pts)
                
                segments(new_x, new_y, og_x, og_y)
                text(x_left_pts, y_left_pts, left_labels, pos=2)
                text(x_right_pts, y_right_pts, right_labels, pos=4)
        }
        
}




coa_plot <- function(x,y,x_txt,y_txt,l_txt,pt_b,pt_c,fb,gr_t,labels,l_color="black", so)
{
        # x,y are vectors; pt_b TRUE/FALSe to show point txt; pt_min: min value to show pt txt
        # fb = file base name; gr_t = mz or lc plot type; labels includes titles, axis labels, and other specs
        
        "//Data/it/DBMS/Integrations/ERPDev/imgs/lcms/" -> md
        dev_fn <- paste(md, fb, "_", gr_t, ".png", sep="")
        
        png(filename=dev_fn, width=1400, height=500) # Open a device: bmp & pdf are other options here
        if (pt_b) {ylim <- c(0,105)} 
        else {
                ylim <- range(y)
                ylim[2] <- ylim[2]*1.05
        }
        
        plot(x, y, 
             main=labels[1], sub=labels[2],xlab=labels[6], ylab=labels[7], type="n", axes=TRUE,
             col.lab="blue", col.axis="blue", ylim = ylim
        )
        lines(x,y, lty=1, type=labels[8], pch='', col=l_color)
        if (pt_b) {
                
                custom.labels(x_txt, y_txt, l_txt, ofs_amt=35)
                text(min(x),max(ylim), labels[3], pos=4, offset=-1)
                mtext(paste("SO:", so), side=3, adj=0, line=2)
                mtext(labels[4], side=3, line=1, adj=1)
                mtext(labels[5], side=3, line=0, adj=1)
        }
        
        if(pt_c) {
                #textxy(x_txt, 0.96*y_txt, l_txt, cex=1, offset=-0.8)
                text(min(x), max(ylim), labels=l_txt, pos=4, offset=-1)
                mtext (labels[9], side=3, adj=0)
                mtext(paste("SO:", so), side=3, adj=0, line=2)
                mtext(labels[3], side=3, line=1, adj=1)
                mtext(labels[4], side=3, line=0, adj=1)
        }
        
        dev.off()
        
        
        path <- paste("//DATA4/Mass Spec DataRrepository/imgs", so, sep="/")
        if (!file.exists(path)) {
                dir.create(path)
        }
        fn <- paste(path, "/", fb,"_",gr_t,".png",sep="")                 # filename
        png(filename=fn, width=680, height=440)        # Open a device: bmp & pdf are other options here
        if (pt_b) {
                ylim <- c(0,105)
        } else {
                ylim <- range(y)
                ylim[2] <- ylim[2]*1.05
        }
        plot(x, y, 
             main=labels[1], sub=labels[2],xlab=labels[6], ylab=labels[7], type="n", axes=TRUE,
             col.lab="blue", col.axis="blue", ylim = ylim
        )
        lines(x,y, lty=1, type=labels[8], pch='', col=l_color)
        if (pt_b) {
                
                custom.labels(x_txt, y_txt, l_txt, ofs_amt=25)
                text(min(x),max(ylim), labels[3], pos=4, offset=-1)
                mtext(paste("SO:", so), side=3, adj=0, line=2)
                mtext(labels[4], side=3, line=0.7, adj=1, cex=0.8)
                mtext(labels[5], side=3, line=0, adj=1, cex=0.8)
        }
        if(pt_c) {
                text(min(x), max(ylim), labels=l_txt, pos=4, offset=-1)
                mtext(paste("SO:", so), side=3, adj=0, line=2)
                mtext (labels[9], side=3, adj=0, cex=0.8)
                mtext(labels[3], side=3, line=0.7, adj=1, cex=0.8)
                mtext(labels[4], side=3, line=0, adj=1, cex=0.8)
        }
        
        dev.off()
        # Close the image device that was opened above; the actions in the meantime have been recorded
        return(fn)
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
                s_int <- type.convert(df[peak_loc[l]+6,2])
                e_int <- type.convert(df[peak_loc[l]+7,1])
                table <- cbind(pl, peak_id, rt, area, peak_s_time, peak_e_time, peak_width, s_int, e_int)
                areatable <- rbind(areatable, table)
                
        }
        return(areatable)
}

for (i in 1:length(file.list)){
        
        filenamefull <- file.list[i]
        
        tempdata <- scan(filenamefull, 
                         what= list(Type="", Value = ""), #Reads into two columns as a list
                         sep = "\t", #determines the two columns by the seperater tab
                         multi.line = TRUE, #does not break 
                         fill=TRUE,
                         comment.char = "{",
                         allowEscapes = FALSE)
        #Converts the list to a data.frame
        data <- data.frame(do.call('cbind', tempdata),
                               stringsAsFactors = FALSE)
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
        
        #Generates empty vectors for the results to go into
        Results <- vector()
        Area <- vector()
        FlrArea <- vector()
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
                        if (any(grepl("% BPI", functionsplit[f]))){ #Is this MS data?
                                ms_data <- as.data.frame(functionsplit[f]) #put in data frame
                                colnames(ms_data) <- c("Type", "Value")
                                next
                        }
                        if (any(grepl("DAD", functionsplit[f]))){ #Is it LC data?
                                dad_data <- as.data.frame(functionsplit[f])
                                colnames(dad_data) <- c("Type", "Value")
                                next
                        }
                        if (any(grepl("FLR", functionsplit[f]))){ #Is it FLR data?
                                flr_data2 <- as.data.frame(functionsplit[f])
                                colnames(flr_data2) <- c("Type", "Value")
                                next
                        }
                }
                
                ## Find Mass value:
                if (exists("ms_data")) { #If Ms data exists, find the mass
                        Foundmass <- ms_data$Value[grep("BPM", ms_data$Type)]
                } else { #Otherwise, put 0 for the mass and report missing data
                        Foundmass <- 0
                        msg <- paste(ProdNum[ss], "is Missing MS data")
                        print(msg)
                } #append the sample result to all results
                Results <- append(Results, Foundmass) 
                
                ## Find LC Area
                if (exists("dad_data")) { #If LC data exists, find Max Area %
                        percentareas <- as.numeric(dad_data$Value[grep("Area %Total", dad_data$Type)])
                        Foundarea <- max(percentareas)
                } else { #Otherwise put 0 for Area and report missing data
                        Foundarea <- 0
                        msg <- paste(ProdNum[ss], "is missing LC data")
                        print(msg)
                } #Append result to all results
                Area <- append(Area, Foundarea)
                
                ## Find FLR Area
                if (exists("flr_data2")) { #If FLR data exits, find FLR area
                        percentareas2 <- as.numeric(flr_data2$Value[grep("Area %Total", flr_data2$Type)])
                        FoundFLRarea <- max(percentareas2)
                } else {#Otherwise put 0 for FLR area
                        FoundFLRarea <- 0
                        #Don't bother reporting as most samples won't have FLR data.
                } #Append FLR area of sample to overall FlrArea vector
                FlrArea <- append(FlrArea, FoundFLRarea)
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
        columnnames <- c("Well", "ProdNum", "Date", "Time", "UserName", "Instrument", 
                         "Results", "PlateLoc", "Row", "Col", "SSID", "ExpectMass", 
                         "TestType", "SO", "ESBarCode", "Area")
        #Bind together the different vectors into a data frame.
        msESMS <- data.frame(cbind(Well, ProdNum, Date, Time, UserName, Instrument, 
                                   Results, PlateLoc, Row, Col, SSID, ExpectMass, 
                                   TestType_MS, SO, ESBarCode, Area_null))
        
        areaESMS <- data.frame(cbind(Well, ProdNum, Date, Time, UserName, Instrument,
                                     Results_null, PlateLoc, Row, Col, SSID,
                                     ExpectMass_null, TestType_HPLC, SO, ESBarCode, Area))
        
        #Append the new column names, removing the Area_null name to Area
        colnames(msESMS) <- columnnames
        colnames(areaESMS) <- columnnames
        #puts the two TestType data frames together
        rbind(msESMS, areaESMS) -> likeESMS
        #Writes out the now ESMS replica data to the outputnume
        
        
        #####Graphing functions####
       
        
        
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
                        i_txt <- localMaxima(y_mzi)
                        x_LM <- x_mzi[i_txt]
                        y_LM <- y_mzi[i_txt]
                        y_min <- 0.25*max(y_mzi)
                        i_min <- which(y_LM>y_min)
                        x_txt <- x_LM[i_min]
                        y_txt <- y_LM[i_min]
                        if (length(y_mzi)<20){
                                centroid <- "yes"
                                x_txt <- x_mzi
                                y_txt <- y_mzi
                        }
                        #IF an FLR trace exists do the following:
                        if (length(g_lci) > 1) { #If an FLR trace exists, do the following
                                
                                #Find Range values for FLR and LC
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
                                #Convert Y values to AU units using the range.
                                y_lci_LC_au <- (y_lci_LC/100)*lc_range
                                #Generate labels for HPLC data
                                lc_data <- txti[g_lci[1]:(g_lci[2]-1),]
                                lc_areatable <- peak_table(lc_data)
                                peak_lc <- lc_areatable[which(lc_areatable$area==max(lc_areatable$area)),]
                                x_peak_lc <- peak_lc$rt
                                y_peak_lc <- peak_lc$area
                                lc_peak <- peak_lc$pl
                                y_peak_txt_lc <- paste("Area % Purity:", y_peak_lc, "%")
                                
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
                                mzlabs  <- c(paste("LGC Biosearch Technologies", "\n", "MS:", dfb, "\n"),"",mzl_1,mzl_2,mzl_3,"m/z","%","h")
                                m_y_txt <- y_txt[which(y_txt>0.25*max(y_txt))]
                                x_txt <- x_txt[which(y_txt==m_y_txt)]
                                y_txt <- m_y_txt
                                rm(centroid)
                        } else {
                                mzlabs  <- c(paste("LGC Biosearch Technologies", "\n", "MS:", dfb, "\n"),"",mzl_1,mzl_2,mzl_3,"m/z","%","l")
                        }
                        
                        
                        lcl_1 <- grep("DAD:", txti[,2], value=TRUE) 	# hplc conditions
                        lcl_1 <- paste("Diode Array Detector:", str_sub(lcl_1, start=6))
                        buf_col <- t_hdr[grep("^Conditions", t_hdr[,1]),2]
                        lcl_2 <- buf_col  # hplc conditions
                        lcl_3 <- ""  # hplc conditions
                        lclabs  <- c(paste("LGC Biosearch Technologies", "\n", "UPLC:",dfb, "\n"),"",lcl_1,lcl_2,lcl_3,"Time (min)","AU","o", lc_dis_range)
                        
                        
                        fn_mz <- coa_plot(x_mzi,y_mzi,x_txt,y_txt,x_txt,TRUE,FALSE,dfb,"mz",mzlabs, l_col="red", so=SO[i])        # mass spec
                        
                        if (length(g_lci) == 1) {
                                fn_lc <- coa_plot(x_lci,y_lci_au,x_peak,lc_range,y_peak_txt,FALSE,TRUE,dfb,"lc",lclabs, so=SO[i])	# hplc
                        } else {
                                fn_lc <- coa_plot(x_lci_LC, y_lci_LC_au, x_peak_lc, lc_range, y_peak_txt_lc,FALSE, TRUE, dfb, "lc", lclabs, so=SO[i])
                                
                                flrl_1 <- grep("ACQUITY FLR", txti[,2], value=TRUE)
                                flrl_2 <- buf_col
                                flrl_3 <- ""
                                flrlabs <- c(paste("FLR", dfb),"",flrl_1, flrl_2, flrl_3, "Time", "Emission %", "o", flr_dis_range)
                                
                                fn_flr <- coa_plot(x_lci_flr, y_lci_flr_au, x_peak_flr, flr_range, y_peak_txt_flr,FALSE,TRUE,dfb,"fl",flrlabs, so=SO[i])
                        }
                        
                }
        }
        
        
        
}