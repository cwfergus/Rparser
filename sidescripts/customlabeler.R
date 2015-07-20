
og_x <- x_txt; og_y <- y_txt; og_labels <- x_txt;  ofs_amt <- 15; linecol = par("fg"); srt <- 0
x <- x_mzi; y <- y_mzi; x_txt <- x_txt; y_txt <- y_txt; l_txt <- x_txt; labels <- mzlabs; fb <- dfb; gr_t <- "mz" ; l_color <- "red"; so <- SO[i]


custom.labels <- function (og_x, og_y, og_labels = NULL, ofs_amt=15,
                           linecol = par("fg"), srt = 0, ...) 
{
        
        
        rounder <- function(x) {round(x+10^-9)}
        
        posfinder <- function(data) {
                l_boxY <- ceiling(l_w)+(l_w/8)
                
                tens <- seq(0, 100, l_boxY)
                name_list <- vector("character")
                for (i in unique(data$y_group)){
                        name <- paste("tens", i, sep="_")
                        assign(name, tens)
                        name_list <- append(name_list, name)
                }
                for ( i in 1:nrow(data)){
                        if (data$pos[i]==3){
                                data$order[i] <- 1
                        } else {
                                data$order[i] <- 2
                        }
                }
                
                data<-arrange(data, order, desc(og_y))
                for (i in 1:nrow(data)){
                        group <- data$y_group[i]
                        ten_grp_name <- paste("tens", group, sep="_")
                        temp_tens <- get(ten_grp_name)
                        value <- data$og_y[i]
                        temp_pos <- data$pos[i]
                        if (temp_pos==3){
                                data$new_y[i] <- value
                                y <- temp_tens[which.min(abs(temp_tens-(value+3)))]
                                temp_tens <- temp_tens[which(temp_tens!=y)]
                        } else {
                                data$new_y[i]<- temp_tens[which.min(abs(temp_tens-value))]
                                temp_tens <- temp_tens[which(temp_tens!=data$new_y[i])]
                        }
                        assign(ten_grp_name, temp_tens)
                }
                data
        }
        
        limit_finder <- function(og_y) {x <- seq(0, 1, .02)
                                        
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
                #Determine which points to label, labels largest 10.
                limit <- limit_finder(og_y)
                min <- limit*max(og_y)
                y_min <- which(og_y>min)
                og_x <- og_x[y_min]
                og_y <- og_y[y_min]
                og_labels <- og_labels[y_min]
                
                #create a label table to facilitate "smart" labeling.
                x_group <- vector(mode="numeric", length=length(og_x))
                y_group <- vector(mode="numeric", length=length(og_x))
                near_left <- vector(mode="logical", length=length(og_x))
                near_right <- vector(mode="logical", length=length(og_x))
                new_x <- vector(mode="numeric", length=length(og_x))
                new_y <- vector(mode="numeric", length=length(og_x))
                pos <- vector(mode="numeric", length=length(og_x))
                order <- vector(mode="numeric", length=length(og_x))
                label_table <- data.frame(og_x, og_y, og_labels, x_group, y_group, near_left, near_right, new_x, new_y, pos, order)
                
                plot_boundries <- par("usr")[1:2]
                x_offset <- diff(par("usr")[1:2])/ofs_amt
                
                l_l<- max(strwidth(og_labels))
                l_w <- max(strheight(og_labels))
                
                margin_s <- x_offset+l_l + 5
                
                
                #group values by X value
                diffs <- diff(og_x)
                grp_borders <- which(diffs>(l_l/2))
                
                if (length(grp_borders)==0){
                        label_table$x_group <- NA
                } else {
                        for (i in 1:length(grp_borders)){
                                label_table[grp_borders[i], "x_group"] <- i
                        }
                }
                
                
                label_table[which(label_table$x_group==0),"x_group"] <- NA
                label_table$x_group <- na.locf(label_table$x_group, na.rm=FALSE, fromLast = TRUE)
                
                label_table[is.na(label_table$x_group), "x_group"] <- (length(grp_borders)+1)
                
                #Determine if an original X value in a group is near plot boundries
                for (i in unique(label_table$x_group)){
                        x_group <- i
                        smallest <- min(label_table$og_x[which(label_table$x_group==x_group)])
                        largest <- max(label_table$og_x[which(label_table$x_group==x_group)])
                        if (smallest-plot_boundries[1] <= margin_s){
                                label_table$near_left[which(label_table$x_group==x_group)] <- TRUE
                        }
                        if (plot_boundries[2]-largest <= margin_s){
                                label_table$near_right[which(label_table$x_group==x_group)] <- TRUE
                        }
                }
                
                
                #Determine new X values based on grouping/plot boundries
                for (i in unique(label_table$x_group)){
                        x_group <- i
                        tempdata <- label_table[which(label_table$x_group==x_group),]
                        
                        ny <- nrow(tempdata)
                        if (tempdata$near_left[1]) {
                                
                                x_offsets <- rep(x_offset, ny)[1:ny]
                                new_x <- tempdata$og_x + x_offsets
                                #new_y <- posfinder(tempdata$og_y)
                                tempdata$pos <- 4
                        }
                        if (tempdata$near_right[1]){
                                
                                x_offsets <- rep(-x_offset, ny)[1:ny]
                                new_x <- tempdata$og_x + x_offsets
                                #new_y <- posfinder(tempdata$og_y)
                                tempdata$pos <- 2
                        }
                        if (!exists("x_offsets")){
                                if(nrow(tempdata)==1){
                                        tempdata$pos <- 3
                                        new_x <- tempdata$og_x
                                        new_y <- tempdata$og_y
                                        x_offsets <- NULL
                                } else {
                                        divide <- rounder(ny/2)
                                        
                                        x_offsets <- rep(c(rep(-x_offset, divide), rep(x_offset, divide)), ny/2 + 1)[1:ny]
                                        new_x <- tempdata$og_x + x_offsets
                                        
                                        y_left_pts <- tempdata$og_y[1:divide]
                                        tempdata$pos[1:divide] <- rep(2, length(y_left_pts))
                                        y_right_pts <- tempdata$og_y[(divide+1):length(tempdata$og_y)]
                                        tempdata$pos[(divide+1):length(tempdata$og_y)] <- rep(4, length(y_right_pts))
                                        
                                }
                        }
                        label_table$new_x[which(label_table$x_group==x_group)] <- new_x
                        #                         label_table$new_y[which(label_table$x_group==x_group)] <- new_y
                        label_table$pos[which(label_table$x_group==x_group)] <- tempdata$pos
                        rm(x_offsets)
                }
                
                #Determine Y groups by new X value locations
                diffs <- diff(label_table$new_x)
                grp_borders <- which(diffs>(2*margin_s))
                
                #Need to determine this much better. Need Pos to take a roll (see 1000719168, i=3)
        
                
                if (length(grp_borders)==0){
                        label_table$y_group <- NA
                } else {
                        for (i in 1:length(grp_borders)){
                                label_table[grp_borders[i], "y_group"] <- i
                        }
                }
                label_table[which(label_table$y_group==0),"y_group"] <- NA
                label_table$y_group <- na.locf(label_table$y_group, na.rm=FALSE, fromLast = TRUE)
                
                label_table[is.na(label_table$y_group), "y_group"] <- (length(grp_borders)+1)
                
                #Determine new Y values based on y groups.
                
                label_table <- posfinder(label_table)
                
                
                for (i in 1:nrow(label_table)) {
                        label <- label_table[i,]
                        segments(label$new_x, label$new_y, label$og_x, label$og_y)
                        text(label$new_x, label$new_y, label$og_label, srt=srt, pos=label$pos)
                }
                
        }
        
}
