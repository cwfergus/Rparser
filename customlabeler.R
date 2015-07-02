
custom.labels <- function (og_x, og_y, og_labels = NULL, x_offsets = NA, ofs_amt=15,
                           linecol = par("fg"), srt = 0, ...) 
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
                text(x_left_pts, y_left_pts, left_labels, srt = srt, pos=2)
                text(x_right_pts, y_right_pts, right_labels, srt=srt, pos=4)
        }
        
}
