
custom.labels <- function (sl_x, sl_y, sl_labels = NULL, x_offsets = NA, offset_amt=15
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
        
        ny <- length(sl_y)
        divide <- rounder(ny/2)
        
        if (is.na(x_offsets)) {
                x_offset <- diff(par("usr")[1:2])/offset_amt
                x_offsets <- rep(c(rep(-x_offset, divide), rep(x_offset, divide)), ny/2 + 1)[1:ny]
        }
        
        new_x <- sl_x + x_offsets
        
        x_left_pts <- new_x[1:divide]
        left_labels <- sl_labels[1:divide]
        
        x_right_pts <- new_x[(divide+1):length(new_x)]
        right_labels <- sl_labels[(divide+1):length(new_x)]
        
        
        
        y_left_pts <- sl_y[1:divide]
        y_right_pts <- sl_y[(divide+1):length(sl_y)]
        
        y_left_pts <- posfinder(y_left_pts)
        y_right_pts <- posfinder(y_right_pts)
        
        new_y <- c(y_left_pts, y_right_pts)
        #                         sl_x <- sl_x + x_offsets
        #                         sort.index <- sort.list(sl_x)
        #                         sl_x <- sl_x[sort.index]
        #                         sl_y <- sl_y[sort.index]
        #                         nx <- length(sl_x)
        #                         newx <- seq(sl_x[1], sl_x[nx], length = length(sl_labels))
        
        segments(new_x, new_y, sl_x, sl_y)
        
        text(x_left_pts, y_left_pts, left_labels, srt = srt, pos=2)
        text(x_right_pts, y_right_pts, right_labels, srt=srt, pos=4)
        
}
