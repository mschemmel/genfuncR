genearrow <- function(x1, y1, x2, y2, direction, ...){
    arrsize = x2 - x1
    scalefactor <- ifelse(arrsize < 0.4, 0.4, 0.05)
    arrw = arrsize * scalefactor
    
    if(direction == "downstream"){
         garrow <- grid.polygon(c(x1, x1, x2 - arrw, x2, x2 - arrw),
                        c(y1, y2, y2, mean(c(y1, y2)), y1), ...)
    }
    else if(direction == "upstream"){
        garrow <- grid.polygon(c(x1 + arrw, x1, x1 + arrw, x2, x2),
                        c(y1, mean(c(y2, y1)), y2, y2, y1), ...)
    }
    else {
        stop("Invalid 'direction' statement")
    }
    return(garrow)
}

