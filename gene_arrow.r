# random gff file
gff <- data.frame(seqname = c("Chr1", "Chr1", "Chr2", "Chr2", "Chr2", "Chr3", "Chr3"),
                    source = rep("Genome",7),
                    feature = rep("Gene",7),
                    start = c(34, 370, 800, 1100, 1500, 2020, 2500),
                    end = c(364,700, 950, 1250, 2000, 2200, 2700),
                    score = rep(".", 7),
                    strand = c("+","+","-","+","+","-","+"),
                    frame = c(rep(0,7)),
                    attribute = paste0("seq", seq(1,7,1))
                    )

#' plots an arrow in an defined direction representing a single gene

#' @param x1 x1 coordinate
#' @param y1 y1 coordinate
#' @param x2 x2 coordinate
#' @param y2 y2 coordinate
#' @param direction Direction of the gene (upstream | downstream)
#' return grid polygon object

#' @examples
#' genearrow(1,5,4,6,"downstream")

genearrow <- function(x1, y1, x2, y2, direction, ...){
    arrsize <- x2 - x1
    arrw <- arrsize * 0.2 # TODO: make size dynamic

    if(direction == "downstream"){
         grid.polygon(c(x1, x1, x2 - arrw, x2, x2 - arrw),
                      c(y1, y2, y2, mean(c(y1, y2)), y1), ...)
    }
    else if(direction == "upstream"){
        grid.polygon(c(x1 + arrw, x1, x1 + arrw, x2, x2),
                     c(y1, mean(c(y2, y1)), y2, y2, y1), ...)
    }
    else {
        stop("Invalid 'direction' statement")
    }
}

geneset <- function(gff_file) {
    # start plotting grids
    # create newpage to draw on
    grid::grid.newpage()

    # outer viewport
    grid::pushViewport(grid::viewport(name = "outer",
                        x = grid::unit(0.5, "npc"),
                        y = grid::unit(0.5, "npc"),
                        width = 1,
                        height = 1))

    # main viewport
    grid::pushViewport(grid::viewport(name = "main",
                        x = grid::unit(0.5, "npc"),
                        y = grid::unit(0.5, "npc"),
                        width = 0.7,
                        height = 0.3,
                        just = c("center")))
                        


    # store some constants
    max_value <- max(gff_file$start, gff_file$end)
    s1_height <- 0.8
    s2_height <- 0.2
    genomic_vp_width_x0 <- 0
    genomic_vp_width_x1 <- 1

    # helper function
    relative <- function(x, max_x) ((x*genomic_vp_width_x1)/max_x)
    
    # genomic viewport
    grid::pushViewport(grid::viewport(name = "genomic",
                        x = grid::unit(0.5, "npc"),
                        y = grid::unit(0.5, "npc"),
                        width = 0.9,
                        height = 1,
                        just = c("center")))

    # forward and reverse direction
    grid::grid.segments(x0 = unit(genomic_vp_width_x0, "npc"), y0 = unit(s1_height, "npc"), x1 = unit(genomic_vp_width_x1, "npc"), y1 = unit(s1_height, "npc"))   
    grid::grid.segments(x0 = unit(genomic_vp_width_x0, "npc"), y0 = unit(s2_height, "npc"), x1 = unit(genomic_vp_width_x1, "npc"), y1 = unit(s2_height, "npc"))   
    #grid::grid.xaxis(label = seq(0,max_value, max_value/10), at = seq(0, 1, 0.1))

    # add genes
    for(i in seq(1:nrow(gff_file))) {
        # downstream
        if(gff_file$strand[i] == "+") {
            genearrow(x1 = unit(relative(gff_file$start[i], max_value), "npc"),
                      y1 = unit(0.7, "npc"),
                      x2 = unit(relative(gff_file$end[i], max_value), "npc"),
                      y2 = unit(0.9, "npc"), direction = "downstream", gp = gpar(fill = "darkslategray"))
        }

        # upstream
        else if (gff_file$strand[i] == "-") {
            genearrow(x1 = unit(relative(gff_file$start[i], max_value), "npc"),
                      y1 = unit(0.1, "npc"),
                      x2 = unit(relative(gff_file$end[i], max_value), "npc"),
                      y2 = unit(0.3, "npc"), direction = "upstream", gp = gpar(fill = "navajowhite3"))
        }
        else {
           warning("Start and end position are the same.")
        }
    }
    grid::popViewport(1)
}



