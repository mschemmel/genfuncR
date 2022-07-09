#' plots an arrow in an defined direction representing a single gene
#' @param x1 x1 coordinate
#' @param y1 y1 coordinate
#' @param x2 x2 coordinate
#' @param y2 y2 coordinate
#' @param direction Direction of the gene (upstream | downstream)
#' @param arr_type Appearance of arrow head (arrow | barrow | box)
#' @examples
#' genearrow(x1 = 1, y1 = 5, x2 = 4, y2 = 6, arrow_type = "arrow", direction = "downstream")

# TODO: add multiple arrow head types
genearrow <- function(x1, y1, x2, y2, direction, arr_type = "arrow", ...) {
    arrsize <- x2 - x1
    arrw <- arrsize * 0.2 # TODO: make size dynamic
    barrow_head_size <- 0.05
    if (direction == "downstream") {
        if (arr_type == "arrow") {
            grid.polygon(c(x1, x1, x2 - arrw, x2, x2 - arrw),
                         c(y1, y2, y2, mean(c(y1, y2)), y1), ...)
        }
        else if (arr_type == "barrow") {
            grid.polygon(c(x1, x1, x2 - arrw, x2 - arrw, x2, x2 - arrw, x2 - arrw),
                         c(y1, y2, y2, y2 + unit(barrow_head_size, "npc"), mean(c(y1, y2)), y1 - unit(barrow_head_size, "npc"), y1), ...)
        }
        else if (arr_type == "box") {
            grid.polygon(c(x1, x1, x2, x2),
                         c(y1, y2, y2, y1), ...)
        }
    }
    else if (direction == "upstream") {
        if (arr_type == "arrow") {
            grid.polygon(c(x1 + arrw, x1, x1 + arrw, x2, x2),
                         c(y1, mean(c(y2, y1)), y2, y2, y1), ...)
        }
        else if (arr_type == "barrow") {
            grid.polygon(c(x1 + arrw, x1 + arrw, x1, x1 + arrw, x1 + arrw, x2, x2),
                         c(y1, y1 - unit(barrow_head_size, "npc"), mean(c(y1, y2)), y2 + unit(barrow_head_size, "npc"), y2, y2, y1), ...)
        }
        else if (arr_type == "box") {
            grid.polygon(c(x1, x1, x2, x2),
                         c(y1, y2, y2, y1), ...)
        }
    }
    else {
        stop("Invalid 'direction' statement")
    }
}

#' puts a text label on a specific position
#' @param vp_name Name of the viewport
#' @param x_ x coordinate
#' @param y_ y coordinate
#' @param w_ w coordinate
#' @param h_ h coordinate
#' @param label_txt Text to be printed on that position
#' @param angle Angle of text label
#' @param gp_ Grid parameter like font, style, color, ...
#' @examples
#' text_label(x_ = 1, y_ = 3, w_ = 1, h_ = 1, "Test", angle = 45)

text_label <- function(vp_name = NULL, x_, y_, w_, h_, label_txt, angle = 0, gp_ = NULL) {
    grid::pushViewport(grid::viewport(name = vp_name,
                                      x = grid::unit(x_, "npc"),
                                      y = grid::unit(y_, "npc"),
                                      width = w_,
                                      height = h_,
                                      just = c("center")))
    grid::grid.text(label_txt, gp = gp_, rot = angle)
    grid::popViewport(1)
}

#' draw a set of genes based on data.frame or gff file
#' @param gff_file data.frame of gff file or pure own data.frame
#' @param forward_color color of genes in forward direction
#' @param reverse_color color of genes in reverse direction
#' @param arrow_type shape of arrow head
#' @examples
#' geneset(gff)

geneset <- function(gff_file,
                    forward_color = "darkslategray",
                    reverse_color = "navajowhite3",
                    arrow_type = "arrow") {

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
    min_value <- min(gff_file$start, gff_file$end)
    s1_pos <- 0.8
    s2_pos <- 0.2
    genomic_vp_width_x0 <- 0
    genomic_vp_width_x1 <- 1

    # helper function
    relative <- function(x, max_x) ((x * genomic_vp_width_x1) / max_x)

    # genomic viewport
    grid::pushViewport(grid::viewport(name = "genomic",
                        x = grid::unit(0.5, "npc"),
                        y = grid::unit(0.5, "npc"),
                        width = 0.9,
                        height = 1,
                        just = c("center")))

    # forward and reverse direction
    grid::grid.segments(x0 = unit(genomic_vp_width_x0, "npc"),
                        y0 = unit(s1_pos, "npc"),
                        x1 = unit(genomic_vp_width_x1, "npc"),
                        y1 = unit(s1_pos, "npc"))
    # reverse direction
    grid::grid.segments(x0 = unit(genomic_vp_width_x0, "npc"),
                        y0 = unit(s2_pos, "npc"),
                        x1 = unit(genomic_vp_width_x1, "npc"),
                        y1 = unit(s2_pos, "npc"))

    # TODO: add proper axis plus label
    grid::grid.xaxis(label = round(seq(min_value, max_value, (max_value - min_value) / 10), 0), at = seq(0, 1, 0.1))
    text_label(x_ = 0.5, y_ = -0.5, w_ = 0.1, h_ = 0.2, label_txt = "Region (bp)")

    # add features of gff (or dataframe) file
    for (i in seq_len(nrow(gff_file))) {
        # downstream
        if (gff_file$strand[i] == "+") {
            genearrow(x1 = unit(relative(gff_file$start[i], max_value), "npc"),
                      y1 = unit(0.7, "npc"),
                      x2 = unit(relative(gff_file$end[i], max_value), "npc"),
                      y2 = unit(0.9, "npc"),
                      direction = "downstream",
                      gp = gpar(fill = forward_color),
                      arr_type = arrow_type)
        }

        # upstream
        else if (gff_file$strand[i] == "-") {
            genearrow(x1 = unit(relative(gff_file$start[i], max_value), "npc"),
                      y1 = unit(0.1, "npc"),
                      x2 = unit(relative(gff_file$end[i], max_value), "npc"),
                      y2 = unit(0.3, "npc"),
                      direction = "upstream",
                      gp = gpar(fill = reverse_color),
                      arr_type = arrow_type)
        }
        else {
           warning("Unrecognized 'strand' symbol.")
        }
    }
    grid::popViewport(1)
}