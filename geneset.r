#' plots an arrow in an defined direction representing a single gene
#' @param x1 x1 coordinate
#' @param y1 y1 coordinate
#' @param x2 x2 coordinate
#' @param y2 y2 coordinate
#' @param direction Direction of the gene (upstream | downstream)
#' @param arr_type Appearance of arrow head (arrow | barrow | box)
#' @examples
#' genearrow(x1 = 1, y1 = 5, x2 = 4, y2 = 6, arrow_type = "arrow", direction = "downstream")

genearrow <- function(x1, y1, x2, y2, direction, arr_type = "arrow", gp_) {
    arrsize <- x2 - x1
    arrw <- arrsize * 0.2 # TODO: make size dynamic
    barrow_head_size <- 0.05
    if (direction == "downstream") {
        if (arr_type == "arrow") {
            grid::grid.polygon(c(x1, x1, x2 - arrw, x2, x2 - arrw),
                               c(y1, y2, y2, mean(c(y1, y2)), y1),
                               gp = gp_)
        } else if (arr_type == "barrow") {
            grid::grid.polygon(c(x1, x1, x2 - arrw, x2 - arrw, x2, x2 - arrw, x2 - arrw),
                               c(y1, y2, y2, y2 + barrow_head_size, mean(c(y1, y2)), y1 - barrow_head_size, y1),
                               gp = gp_)
        } else if (arr_type == "box") {
            grid::grid.polygon(c(x1, x1, x2, x2),
                               c(y1, y2, y2, y1),
                               gp = gp_)
        } else {
           stop("Invalid 'arrow_type' parameter. Please choose 'arrow', 'barrow' or 'box'.")
        }
    } else if (direction == "upstream") {
        if (arr_type == "arrow") {
            grid::grid.polygon(c(x1 + arrw, x1, x1 + arrw, x2, x2),
                               c(y1, mean(c(y2, y1)), y2, y2, y1),
                               gp = gp_)
        } else if (arr_type == "barrow") {
            grid::grid.polygon(c(x1 + arrw, x1 + arrw, x1, x1 + arrw, x1 + arrw, x2, x2),
                               c(y1, y1 - barrow_head_size, mean(c(y1, y2)), y2 + barrow_head_size, y2, y2, y1),
                               gp = gp_)
        } else if (arr_type == "box") {
            grid::grid.polygon(c(x1, x1, x2, x2),
                               c(y1, y2, y2, y1),
                               gp = gp_)
        } else {
           stop("Invalid 'arrow_type' parameter. Please choose 'arrow', 'barrow' or 'box'.")
        }
    } else {
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
#' @param gff_file data.frame of gff file
#' @param forward_color color of genes in forward direction (default = "darkslategray")
#' @param reverse_color color of genes in reverse direction (default = "darkslategray")
#' @param transparency alpha value of gene annotation box
#' @param arrow_type shape of arrow head (default = "arrow")
#' @param gene_height height in percent (0-1) of gene box (default = 1)
#' @param distance distance between forward and reverse strand in percent (0-1) (default = 1)
#' @param show_axis show axis or not (default = TRUE)
#' @param axis_label_text text of x axis label (default: "Region (bp)")
#' @param axis_interval numerical interval of axis (default = NULL)
#' @param range vector of start and end of region of interest (default = NULL)
#' @param border boolean if border visible (default = FALSE)
#' @param show_values boolean if displayed range should also be printed
#' @examples
#' geneset(gff)

geneset = setClass("geneset",
                   slots = list(
                        gff_file = "ANY",
                        gene_param = "list",
                        plot_param = "list"
                   ))

geneset <- function(gff_file = NULL,
                    forward_color = "darkslategray",
                    reverse_color = "darkslategray",
                    transparency = 1,
                    arrow_type = "arrow",
                    gene_height = 1,
                    distance = 1,
                    show_axis = TRUE,
                    axis_label_text = "Region (bp)",
                    axis_interval = NULL,
                    range = NULL,
                    border = FALSE,
                    show_values = FALSE) {

    .Object = new("geneset")

    # check if input file is valid
    if (!is.null(gff_file)) {
        if (nrow(gff_file) == 0) {
            stop("Input data frame is empty")
        }
        if (!all(c("chr", "start", "end", "strand") %in% colnames(gff_file))) {
            cat("Input has to contain at least 'start', 'end' and 'strand' column.\n")
            cat("Found colnames: ", colnames(gff_file), "\n")
            stop()
        }
    }

    # order input by start and end column and filter
    gff_file <- gff_file[with(gff_file, order(gff_file$start, gff_file$end)), ]
    gff_file <- gff_file[gff_file$start > range[1] & gff_file$end < range[2], ]
    .Object@gff_file = gff_file

    if(show_values) print(gff_file[which(gff_file$start > range[1] & gff_file$end < range[2]), ])

    # store some constants
    min_value <- range[1]
    max_value <- range[2]

    # TODO: solve axis intervals properly
    # check for proper axis interval value
    interval <- max_value - min_value
    if (interval %% axis_interval != 0) {
        stop("Please provide a multiple of your specified region for proper axis style.")
    }
    if (interval / axis_interval > 100) {
        stop("Provided axis interval will result in more than 100 labels.")
    }

    # overall size of genebox
    gene_box_height <- ifelse(gene_height > 0 & gene_height <= 1,
                              0.2 * gene_height,
                              stop("Height of gene box (gene_height) have to be between 0 and 1."))
    .Object@gene_param$gene_box_height = gene_box_height
    
    # position of forward and reverse strand
    forward_strand_pos <- 0.2
    strand_gap <- ifelse(distance >= 0 & distance <= 1,
                                     gene_box_height + ((forward_strand_pos + gene_box_height) * distance),
                                     stop("Distance between forward and reverse strand have to be between 0 and 1."))
    reverse_strand_pos <- forward_strand_pos + strand_gap

    .Object@gene_param$forward_strand_pos = forward_strand_pos
    .Object@gene_param$strand_gap = strand_gap
    .Object@gene_param$reverse_strand_pos = reverse_strand_pos
    .Object@gene_param$distance = distance

    # define default params
    .Object@plot_param$min_value = min_value
    .Object@plot_param$max_value = max_value
    .Object@gene_param$forward_color = forward_color
    .Object@gene_param$reverse_color = reverse_color
    .Object@gene_param$transparency = transparency
    .Object@gene_param$arrow_type = arrow_type
    .Object@gene_param$gene_height = gene_height
    .Object@plot_param$interval = interval
    .Object@plot_param$show_axis = show_axis
    .Object@plot_param$axis_label_text = axis_label_text
    .Object@plot_param$axis_interval = axis_interval
    .Object@plot_param$border = border
    .Object@plot_param$show_values = show_values

    return(.Object)
}

setMethod(f = "show",
          signature = "geneset",
          definition = function(object) {
            # helper function
            relative <- function(x) (x - object@plot_param$min_value) / (object@plot_param$max_value - object@plot_param$min_value)
            # create new device and newpage
            dev.new(width = 12, height = 6, unit = "in")
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

            # draw border if requested
            if (object@plot_param$border) grid::grid.rect()


            # forward direction
            grid::grid.segments(x0 = grid::unit(0, "npc"),
                                y0 = grid::unit(object@gene_param$reverse_strand_pos, "npc"),
                                x1 = grid::unit(1, "npc"),
                                y1 = grid::unit(object@gene_param$reverse_strand_pos, "npc"))
            # reverse direction
            grid::grid.segments(x0 = grid::unit(0, "npc"),
                                y0 = grid::unit(object@gene_param$forward_strand_pos, "npc"),
                                x1 = grid::unit(1, "npc"),
                                y1 = grid::unit(object@gene_param$forward_strand_pos, "npc"))

            # add features of gff
            for (i in seq_len(nrow(object@gff_file))) {
                genearrow(x1 = grid::unit(relative(object@gff_file$start[i]), "npc"),
                          y1 = ifelse(object@gff_file$strand[i] == "+",
                                      grid::unit(object@gene_param$reverse_strand_pos - (object@gene_param$gene_box_height / 2), "npc"),
                                      grid::unit(object@gene_param$forward_strand_pos - (object@gene_param$gene_box_height / 2), "npc")),
                          x2 = grid::unit(relative(object@gff_file$end[i]), "npc"),
                          y2 = ifelse(object@gff_file$strand[i] == "+",
                                      grid::unit(object@gene_param$reverse_strand_pos + (object@gene_param$gene_box_height / 2), "npc"),
                                      grid::unit(object@gene_param$forward_strand_pos + (object@gene_param$gene_box_height / 2), "npc")),
                          direction = ifelse(object@gff_file$strand[i] == "+",
                                             "downstream",
                                             "upstream"),
                          gp_ = grid::gpar(fill = ifelse(object@gff_file$strand[i] == "+",
                                                         object@gene_param$forward_color,
                                                         object@gene_param$reverse_color),
                                          alpha = object@gene_param$transparency),
                          arr_type = object@gene_param$arrow_type)
            }

            # add axis label
            if (object@plot_param$show_axis) {
                # set axis label
                axis_label <- round(seq(object@plot_param$min_value, object@plot_param$max_value, object@plot_param$axis_interval), 0)

                grid::grid.xaxis(label = axis_label,
                                 at = seq(0, 1, 1 / (length(axis_label) - 1)))

                # add axis label text
                text_label(x_ = 0.5,
                           y_ = -0.3,
                           w_ = 0.1,
                           h_ = 0.2,
                           label_txt = object@plot_param$axis_label_text)
            }
            grid::popViewport(1)
            })