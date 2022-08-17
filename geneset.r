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
#' @param gff_file data.frame of gff file or pure own data.frame
#' @param cov_file data.frame of genomice coverage
#' @param forward_color color of genes in forward direction (default = "darkslategray")
#' @param reverse_color color of genes in reverse direction (default = "darkslategray")
#' @param transparency alpha value of gene annotation box
#' @param arrow_type shape of arrow head (default = "arrow")
#' @param gene_height height in percent (0-1) of gene box (default = 1)
#' @param distance distance between forward and revers strand in percent (0-1) (default = 1)
#' @param show_axis show axis or not (default = TRUE)
#' @param axis_label_text text of x axis label (default: "Region (bp)")
#' @param axis_interval numerical interval of axis (default = NULL)
#' @param range vector of start and end of region of interest (default = NULL)
#' @param border boolean if border visible (default = FALSE)
#' @param show_values boolean if displayed range should also be printed
#' @param annotation vector of locations where to draw line annotations
#' @examples
#' geneset(gff)

geneset <- function(gff_file,
                    cov_file = NULL,
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
                    show_values = FALSE,
                    annotation = NULL) {

    # check if input file is valid
    if (nrow(gff_file) == 0) {
        stop("Input data frame is empty")
    }
    if (!all(c("start", "end", "strand") %in% colnames(gff_file))) {
        cat("Input has to contain at least 'start', 'end' and 'strand' column.\n")
        cat("Found colnames: ", colnames(gff_file), "\n")
        stop()
    }

    # order input by start and end column and filter
    gff_file <- gff_file[with(gff_file, order(gff_file$start, gff_file$end)), ]
    gff_file <- gff_file[gff_file$start > range[1] & gff_file$end < range[2], ]

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
                        y = grid::unit(0.3, "npc"),
                        width = 0.7,
                        height = 0.3,
                        just = c("center")))

    # draw border if requested
    if (border) grid::grid.rect()

    # store some constants
    min_value <- range[1]
    max_value <- range[2]

    gene_box_height <- ifelse(gene_height > 0 & gene_height <= 1,
                              0.2 * gene_height,
                              stop("Height of gene box (gene_height) have to be between 0 and 1."))
    s2_pos <- 0.2
    gap <- ifelse(distance >= 0 & distance <= 1,
                  gene_box_height + ((s2_pos + gene_box_height) * distance),
                  stop("Distance between forward and reverse strand have to be between 0 and 1."))

    s1_pos <- s2_pos + gap
    genomic_vp_width_x0 <- 0
    genomic_vp_width_x1 <- 1

    # helper function
    relative <- function(x) (x - min_value) / (max_value - min_value)

    # forward direction
    grid::grid.segments(x0 = grid::unit(genomic_vp_width_x0, "npc"),
                        y0 = grid::unit(s1_pos, "npc"),
                        x1 = grid::unit(genomic_vp_width_x1, "npc"),
                        y1 = grid::unit(s1_pos, "npc"))
    # reverse direction
    grid::grid.segments(x0 = grid::unit(genomic_vp_width_x0, "npc"),
                        y0 = grid::unit(s2_pos, "npc"),
                        x1 = grid::unit(genomic_vp_width_x1, "npc"),
                        y1 = grid::unit(s2_pos, "npc"))


    # add features of gff
    for (i in seq_len(nrow(gff_file))) {
        genearrow(x1 = grid::unit(relative(gff_file$start[i]), "npc"),
                  y1 = ifelse(gff_file$strand[i] == "+", grid::unit(s1_pos - (gene_box_height / 2), "npc"),
                                                         grid::unit(s2_pos - (gene_box_height / 2), "npc")),
                  x2 = grid::unit(relative(gff_file$end[i]), "npc"),
                  y2 = ifelse(gff_file$strand[i] == "+", grid::unit(s1_pos + (gene_box_height / 2), "npc"),
                                                         grid::unit(s2_pos + (gene_box_height / 2), "npc")),
                  direction = ifelse(gff_file$strand[i] == "+", "downstream", "upstream"),
                  gp_ = grid::gpar(fill = ifelse(gff_file$strand[i] == "+", forward_color, reverse_color),
                                   alpha = transparency),
                  arr_type = arrow_type)
    }

    # add axis label
    if (show_axis) {
        # set axis label
        axis_label <- round(seq(min_value, max_value, axis_interval), 0)

        grid::grid.xaxis(label = axis_label,
                         at = seq(0, 1, 1 / (length(axis_label) - 1)))

        # add axis label text
        text_label(x_ = 0.5,
                   y_ = -0.3,
                   w_ = 0.1,
                   h_ = 0.2,
                   label_txt = axis_label_text)
    }

    # add line annotations
    for (annot in annotation) {
        grid::grid.segments(x0 = grid::unit(relative(annot), "npc"),
                            y0 = grid::unit(s1_pos + (gene_box_height / 2), "npc"),
                            x1 = grid::unit(relative(annot), "npc"),
                            y1 = grid::unit(s2_pos - (gene_box_height / 2), "npc"),
                            gp = grid::gpar(fill = "red", col = "red", lwd = 1))
    }

    grid::popViewport(1)

    if (!is.null(cov_file)) {
        # coverage viewport
        grid::pushViewport(grid::viewport(name = "main",
                            x = grid::unit(0.5, "npc"),
                            y = grid::unit(0.65, "npc"),
                            width = 0.7,
                            height = 0.3,
                            just = c("center")))
        grid::grid.rect()

        # add coverage of gff
        cov_file <- cov_file[cov_file$start >= range[1] & cov_file$start <= range[2], ]
        max_in_range <- max(cov_file$cov)

        for (i in seq_len(nrow(cov_file))) {
            if(cov_file$cov[i] != 0) {
                grid::grid.segments(x0 = grid::unit(relative(cov_file$start[i]), "npc"),
                                    y0 = grid::unit(0, "npc"),
                                    x1 = grid::unit(relative(cov_file$start[i]), "npc"),
                                    y1 = grid::unit((cov_file$cov[i] / max_in_range), "npc"),
                                    gp = grid::gpar(fill = "gold3", col = "gold3", lwd = 0.1))
            }    
        }
        grid::popViewport(1)
    }
}