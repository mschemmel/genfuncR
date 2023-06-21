#' genearrow plots an arrow in an defined direction representing a single gene
#' @param x1 x1 coordinate
#' @param x2 x2 coordinate
#' @param pos positional coordinates of forward and reverse strand
#' @param direction direction of the gene (upstream | downstream)
#' @param gd_ grid gpar object to edit the arrow appearance
#' @examples
#' drawGene(x1 = 1, x2 = 4, pos = 0.4, direction = "forward")
drawGene <- function(x1, x2, pos, direction, forward_color, reverse_color, gene_height = 0.15) {
    arrow_head_width <- (x2 - x1) * 0.2 # TODO: make size dynamic
    y1 <- pos - gene_height
    y2 <- pos + gene_height

    ifelse(direction == "+",
           grid::grid.polygon(c(x1, x1, x2 - arrow_head_width, x2, x2 - arrow_head_width),
                              c(y1, y2, y2, pos, y1),
                              gp = grid::gpar(fill = forward_color)),
           grid::grid.polygon(c(x1 + arrow_head_width, x1, x1 + arrow_head_width, x2, x2),
                              c(y1, pos, y2, y2, y1),
                              gp = grid::gpar(fill = reverse_color))
    )
}

#' draw strand (5' ------ 3') on specified position
#' @param direction character string specifying forward or reverse strand
#' @param y_ y position to draw
#' @param show_direction boolean if direction labels should be displayed
#' @examples
#' drawStrand(direction = "forward", y_ = 0.4)
drawStrand <- function(direction = "+", y_, show_direction = TRUE) {
  direction_label <- if (direction == "+") c("5'", "3'") else c("3'", "5'")

  grid::grid.segments(x0 = grid::unit(0, "npc"),
                      y0 = grid::unit(y_, "npc"),
                      x1 = grid::unit(1, "npc"),
                      y1 = grid::unit(y_, "npc"))
  if (identical(show_direction, TRUE)) {
    grid::grid.text(x = c(-0.05, 1.05), y = y_, label = direction_label, gp = grid::gpar(fontsize = 16))
  }
}

#' draw Feature on specified position
#' @param pos position to put
#' @param strand_pos y position to draw
#' @param value numerical value of feature
#' @param col color of feature
drawFeature <- function(pos, strand_pos, value, col, type = "circle") {
  if(type == "circle") {
    grid::grid.circle(x = grid::unit(pos, "npc"),
                      y = grid::unit(strand_pos, "npc"),
                      r = grid::unit(abs(value), "npc"),
                      gp = grid::gpar(fill = col))
  }
  else if (type == "triangle") {
    sizefactor <- 0.05
    grid::grid.polygon(c(pos-(pos*sizefactor),pos,pos+(pos*sizefactor)),
                       c(strand_pos-(strand_pos*sizefactor),strand_pos+(strand_pos*sizefactor), strand_pos-(strand_pos*sizefactor)),
                       gp = grid::gpar(fill = col))
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
#' drawText(x_ = 1, y_ = 3, w_ = 1, h_ = 1, "Test", angle = 45)
drawText <- function(vp_name = NULL, x_, y_, w_, h_, label_txt = NULL, angle = 0, gp_ = NULL) {
    grid::pushViewport(grid::viewport(name = vp_name,
                                      x = grid::unit(x_, "npc"),
                                      y = grid::unit(y_, "npc"),
                                      width = w_,
                                      height = h_,
                                      just = c("bottom")))
    grid::grid.text(label_txt, gp = gp_, rot = angle)
    grid::popViewport(1)
}

#' prepare (order and subset) input data
#' @param dataset data.frame of input
#' @param chromosome target chromosome
#' @param begin start coordinate of target region
#' @param stop end coordinate of target region
#' @param strand strand direction
#' @examples
#' prepare(dataset, "Chr1A", 1000, 2000, "both")
prepare <- function(dataset,
                             chromosome,
                             begin,
                             stop,
                             strand = "both") {

    # check if input file is valid
    if (!is.null(dataset)) {
        needed_columns <- c("chr", "start", "end", "strand")
        if (!all(needed_columns %in% colnames(dataset))) {
            cat("Input has to contain at least ", needed_columns, "\n")
            cat("Found colnames: ", colnames(dataset), "\n")
            stop()
        }
    }

    # filter by start and end column and order dataset
    dataset <- dataset[which(dataset$chr == chromosome & dataset$start >= begin & dataset$end <= stop), ]
    dataset <- dataset[with(dataset, order(dataset$start, dataset$end)), ]

    # check strand information
    if (strand != "both") {
      dataset <- dataset[dataset$strand == strand, ]
    }
    # check if data frame not empty
    if (nrow(dataset) == 0) {
        cat("Input data frame is empty after filtering (", chromosome, ":", begin, "-", stop, ")\n", sep = "")
    }
    return (dataset)
}

#' get coordinates of viewports to draw on
#' @param x number of tracks to draw
#' @examples
#' getLayout(list(annoTrack, ...))
getLayout <- function(x) {
  coordinates <- list("height_of_vp" = 0.4,
                      "y_position_of_vp" = 0.3)

  if (x != 1) {
    coordinates$height_of_vp <- 0.95 / x
    coordinates$y_position_of_vp <- dropLast(seq(0.05, 1, coordinates$height_of_vp))
  }
  return (coordinates)
}

#' helper to retrieve last element of vector
#' @param x vector to get the last element from
#' @examples
#' last(c(1, 2, 3))
last <- function(x) return (tail(x, n = 1))

#' helper to retrieve first element of vector
#' @param x vector to get the first element from
#' @examples
#' first(c(1, 2, 3))
first <- function(x) return (x[1])

#' helper to drop last element of vector
#' @param x vector to drop last element
#' @examples
#' dropLast(c(1, 2, 3))
dropLast <- function(x) return (x[-length(x)])

#' calculate y scale labels of annotation tracks
#' @param x data frame to infer max and min values from
#' @param threshold maximal y value for reference
#' @examples
#' getAnnoYScale(dat, threshold = 10)
getAnnoYScale <- function(x, range_ = NULL) {
  if (identical(x, numeric(0))) return(c(0, 0.5, 1))
  # get min and max of y scale of annoTrack
  interval <- c(min(x), max(x))
  if (length(x) == 1) interval <- c(0, x)
  
  # test if user provided specific range
  if (!is.null(range_)) interval <- range_
  return (pretty(interval))
}

#' calculate y scale breaks of annotation tracks
#' @param x vector of y scale labels to infer breaks
#' @examples
#' getAnnoYBreaks(c(0:10))
getAnnoYBreaks <- function(x) {
  if (identical(x, numeric(0))) return(c(0, 0.5, 1))
  no_of_breaks <- ifelse(length(x) < 3, 3, length(x))
  return (seq(0, 1, 1 / (no_of_breaks - 1)))
}

#' calculate x scale label
#' @param start vector of x start values
#' @param end vector of x end values
#' @param upstream integer of bp upstream
#' @param downstream integer of bp downstream
#' @examples
#' getXLabel(10, 100, 20, 20)
getXLabel <- function(start, end, upstream, downstream) {
  min_ <- min(start) - upstream
  min_ <- ifelse(min_ < 0, 0, min_)
  max_ <- max(end) + downstream
  return (pretty(c(min_:max_)))
}

#' test if value is in specific range
#' @param x single value to test
#' @param min_ minimal allowed value
#' @param max_ maximal allowed value
#' @examples
#' inRange(10,5,15)
inRange <- function(x, min_, max_) { return(all(x >= min_ & x <= max_)) }

#' print values of tracks
#' @param object data.frame containing data
#' @examples
#' showValues(object)
showValues <- function(object) {
  ifelse(nrow(object) > 100,
         print(head(object)),
         print(object)
  )
}

#' check provided chromosomes
#' @param chromosomes vector of chromosomes
#' @examples
#' checkChromosomes(c("Chr1", "Chr1"))

checkChromosomes <- function(chromosomes) {
  found_chromosomes <- unique(chromosomes)
  if (length(found_chromosomes) > 1) {
    warning("Found more than one chromosome identifier, use first.")
    return(first(found_chromosomes))
  } else {
    return (found_chromosomes)
  }
}

#' calculate relative x position of feauture
#' @param x numeric absolute numeric value
#' @param xmin minimal value of scale
#' @param xmax maximal value of scale
#' @examples
#' relativePosition(100, 50, 150)
relativePosition <- function(x, xmin, xmax) {
  return ((x - xmin) / (xmax - xmin))
}

#' set track specific layout parameter
#' @param vp_y_position y position of viewport
#' @param vp_height height of viewport
#' @examples
#' track(...)
track = setClass("track",
                 slots = list(layout = "list",
                              xmin = "numeric",
                              xmax = "numeric",
                              upstream = "numeric",
                              downstream = "numeric"),
                 prototype = list(
                      layout = list(
                                  vp_y_position = 0.5,
                                  vp_height = 0.3,
                                  vp_width = 0.7
                      ),
                      xmin = 0,
                      xmax = 0,
                      upstream = 10,
                      downstream = 10
                 )
)

# create new environment to exchange variables
shared <- new.env()

#' draw all tracks
#' @param tracks list of tracks to plot
#' @examples
#' geneset(list(geneTrack(...), annoTrack(anno)))
geneset = setClass("geneset", slots = list(tracks = "list"))

# constructor method
geneset <- function(track) {
  .Object <- new("geneset")
  if(typeof(track) != "list") track <- list(track)
  .Object@tracks <- track
  return (.Object)
}

setMethod(f = "show",
          signature = "geneset",
          definition = function(object) {
            # get layout of plot
            layout <- getLayout(length(object@tracks))

            # create new device and newpage
            grid::grid.newpage()

            # outer viewport
            grid::pushViewport(grid::viewport(name = "outer",
                                              x = grid::unit(0.5, "npc"),
                                              y = grid::unit(0.5, "npc"),
                                              width = 1,
                                              height = 1))

            lapply(c(1:length(object@tracks)), function(x) {
                              object@tracks[[x]]@layout["vp_y_position"] <- layout$y_position_of_vp[[x]]
                              object@tracks[[x]]@layout["vp_height"] <- layout$height_of_vp
                              print(object@tracks[[x]])
            })
})

#' draw a data track based on data.frame
#' @param track_file data.frame with chr, start, stop and strand information
#' @param label name of the track
#' @param values column where the actual data is stored (default: value)
#' @param upstream distance in bp upstream
#' @param downstream distance in bp downstream
#' @param border draw border around annoTrack viewport
#' @param yrange range of y scale
#' @param label_gp gp object to edit label style
#' @param track_gp gp object to edit track appearence
#' @param label_orientation text orientation of track label
#' @param show_values logical if underlying values should be displayed
#' @examples
#' annoTrack(gff, "Coverage", "line", "firebrick")
annoTrack = setClass("annoTrack",
                     slots = list(
                       track_param = "list"
                     ),
                     contains = "track"
)

annoTrack <- function(track_file = NULL,
                      label = "Track",
                      values = "value",
                      upstream = 0,
                      downstream = 0,
                      border = TRUE,
                      yrange = NULL,
                      label_gp = grid::gpar(fontsize = 20, col = "black"),
                      track_gp = grid::gpar(col = "gray40", lwd = 1),
                      label_orientation = "h",
                      show_values = FALSE
                      ) {

    .Object <- new("annoTrack")

    # get environment variables
    xmin <- ifelse(exists("xmin", shared), get("xmin", shared), first(xLabel))
    xmax <- ifelse(exists("xmax", shared), get("xmax", shared), last(xLabel))
    chromosome <- ifelse(exists("chromosome", shared), get("chromosome", shared), checkChromosomes(track_file$chr))

    # prepare track file
    track_file <- prepare(track_file, chromosome, xmin, xmax)

    # get column with data
    if (values %in% names(track_file)) {
        names(track_file)[names(track_file) == values] <- "value"
    } else {
        cat("Column ", values, " not present in track file.\n")
        stop()
    }

    # get min and max of scale
    yscale_label <- getAnnoYScale(track_file$value, yrange)
    yscale_at <- getAnnoYBreaks(yscale_label)
    yinterval <- diff(range(yscale_label))
    # assign parameter to object
    .Object@track_param$track_file <- track_file
    if(show_values) showValues(.Object@track_param$track_file)
    .Object@xmin <- xmin
    .Object@xmax <- xmax
    .Object@upstream = upstream
    .Object@downstream = downstream
    .Object@track_param$yscale_label <- yscale_label
    .Object@track_param$yscale_at <- yscale_at
    .Object@track_param$yinterval <- ifelse(yinterval < 1, 1, yinterval)
    .Object@track_param$ymax <- max(yscale_label)
    .Object@track_param$ymin <- min(yscale_label)
    .Object@track_param$start_y <- ifelse(0 %in% yscale_label, yscale_at[which(yscale_label == 0)], 0)
    .Object@track_param$label <- label
    .Object@track_param$label_gp <- label_gp
    .Object@track_param$track_gp <- track_gp
    .Object@track_param$border <- border
    .Object@track_param$label_orientation <- ifelse(!(label_orientation %in% c("h", "v")),
                                                    "h",
                                                    label_orientation)
    return (.Object)
}

setMethod(f = "show",
          signature = "annoTrack",
          definition = function(object) {
            grid::pushViewport(grid::viewport(x = grid::unit(0.5, "npc"),
                                              y = grid::unit(object@layout["vp_y_position"], "npc"),
                                              width = grid::unit(object@layout["vp_width"], "npc"),
                                              height = as.numeric(object@layout["vp_height"]) - 0.025,
                                              just = c("bottom")))

            # draw border if requested (default)
            if (object@track_param$border) grid::grid.rect()

            # add track label
            grid::grid.text(object@track_param$label,
                            x = -0.075,
                            y = 0.5,
                            just = ifelse(object@track_param$label_orientation == "h", "right", "center"),
                            rot = ifelse(object@track_param$label_orientation == "h", 0, 90),
                            gp = object@track_param$label_gp)

            # add yaxis
            grid::grid.yaxis(label = object@track_param$yscale_label,
                             at = object@track_param$yscale_at,
                             gp = grid::gpar(fontsize = 16))

        if (nrow(object@track_param$track_file) != 0) {
            # add middle line
            grid::grid.segments(x0 = grid::unit(0, "npc"),
                                y0 = grid::unit(object@track_param$start_y, "npc"),
                                x1 = grid::unit(1, "npc"),
                                y1 = grid::unit(object@track_param$start_y,  "npc"),
                                gp = grid::gpar(col = "gray30"))

            left_border <- relativePosition(object@track_param$track_file$start, object@xmin, object@xmax)
            right_border <- relativePosition(object@track_param$track_file$end, object@xmin, object@xmax)
            width_ <- right_border - left_border
            height_ <- object@track_param$track_file$value / object@track_param$yinterval

            # which region should be displayed
            grid::grid.rect(x = grid::unit(left_border, "npc"),
                            y = grid::unit(object@track_param$start_y, "npc"),
                            width = grid::unit(width_, "npc"),
                            height = grid::unit(height_, "npc"),
                            just = c("left", "bottom"),
                            gp = object@track_param$track_gp)
        }
        grid::popViewport(1)
    })

#' draw a set of genes based on data.frame
#' @param track_file data.frame with chr, start, end and strand column
#' @param forward_color color of genes in forward direction (default = "darkslategray")
#' @param reverse_color color of genes in reverse direction (default = "darkslategray")
#' @param upstream distance in bp upstream
#' @param downstream distance in bp downstream
#' @param features data.frame of features to add to strand position
#' @param show_axis show axis or not (default = TRUE)
#' @param axis_label_text text of x axis label (default: "Region (bp)")
#' @param axis_label_offset offset of label (default: -0.5)
#' @param axis_label_gp gp object of x axis label
#' @param border boolean if border visible (default = FALSE)
#' @param show_values boolean if displayed range should also be printed
#' @examples
#' geneTrack(track_file)
geneTrack = setClass("geneTrack",
                   slots = list(
                        track_file = "ANY",
                        gene_param = "list",
                        plot_param = "list"
                   ),
                   prototype = list(
                    plot_param = list(
                      forward_strand_pos = 0.8,
                      single_strand_pos = 0.5,
                      reverse_strand_pos = 0.2
                    )
                   ),
                   contains = "track"
)
# constructor method
geneTrack <- function(track_file,
                      forward_color = "navajowhite3",
                      reverse_color = "darkslategray",
                      upstream = 0,
                      downstream = 0,
                      strand = "both",
                      features = NULL,
                      show_axis = TRUE,
                      show_direction_label = TRUE,
                      axis_label_text = NULL,
                      axis_label_offset = -0.8,
                      axis_label_gp = NULL,
                      border = FALSE,
                      show_values = FALSE
                      ) {

    .Object <- new("geneTrack")

    # assign data to object
    if (strand != "both") {
      track_file <- track_file[track_file$strand == strand, ]
    }
    .Object@track_file <- track_file
    if (show_values) showValues(.Object@track_file)

    # check for chromosome label
    .Object@gene_param$chromosome <- checkChromosomes(.Object@track_file$chr)

    # determine axis label
    .Object@upstream = upstream
    .Object@downstream = downstream
    axis_label <- getXLabel(.Object@track_file$start, .Object@track_file$end, .Object@upstream, .Object@downstream)
    .Object@xmin <- first(axis_label)
    .Object@xmax <- last(axis_label)

    assign("xmin", .Object@xmin, shared)
    assign("xmax", .Object@xmax, shared)
    assign("chromosome", .Object@gene_param$chromosome, shared)

    .Object@plot_param$axis_label_text <- ifelse(!is.null(axis_label_text),
                                                 axis_label_text,
                                                 paste(.Object@gene_param$chromosome, "(bp)"))
    if(!is.null(features)) {
      .Object@gene_param$features <- features
    }

    # set defaults
    .Object@gene_param$forward_color <- forward_color
    .Object@gene_param$reverse_color <- reverse_color
    .Object@plot_param$strand <- strand
    .Object@plot_param$show_axis <- show_axis
    .Object@plot_param$show_direction_label <- show_direction_label
    .Object@plot_param$axis_label <- axis_label
    .Object@plot_param$axis_label_offset <- axis_label_offset
    .Object@plot_param$axis_label_gp <- axis_label_gp
    .Object@plot_param$border <- border
    .Object@plot_param$show_values <- show_values

    return (.Object)
}

setMethod(f = "show",
          signature = "geneTrack",
          definition = function(object) {
            # main viewport
            grid::pushViewport(grid::viewport(name = "main_outer",
                                              x = grid::unit(0.5, "npc"),
                                              y = grid::unit(object@layout["vp_y_position"], "npc"),
                                              width = grid::unit(object@layout["vp_width"], "npc"),
                                              height = object@layout["vp_height"],
                                              just = c("bottom")))

            grid::pushViewport(grid::viewport(name = "main",
                                              x = grid::unit(0.5, "npc"),
                                              y = grid::unit(0.8, "npc"),
                                              width = 1,
                                              height = 0.5,
                                              just = c("top")))

            # draw border if requested
            if (object@plot_param$border) grid::grid.rect()

            # plot strand
            if (object@plot_param$strand == "both") {
              drawStrand(direction = "+",
                         y_ = object@plot_param$forward_strand_pos,
                         show_direction = object@plot_param$show_direction_label)
              drawStrand(direction = "-",
                         y_ = object@plot_param$reverse_strand_pos,
                         show_direction = object@plot_param$show_direction_label)

              strand_positions <- ifelse(object@track_file$strand == "+",
                                  object@plot_param$forward_strand_pos,
                                  object@plot_param$reverse_strand_pos)
            } else {
              drawStrand(direction = object@plot_param$strand,
                         y_ = object@plot_param$single_strand_pos,
                         show_direction = object@plot_param$show_direction_label)

              strand_positions <- object@plot_param$single_strand_pos
            }

            # add all genes/transcripts
            mapply(drawGene,
                   x1 = grid::unit(relativePosition(object@track_file$start, object@xmin, object@xmax), "npc"),
                   x2 = grid::unit(relativePosition(object@track_file$end, object@xmin, object@xmax), "npc"),
                   pos =  strand_positions,
                   direction = object@track_file$strand,
                   forward_color = object@gene_param$forward_color,
                   reverse_color = object@gene_param$reverse_color)

            # add strand specific features
            if (!is.null(object@gene_param$features)) {
              drawFeature(relativePosition(object@gene_param$features$pos, object@xmin, object@xmax),
                          ifelse(object@gene_param$features$strand == "+", object@plot_param$forward_strand_pos, object@plot_param$reverse_strand_pos),
                          object@gene_param$features$value,
                          object@gene_param$features$col,
                          object@gene_param$features$type)
            }

            # add axis label
            if (object@plot_param$show_axis) {
                grid::grid.xaxis(label = object@plot_param$axis_label,
                                 at = seq(0, 1, 1 / (length(object@plot_param$axis_label) - 1)))

                # add axis label text
                drawText(x_ = 0.5,
                         y_ = object@plot_param$axis_label_offset,
                         w_ = 0.1,
                         h_ = 0.2,
                         label_txt = object@plot_param$axis_label_text,
                         gp_ = object@plot_param$axis_label_gp)
            }
            grid::popViewport(2)
})
