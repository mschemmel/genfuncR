#' genearrow plots an arrow in an defined direction representing a single gene
#' @param x1 x1 coordinate
#' @param x2 x2 coordinate
#' @param pos positional coordinates of forward and reverse strand
#' @param direction direction of the gene (upstream | downstream)
#' @param gd_ grid gpar object to edit the arrow appearance
#' @examples
#' genearrow(x1 = 1, y1 = 5, x2 = 4, y2 = 6, arrow_type = "arrow", direction = "downstream")

genearrow <- function(x1, x2, pos, direction, forward_color = "darkslategray", reverse_color = "navajowhite3", gene_height = 0.15) {
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
#' textLabel(x_ = 1, y_ = 3, w_ = 1, h_ = 1, "Test", angle = 45)

textLabel <- function(vp_name = NULL, x_, y_, w_, h_, label_txt = NULL, angle = 0, gp_ = NULL) {
    grid::pushViewport(grid::viewport(name = vp_name,
                                      x = grid::unit(x_, "npc"),
                                      y = grid::unit(y_, "npc"),
                                      width = w_,
                                      height = h_,
                                      just = c("bottom")))
    grid::grid.text(label_txt, gp = gp_, rot = angle)
    grid::popViewport(1)
}

#' prepareAndFilter (order and subset) input data
#' @param dataset data.frame of input
#' @param chromosome target chromosome
#' @param st start coordinate of target region
#' @param en end coordinate of target region
#' @examples
#' prepareAndFilter(gff, "Chr1A", 1000, 2000)

prepareAndFilter <- function(dataset,
                             chromosome,
                             st,
                             en) {

    # check if input file is valid
    if (!is.null(dataset)) {
        needed_columns <- c("chr", "start", "end", "strand")
        if (!all(needed_columns %in% colnames(dataset))) {
            cat("Input has to contain at least ", needed_columns, "\n")
            cat("Found colnames: ", colnames(dataset), "\n")
            stop()
        }
    }

    # order input by start and end column and filter
    dataset <- dataset[with(dataset, order(dataset$start, dataset$end)), ]
    dataset <- dataset[dataset$chr == chromosome & dataset$start >= st & dataset$end <= en, ]

    # check if data frame not empty
    if (nrow(dataset) == 0) {
        cat("Input data frame is empty after filtering (", chromosome, ":", st, "-", en, ")\n", sep = "")
    }
    return (dataset)
}

#' get coordinates of viewports to draw on
#' @param x number of tracks to draw
#' @examples
#' getLayout(list(annoTrack, ...))
getLayout <- function(x) {
  coordinates <- list("size_per_vp" = 0.3,
                      "places_of_vp" = 0.5)

  if (x != 1) {
    coordinates$size_per_vp <- 0.85 / x
    coordinates$places_of_vp <- head(seq(0.15, 1, coordinates$size_per_vp), n = -1)
  }
  return(coordinates)
}

#' helper to retrieve last element of vector
#' @param x vector to get the last element from
#' @examples
#' last(c(1,2,3))
last <- function(x) return (tail(x, n = 1))

#' helper to retrieve first element of vector
#' @param x vector to get the first element from
#' @examples
#' first(c(1,2,3))
first <- function(x) return (x[1])

#' calculate y scale labels of annotation tracks
#' @param x data frame to infer max and min values from
#' @param threshold maximal y value for reference
#' @examples
#' getAnnoYScale(dat, threshold = 10)

getAnnoYScale <- function(x, threshold = NULL) {
  # get min and max of y scale of annoTrack
  yMinScale <- ifelse(min(x$value) >= 0, 0, min(x$value))
  yMaxScale <- ifelse(max(x$value) >= 0, max(x$value), 0)
  if (!is.null(threshold)) {
    yMaxScale <- ifelse(threshold < yMaxScale, yMaxScale,threshold)
  }
  yScaleLabel <- pretty(c(yMinScale:yMaxScale))
  return(yScaleLabel)
}

#' calculate y scale breaks of annotation tracks
#' @param x vector of y scale labels to infer breaks
#' @examples
#' getAnnoYBreaks(c(0:10))
getAnnoYBreaks <- function(x) {
  return(seq(0, 1, 1 / (length(x) - 1)))
}

#' set track specific layout parameter
#' @param vp_y_position y position of viewport
#' @param vp_height height of viewport
#' @examples
#' track(...)

track = setClass("track",
                 slots = list(vp_y_position = "numeric",
                              vp_height = "numeric"),
                 prototype = list(vp_y_position = 0.5,
                                  vp_height = 0.3)
)

# create new environment to exchange variables
shared <- new.env()

#' draw all tracks
#' @param locus geneTrack object to draw
#' @param annotation single or list of annoTracks to draw
#' @examples
#' geneset(locus, annotation)
geneset = setClass("geneset", slots = list(tracks = "ANY"))

# constructor method
geneset <- function(locus, annotation = NULL) {
  .Object <- new("geneset")
  .Object@tracks$locus <- locus
  .Object@tracks$annotation <- annotation
  .Object@tracks$no_of_tracks <- sum(length(locus), length(annotation))

  return(.Object)
}

setMethod(f = "show",
          signature = "geneset",
          definition = function(object) {
            layout <- getLayout(object@tracks$no_of_tracks)

            # draw geneTrack
            object@tracks$locus@vp_y_position <- first(layout$places_of_vp)
            object@tracks$locus@vp_height <- layout$size_per_vp
            print(object@tracks$locus)

            # draw annoTrack(s)
            if (!is.null(object@tracks$annotation)) {
              if (typeof(object@tracks$annotation) != "list") {
                object@tracks$annotation@vp_y_position <- last(layout$places_of_vp)
                object@tracks$annotation@vp_height <- layout$size_per_vp
                print(object@tracks$annotation)
              }
              else {
                lapply(seq_along(object@tracks$annotation), function(x) {
                  object@tracks$annotation[[x]]@vp_y_position <- layout$places_of_vp[-1][[x]]
                  object@tracks$annotation[[x]]@vp_height <- layout$size_per_vp
                  print(object@tracks$annotation[[x]])
                })
              }
            }
})

#' draw a data track based on data.frame
#' @param track_file data.frame with chr, start, stop and strand information
#' @param label name of the track
#' @param values column where the actual data is stored (default: value)
#' @param border draw border around annoTrack viewport
#' @param ymax max data value for y axis label
#' @param label_gp gp object to edit label style
#' @param track_gp gp object to edit track appearence
#' @param label_orientation text orientation of track label
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
                      border = TRUE,
                      ymax = NULL,
                      label_gp = grid::gpar(fontsize = 12, col = "black"),
                      track_gp = grid::gpar(col = "gray40", lwd = 1),
                      label_orientation = "horizontal"
                      ) {

    .Object <- new("annoTrack")

    # get column with data
    if (values %in% names(track_file)) {
        names(track_file)[names(track_file) == values] <- "value"
    } else {
        cat("Column ", values, " not present in track file.\n")
        stop()
    }

    # get min and max of scale
    yscale_label <- getAnnoYScale(track_file, ymax)
    yscale_at <- getAnnoYBreaks(yscale_label)

    # set environment variables
    xmin <- shared$min_value
    xmax <- shared$max_value
    chromosome <- shared$chromosome

    # assign parameter to object
    .Object@track_param$track_file <- prepareAndFilter(track_file, chromosome, xmin, xmax)
    .Object@track_param$xmin <- xmin
    .Object@track_param$xmax <- xmax
    .Object@track_param$yscale_label <- yscale_label
    .Object@track_param$yscale_at <- yscale_at
    .Object@track_param$ymax <- last(yscale_label)
    .Object@track_param$start_y <- 0
    .Object@track_param$label <- label
    .Object@track_param$label_gp <- label_gp
    .Object@track_param$track_gp <- track_gp
    .Object@track_param$border <- border

    if (!(label_orientation %in% c("vertical", "horizontal"))) {
        .Object@track_param$label_orientation <- "horizontal"
        cat("DataTrack: ", label, " Unknown label_orientation value. Set to 'horizontal'")
    } else {
        .Object@track_param$label_orientation <- label_orientation
    }
    return(.Object)
}

setMethod(f = "show",
          signature = "annoTrack",
          definition = function(object) {
            relative <- function(x) (x - object@track_param$xmin) / (object@track_param$xmax - object@track_param$xmin)
            grid::pushViewport(grid::viewport(x = grid::unit(0.5, "npc"),
                                              y = grid::unit(object@vp_y_position, "npc"),
                                              width = 0.7,
                                              height = object@vp_height - 0.02,
                                              just = c("bottom")))

            # draw border if requested (default)
            if (object@track_param$border) grid::grid.rect()

            # add track label
            grid::grid.text(object@track_param$label,
                            x = -0.1,
                            y = 0.5,
                            just = ifelse(object@track_param$label_orientation == "horizontal", "right", "center"),
                            rot = ifelse(object@track_param$label_orientation == "horizontal", 0, 90),
                            gp = object@track_param$label_gp)

            # add yaxis
            grid::grid.yaxis(label = object@track_param$yscale_label,
                             at = object@track_param$yscale_at,
                             gp = grid::gpar(fontsize = 8))

        if (nrow(object@track_param$track_file) != 0) {
            if (!all(object@track_param$track_file$value >= 0)) {
                object@track_param$start_y <- 0.5
                # add middle line
                grid::grid.segments(x0 = grid::unit(0, "npc"),
                                    y0 = grid::unit(object@track_param$start_y, "npc"),
                                    x1 = grid::unit(1, "npc"),
                                    y1 = grid::unit(object@track_param$start_y,  "npc"),
                                    gp = grid::gpar(col = "gray30"))
            }

            # which region should be displayed
            grid::grid.segments(x0 = grid::unit(relative(object@track_param$track_file$start), "npc"),
                                y0 = grid::unit(object@track_param$start_y, "npc"),
                                x1 = grid::unit(relative(object@track_param$track_file$end), "npc"),
                                y1 = grid::unit(object@track_param$start_y + (object@track_param$track_file$value / object@track_param$ymax), "npc"),
                                gp = object@track_param$track_gp)
        }
        grid::popViewport(1)
    })

#' draw a set of genes based on data.frame
#' @param track_file data.frame with chr, start, end and strand column
#' @param forward_color color of genes in forward direction (default = "darkslategray")
#' @param reverse_color color of genes in reverse direction (default = "darkslategray")
#' @param transparency alpha value of gene annotation box
#' @param gene_height height in percent (0-1) of gene box (default = 1)
#' @param distance distance between forward and reverse strand in percent (0-1) (default = 1)
#' @param show_axis show axis or not (default = TRUE)
#' @param axis_label_text text of x axis label (default: "Region (bp)")
#' @param axis_label_offset offset of label (default: -0.5)
#' @param axis_label_gp gp object of x axis label
#' @param border boolean if border visible (default = FALSE)
#' @param show_values boolean if displayed range should also be printed
#' @examples
#' geneTrack(gff)

geneTrack = setClass("geneTrack",
                   slots = list(
                        track_file = "ANY",
                        gene_param = "list",
                        plot_param = "list"
                   ),
                   contains = "track"
)
# constructor method
geneTrack <- function(track_file,
                      forward_color = "darkslategray",
                      reverse_color = "navajowhite3",
                      upstream = 10,
                      downstream = 10,
                      transparency = 1,
                      gene_height = 1,
                      distance = 1,
                      show_axis = TRUE,
                      axis_label_text = NULL,
                      axis_label_offset = -0.5,
                      axis_label_gp = NULL,
                      border = FALSE,
                      show_values = FALSE
                      ) {


    .Object <- new("geneTrack")
    # assign data to object
    .Object@track_file <- track_file
    if (show_values) {
        if (nrow(.Object@track_file) > 100) {
            print(head(.Object@track_file))
        } else {
            print(.Object@track_file)
        }
    }

    # check for unique chromosome label
    given_chromosome <- unique(.Object@track_file$chr)
    if (!is.null(given_chromosome)) {
      if (length(given_chromosome) > 1) {
        warning("Found more than one chromosome identifier.")
        .Object@gene_param$chromosome <- first(given_chromosome)
      } else {
        .Object@gene_param$chromosome <- given_chromosome
      }
    }

    shared$chromosome <- given_chromosome

    # determine axis label
    axis_label <- pretty(c(min(.Object@track_file$start - upstream):max(.Object@track_file$end + downstream)))
    min_value <- first(axis_label)
    max_value <- last(axis_label)
    shared$min_value <- min_value
    shared$max_value <- max_value

    # overall size of genebox
    gene_box_height <- ifelse(gene_height > 0 & gene_height <= 1,
                              0.2 * gene_height,
                              stop("Height of gene box (gene_height) have to be between 0 and 1."))
    .Object@gene_param$gene_box_height <- gene_box_height

    # position of forward and reverse strand
    forward_strand_pos <- 0.2
    strand_gap <- ifelse(distance >= 0 & distance <= 1,
                         gene_box_height + ((forward_strand_pos + gene_box_height) * distance),
                         stop("Distance between forward and reverse strand have to be between 0 and 1."))
    reverse_strand_pos <- forward_strand_pos + strand_gap

    .Object@plot_param$axis_label_text <- ifelse(!is.null(axis_label_text),
                                                 axis_label_text,
                                                 paste(.Object@gene_param$chromosome, "(bp)"))

    # default params
    .Object@gene_param$forward_strand_pos <- forward_strand_pos
    .Object@gene_param$reverse_strand_pos <- reverse_strand_pos
    .Object@gene_param$distance <- distance
    .Object@gene_param$forward_color <- forward_color
    .Object@gene_param$reverse_color <- reverse_color
    .Object@gene_param$transparency <- transparency
    .Object@gene_param$gene_height <- gene_height
    .Object@plot_param$min_value <- min_value
    .Object@plot_param$max_value <- max_value
    .Object@plot_param$show_axis <- show_axis
    .Object@plot_param$axis_label <- axis_label
    .Object@plot_param$axis_label_offset <- axis_label_offset
    .Object@plot_param$axis_label_gp <- axis_label_gp
    .Object@plot_param$border <- border
    .Object@plot_param$show_values <- show_values

    return(.Object)
}

setMethod(f = "show",
          signature = "geneTrack",
          definition = function(object) {
            # helper function
            relative <- function(x) (x - object@plot_param$min_value) / (object@plot_param$max_value - object@plot_param$min_value)

            # create new device and newpage
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
                                              y = grid::unit(object@vp_y_position, "npc"),
                                              width = 0.7,
                                              height = object@vp_height,
                                              just = c("bottom")))

            # draw border if requested
            if (object@plot_param$border) grid::grid.rect()

            # forward direction
            grid::grid.segments(x0 = grid::unit(0, "npc"),
                                y0 = grid::unit(object@gene_param$forward_strand_pos, "npc"),
                                x1 = grid::unit(1, "npc"),
                                y1 = grid::unit(object@gene_param$forward_strand_pos, "npc"))

            # reverse direction
            grid::grid.segments(x0 = grid::unit(0, "npc"),
                                y0 = grid::unit(object@gene_param$reverse_strand_pos, "npc"),
                                x1 = grid::unit(1, "npc"),
                                y1 = grid::unit(object@gene_param$reverse_strand_pos, "npc"))

            # add all genes/transcripts
            mapply(genearrow,
                   x1 = grid::unit(relative(object@track_file$start), "npc"),
                   x2 = grid::unit(relative(object@track_file$end), "npc"),
                   pos = ifelse(object@track_file$strand == "+",
                                object@gene_param$reverse_strand_pos,
                                object@gene_param$forward_strand_pos),
                   direction = object@track_file$strand,
                   forward_color = object@gene_param$forward_color,
                   reverse_color = object@gene_param$reverse_color)

            # add axis label
            if (object@plot_param$show_axis) {
                # add axis label
                grid::grid.xaxis(label = object@plot_param$axis_label,
                                 at = seq(0, 1, 1 / (length(object@plot_param$axis_label) - 1)))

                # add axis label text
                textLabel(x_ = 0.5,
                          y_ = object@plot_param$axis_label_offset,
                          w_ = 0.1,
                          h_ = 0.2,
                          label_txt = object@plot_param$axis_label_text,
                          gp_ = object@plot_param$axis_label_gp)
            }
            grid::popViewport(1)
})