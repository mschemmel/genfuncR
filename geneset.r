#' plots an arrow in an defined direction representing a single gene
#' @param x1 x1 coordinate
#' @param y1 y1 coordinate
#' @param x2 x2 coordinate
#' @param y2 y2 coordinate
#' @param direction direction of the gene (upstream | downstream)
#' @param arr_type appearance of arrow head (arrow | barrow | box)
#' @examples
#' genearrow(x1 = 1, y1 = 5, x2 = 4, y2 = 6, arrow_type = "arrow", direction = "downstream")

genearrow <- function(x1, x2, pos, direction, arr_type = "arrow", gp_) {
    arrow_width <- x2 - x1
    arrow_head_width <- arrow_width * 0.2 # TODO: make size dynamic
    barrow_head_size <- 0.05
    
    y1 <- pos - (pos * 0.1) 
    y2 <- pos + (pos * 0.1) 

    if (direction == "+") {
        if (arr_type == "arrow") {
            grid::grid.polygon(c(x1, x1, x2 - arrow_head_width, x2, x2 - arrow_head_width),
                               c(y1, y2, y2, mean(c(y1, y2)), y1),
                               gp = gp_)
        } else if (arr_type == "barrow") {
            grid::grid.polygon(c(x1, x1, x2 - arrow_head_width, x2 - arrow_head_width, x2, x2 - arrow_head_width, x2 - arrow_head_width),
                               c(y1, y2, y2, y2 + barrow_head_size, mean(c(y1, y2)), y1 - barrow_head_size, y1),
                               gp = gp_)
        } else if (arr_type == "box") {
            grid::grid.polygon(c(x1, x1, x2, x2),
                               c(y1, y2, y2, y1),
                               gp = gp_)
        } else {
           stop("Invalid 'arr-' parameter. Please choose 'arrow', 'barrow' or 'box'.")
        }
    } else if (direction == "-") {
        if (arr_type == "arrow") {
            grid::grid.polygon(c(x1 + arrow_head_width, x1, x1 + arrow_head_width, x2, x2),
                               c(y1, mean(c(y2, y1)), y2, y2, y1),
                               gp = gp_)
        } else if (arr_type == "barrow") {
            grid::grid.polygon(c(x1 + arrow_head_width, x1 + arrow_head_width, x1, x1 + arrow_head_width, x1 + arrow_head_width, x2, x2),
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
                                      just = c("bottom")))
    grid::grid.text(label_txt, gp = gp_, rot = angle)
    grid::popViewport(1)
}

#' prepare (order and subset) input data
#' @param dataset data.frame of input
#' @param chromosome target chromosome
#' @param st start coordinate of target region
#' @param en end coordinate of target region
#' @param transparency alpha value of gene an
#' @examples
#' prepare(gff, "Chr1A", 1000, 2000)

prepareAndFilter <- function(dataset,
                             chromosome,
                             st,
                             en) {

    # check if input file is valid
    if (!is.null(dataset)) {
        if (!all(c("chr", "start", "end", "strand") %in% colnames(dataset))) {
            cat("Input has to contain at least 'chr', 'start', 'end' and 'strand' column.\n")
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
    return(dataset)
}

#' draw a data track based on data.frame or gff file
#' @param track_file data.frame of gff file
#' @param label name of the track
#' @param type type of the track (default: line)
#' @param values column where the actual data is stored (default: value)
#' @param ymax max data value for y axis label
#' @param label_gp gp object to edit label style
#' @param track_gp gp object to edit track appearence
#' @param label_orientation text orientation of track label
#' @examples
#' annoTrack(gff, "Coverage", "line", "firebrick")

annoTrack = setClass("annoTrack",
                     slots = list(
                       track_param = "list"
                    ))

annoTrack <- function(track_file = NULL,
                      label = "Track",
                      type = "s",
                      values = "value",
                      label_gp = grid::gpar(fontsize = 10, col = "black"),
                      track_gp = grid::gpar(col = "gray40", lwd = 1),
                      label_orientation = "horizontal"
                      ) {

    .Object = new("annoTrack")

    if (values %in% names(track_file)) {
        names(track_file)[names(track_file) == values] <- "value"
    }
    else {
        cat("Column ", values, " not present in track file.\n")
        stop()
    }

    # get min and max of scale
    minscale <- ifelse(min(track_file$value) >= 0, 0, min(track_file$value))
    maxscale <- ifelse(max(track_file$value) >= 0, max(track_file$value), 0)
    yscale_label <- pretty(c(minscale:maxscale))
    scale_interval <- diff(range(yscale_label))
    yscale_at <- seq(0, 1, 1/(length(yscale_label)-1))

    # assign parameter to object
    .Object@track_param$scale_label = yscale_label
    .Object@track_param$scale_at = yscale_at
    .Object@track_param$ymax = tail(yscale_label, n = 1)
    .Object@track_param$track_file = track_file
    .Object@track_param$label = label
    .Object@track_param$type = type
    .Object@track_param$label_gp = label_gp
    .Object@track_param$track_gp = track_gp

    if (!(label_orientation %in% c("vertical", "horizontal"))) {
        .Object@track_param$label_orientation = "horizontal"
        cat("DataTrack: ", label, " Unknown label_orientation value. Set to 'horizontal'")
    }
    else {
        .Object@track_param$label_orientation = label_orientation
    }
    return(.Object)
}

#' draw a set of genes based on data.frame or gff file
#' @param gff_file data.frame of gff file
#' @param forward_color color of genes in forward direction (default = "darkslategray")
#' @param reverse_color color of genes in reverse direction (default = "darkslategray")
#' @param upstream distance (bp) upstream
#' @param downstream distance (bp) downstream
#' @param transparency alpha value of gene annotation box
#' @param arrow_type shape of arrow head (default = "arrow")
#' @param gene_height height in percent (0-1) of gene box (default = 1)
#' @param distance distance between forward and reverse strand in percent (0-1) (default = 1)
#' @param show_axis show axis or not (default = TRUE)
#' @param axis_label_text text of x axis label (default: "Region (bp)")
#' @param axis_label_offset offset of label (default: -0.5)
#' @param axis_label_gp gp object of x axis label
#' @param border boolean if border visible (default = FALSE)
#' @param show_values boolean if displayed range should also be printed
#' @param tracks named list of data tracks (annotations)
#' @param marker position of highlight annotation
#' @examples
#' geneset(gff)

geneset = setClass("geneset",
                   slots = list(
                        gff_file = "ANY",
                        gene_param = "list",
                        plot_param = "list"
                   )
)
# constructor method
geneset <- function(gff_file,
                    forward_color = "darkslategray",
                    reverse_color = "darkslategray",
                    upstream = 0,
                    downstream = 0,
                    transparency = 1,
                    arrow_type = "arrow",
                    gene_height = 1,
                    distance = 1,
                    show_axis = TRUE,
                    axis_label_text = NULL,
                    axis_label_offset = -0.5,
                    axis_label_gp = NULL,
                    border = FALSE,
                    show_values = FALSE,
                    tracks = NULL,
                    marker = NULL) {

    .Object = new("geneset")

    # assign data to object
    .Object@gff_file = gff_file
    if (show_values) {
        if (nrow(.Object@gff_file) > 100) {
            print(head(.Object@gff_file))
        }
        else {
            print(.Object@gff_file)
        }
    }
    
    # check for unique chromosome label
    given_chromosome <- unique(.Object@gff_file$chr)
    if (!is.null(given_chromosome)) {
      if (length(given_chromosome) > 1) {
        warning("Found more than one chromosome identifier, use first provided.")
        .Object@gene_param$chromosome = given_chromosome[1] 
      }
      else {
        .Object@gene_param$chromosome = given_chromosome
      }
    }
    else {
      warning("No chromosome specified in input data.")
    }

    # determine axis label
    axis_label <- pretty(c(min(.Object@gff_file$start-upstream):max(.Object@gff_file$end+downstream)))
    min_value <- axis_label[1]
    max_value <- tail(axis_label, n = 1)
   
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

    .Object@plot_param$axis_label_text = ifelse(!is.null(axis_label_text),
                                                axis_label_text,
                                                paste(.Object@gene_param$chromosome, "(bp)"))

    # default params
    .Object@gene_param$forward_strand_pos = forward_strand_pos
    .Object@gene_param$reverse_strand_pos = reverse_strand_pos
    .Object@gene_param$distance = distance
    .Object@gene_param$forward_color = forward_color
    .Object@gene_param$reverse_color = reverse_color
    .Object@gene_param$transparency = transparency
    .Object@gene_param$arrow_type = arrow_type
    .Object@gene_param$gene_height = gene_height
    .Object@plot_param$min_value = min_value
    .Object@plot_param$max_value = max_value
    .Object@plot_param$show_axis = show_axis
    .Object@plot_param$axis_label = axis_label
    .Object@plot_param$axis_label_offset = axis_label_offset
    .Object@plot_param$axis_label_gp = axis_label_gp
    .Object@plot_param$border = border
    .Object@plot_param$show_values = show_values
    .Object@plot_param$tracks = tracks
    .Object@plot_param$marker = marker

    return(.Object)
}

setMethod(f = "show",
          signature = "geneset",
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
            
            # placement of track if requested
            if (length(object@plot_param$tracks) != 0) {
                object@plot_param$show_tracks = TRUE
                size_per_vp <- 0.85 / (length(object@plot_param$tracks) + 1) # always show the genebox (+1)
                places_of_vp <- head(seq(0.15, 1, size_per_vp), -1)
            }
            else {
                object@plot_param$show_tracks = FALSE
                size_per_vp <- 0.3
                places_of_vp <- 0.5
            }

            # main viewport
            grid::pushViewport(grid::viewport(name = "main",
                                              x = grid::unit(0.5, "npc"),
                                              y = grid::unit(places_of_vp[1], "npc"),
                                              width = 0.7,
                                              height = size_per_vp,
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
            genearrow(x1 = grid::unit(relative(object@gff_file$start), "npc"),
                      x2 = grid::unit(relative(object@gff_file$end), "npc"),
                      pos = ifelse(object@gff_file$strand == "+", object@gene_param$forward_strand_pos, object@gene_param$reverse_strand_pos),
                      direction = object@gff_file$strand, 
                      gp_ = grid::gpar(fill = "red",
                                       alpha = object@gene_param$transparency),
                      arr_type = object@gene_param$arrow_type)

            # add axis label
            if (object@plot_param$show_axis) {
                # set axis label
                axis_label <- object@plot_param$axis_label
                grid::grid.xaxis(label = axis_label,
                                 at = seq(0, 1, 1 / (length(axis_label) - 1)))

                # add axis label text
                text_label(x_ = 0.5,
                           y_ = object@plot_param$axis_label_offset,
                           w_ = 0.1,
                           h_ = 0.2,
                           label_txt = object@plot_param$axis_label_text,
                           gp_ = object@plot_param$axis_label_gp)
            }
                    
            grid::popViewport(1)

            # plot annotation tracks if requested
            # TODO: refactor!
            if (object@plot_param$show_tracks) {
                # how many tracks should be drawn
                for (x in seq_len(length(object@plot_param$tracks))) {
                    grid::pushViewport(grid::viewport(x = grid::unit(0.5, "npc"),
                                                      y = grid::unit(places_of_vp[-1][x], "npc"),
                                                      width = 0.7,
                                                      height = size_per_vp - 0.02,
                                                      just = c("bottom")))

                    dframe <- prepareAndFilter(object@plot_param$tracks[[x]]@track_param$track_file,
                                               object@gene_param$chromosome,
                                               object@plot_param$min_value,
                                               object@plot_param$max_value)
                    
                    if(object@plot_param$show_values) print(dframe)

                    grid::grid.rect()
                    
                    # add track label
                    grid::grid.text(object@plot_param$tracks[[x]]@track_param$label,
                                    x = ifelse(object@plot_param$tracks[[x]]@track_param$label_orientation == "horizontal", -0.05, -0.1),
                                    y = 0.5,
                                    just = ifelse(object@plot_param$tracks[[x]]@track_param$label_orientation == "horizontal", "right", "center"),
                                    rot = ifelse(object@plot_param$tracks[[x]]@track_param$label_orientation == "horizontal", 0, 90),
                                    gp = object@plot_param$tracks[[x]]@track_param$label_gp)

                    # add yaxis
                    grid::grid.yaxis(label = object@plot_param$tracks[[x]]@track_param$scale_label,
                                     at = object@plot_param$tracks[[x]]@track_param$scale_at,
                                     gp = grid::gpar(fontsize = 8))

        if (!is.null(object@plot_param$marker)) {
          mark = object@plot_param$marker
                            grid::grid.rect(x = grid::unit(relative(mark[1]), "npc"),
                                            y = grid::unit(0.5, "npc"),
                                            width = grid::unit(relative(mark[2]) - relative(mark[1]), "npc"),
                                            height = grid::unit(1,  "npc"),
                                            gp = grid::gpar(col = "red", lwd = 1),
                                            just = "left")
      
              }
                    if (nrow(dframe) != 0) {
                        # get maximum value of data as reference
                        start_y <- 0
                        max_in_range <- max(dframe$value)
                        
                        if (!all(dframe$value >= 0)) {
                            max_in_range <- 2
                            start_y <-  0.5
                            # add middle line
                            grid::grid.segments(x0 = grid::unit(0, "npc"),
                                                y0 = grid::unit(start_y, "npc"),
                                                x1 = grid::unit(1, "npc"),
                                                y1 = grid::unit(start_y,  "npc"),
                                                gp = grid::gpar(col = "gray30"))
                        }

                        # which region should be displayed
                        grid::grid.segments(x0 = grid::unit(relative(dframe$start), "npc"),
                                            y0 = grid::unit(start_y, "npc"),
                                            x1 = grid::unit(relative(dframe$end), "npc"),
                                            y1 = grid::unit(start_y + (dframe$value/max_in_range), "npc"),
                                            gp = object@plot_param$tracks[[x]]@track_param$track_gp)
                    }
                    grid::popViewport(1)
                }
            }
})
