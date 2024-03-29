#' puts a text label on a specific position
#' @param vp_name name of the viewport
#' @param x_ x coordinate
#' @param y_ y coordinate
#' @param w_ w coordinate
#' @param h_ h coordinate
#' @param label_txt text to be printed on that position
#' @param angle angle of text label
#' @param gp_ grid parameter like font, style, color, ...
#' @param cex fontsize of nucleotide
#' @param bg background color of nucleotide
#' @param borderCol color of border of nucleotide
#' @examples
#' textLabel(x_ = 1, y_ = 3, w_ = 1, h_ = 1, label_txt = "Test", angle = 45)

textLabel <- function(vp_name = NULL, x_, y_, w_, h_, label_txt = NULL, angle = 0, gp_ = NULL, cex = 0.4, bg, borderCol) {
    grid::pushViewport(grid::viewport(name = vp_name,
                                      x = grid::unit(x_, "npc"),
                                      y = grid::unit(y_, "npc"),
                                      width = w_,
                                      height = h_,
                                      just = c("bottom"),
                                      gp = grid::gpar(fontsize = cex)))
    grid::grid.rect(gp = grid::gpar(fill = bg, col = borderCol))
    grid::grid.text(label_txt, gp = gp_, rot = angle)
    grid::popViewport(1)
}

#' get width and heigth of text label of given fontsize
#' @param fontsize fontsize of nucleotide
#' @examples
#' nuclDims("A")
nuclDims <- function(fontsize) {
    nucleotides <- c("A", "T", "C", "G", "N")
    nucl_width <- max(unlist(lapply(nucleotides, function(x) strwidth(x, cex = fontsize, units = "inches"))))
    nucl_height <- max(unlist(lapply(nucleotides, function(x) strheight(x, cex = fontsize, units = "inches"))))
    return(c(nucl_width, nucl_height))
}

#' calculate x coordinates of nucleotides
#'
#' @param seq_ character string of nucleotides
#' @param letter number of letters per line
#' @examples
#' xCoords("ATG", 10)
xCoords <- function(seq_, letter) {
    return(rep(0.5:letter, ceiling(length(seq_) / letter))[1:length(seq_)] / letter)
}

#' calculate y coordinates of nucleotides
#'
#' @param seq_ character string of nucleotides
#' @param letter number of letters per line
#' @examples
#' yCoords("ATG", 10)
yCoords <- function(seq_, letter) {
    ypositions <- seq(0.9, 0.1, -0.05)
    return(rep(ypositions[1:ceiling(length(seq_) / letter)], each = letter)[1:length(seq_)])
}

#' plots nucleotide sequence as image
#'
#' @param sequence A nucleotide sequence string
#' @param base_colors A list of colors to be used for every nucleotide
#' @param cex size of font of single nucleotide
#' @param letter number of nucleotides per line
#' @param background background color of nucleotide
#'
#' @examples
#' dna_to_img("TATCGATCGATC", list(A = "#9EE362", T = "#00C0D0", G = "#FFD403", C = "#FF9356", U = "#d83131"), 50, "white")

dna2img <- setClass("dna2img",
                    slots = list(sequence = "character",
                                 base_colors = "list",
                                 cex = "numeric",
                                 letter = "numeric",
                                 background = "character",
                                 border = "character"),
                    prototype = list(base_colors = list("A" = "#9EE362",
                                                        "T" = "#00C0D0",
                                                        "G" = "#FFD403",
                                                        "C" = "#FF9356",
                                                        "U" = "#d83131",
                                                        "N" = "gray60"),
                                     cex = 0.4,
                                     letter = 20,
                                     background = "gray80",
                                     border = "white"))

dna2img <- function(sequence,
                    base_colors,
                    cex = 0.4,
                    letter = 20,
                    background = "gray80",
                    border = "white") {
    .Object <- new("dna2img")
    .Object@sequence <- unlist(strsplit(sequence, split = ""))
    .Object@letter <- letter
    .Object@cex <- cex
    .Object@background <- background
    .Object@border <- border
    return(.Object)
}

setMethod(f = "show",
          signature = "dna2img",
          definition = function(object) {

    grid::grid.newpage()

    # get coordinates of every nucleotide
    coordx <- xCoords(object@sequence, object@letter)
    coordy <- yCoords(object@sequence, object@letter)

    # draw all nucleotides
    grid::pushViewport(grid::viewport(x = 0.5, y = 0.5, width = 0.95, height = 0.95))
    lapply(seq_along(object@sequence), function(x) {
        xpos <- coordx[x]
        ypos <- coordy[x]
        textLabel(label_txt = object@sequence[x],
                  x_ = xpos,
                  y_ = ypos,
                  w_ = nuclDims(fontsize = object@cex)[1],
                  h_ = nuclDims(fontsize = object@cex)[2],
                  bg = object@background,
                  cex = object@cex,
                  borderCol = object@border,
                  gp_ = grid::gpar(fontsize = object@cex * 100,
                                   fontface = "bold",
                                   col = object@base_colors[[object@sequence[x]]]))
    })
    grid::popViewport(1)
})
