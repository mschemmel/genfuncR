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
#' textLabel(x_ = 1, y_ = 3, w_ = 1, h_ = 1, label_txt = "Test", angle = 45)

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

#' plots nucleotide sequence as image
#'
#' @param sequence A nucleotide sequence string
#' @param base_colors A list of colors to be used for every nucleotide
#' return Viewport with printed nucleotide sequence
#'
#' @examples
#' dna_to_img("TATCGATCGATC", list(A = "#9EE362", T = "#00C0D0", G = "#FFD403", C = "#FF9356", U = "#d83131"))

dna2img <- setClass("dna2img",
                    slots = list(sequence = "character",
                                 base_colors = "list",
                                 letter_per_line = "numeric"),
                    prototype = list(base_colors = list(A = "#9EE362",
                                                        T = "#00C0D0",
                                                        G = "#FFD403",
                                                        C = "#FF9356",
                                                        U = "#d83131",
                                                        N = "gray60"),
                                     letter_per_line = 50))

dna2img <- function(sequence,
                    base_colors,
                    letter_per_line) {
    .Object <- new("dna2img")
    .Object@sequence <- unlist(strsplit(sequence, split = ""))
    .Object@letter_per_line = letter_per_line
    return(.Object)
}

setMethod(f = "show",
          signature = "dna2img",
          definition = function(object) {

    grid::grid.newpage()

    # get coordinates for every nucleotide
    ypositions <- seq(0.90, 0, -0.025)
    coordx <- rep(0.5:object@letter_per_line, ceiling(length(object@sequence)/object@letter_per_line))[1:length(object@sequence)] / object@letter_per_line
    coordy <- rep(ypositions[1:ceiling(length(object@sequence)/object@letter_per_line)], each = object@letter_per_line)[1:length(object@sequence)]

    # draw all nucleotides
    grid::pushViewport(grid::viewport(x = 0.5, y = 0.5, width = 0.95, height = 0.95))

    lapply(seq_along(object@sequence), function(x){
        xpos <- coordx[x]
        ypos <- coordy[x]
        textLabel(label_txt = object@sequence[x],
                  x_ = xpos,
                  y_ = ypos,
                  w_ = 0.01,
                  h_ = 0.1,
                  gp_ = grid::gpar(fontsize = 12,
                                   fontface = "bold",
                                   col = object@base_colors[[object@sequence[x]]]))
        })
    grid::popViewport(1)
})