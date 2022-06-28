#' plots nucleotide sequence as image
#'
#' @param sequence A nucleotide sequence string
#' @param base_colos A list of colors to be used for every nucleotide
#' return Viewport with printed nucleotide sequence
#'
#' @examples
#' dna_to_img("TATCGATCGATC", list(A = "#9EE362", T = "#00C0D0", G = "#FFD403", C = "#FF9356", U = "#d83131"))

dna_to_img <- function(sequence, base_colors = NULL) {
    seq_ <- unlist(strsplit(sequence, split = ""))

    if(is.null(base_colors)) {
        base_colors <- list(A = "#9EE362", T = "#00C0D0", G = "#FFD403", C = "#FF9356", U = "#d83131")
    }

    # define some constants
    y_ <- 0.95          # starting point
    y_spacer <- 0.05    # spacer in y direction
    x_spacer <- 20      # spacer in x direction

    grid::pushViewport(grid::viewport(x = 0.5, y = 0.5, width = 0.99, height = 0.99))
    count <- 1
    for (i in seq_) {
        bcol <- unlist(unname(base_colors[i]))
        x_ <- count / x_spacer
        if (x_ == 1) {
            count <- 1
            x_ <- count / x_spacer
            y_ <- y_ - y_spacer
        }
        grid::pushViewport(grid::viewport(x = x_, y = y_, width = 0.01, height = 0.1))
        grid::grid.text(i, gp = grid::gpar(fontsize = 12, fontface = "bold", col = bcol))
        grid::popViewport(1)
        count <- count + 1
    }
}
