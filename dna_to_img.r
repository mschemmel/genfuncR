dna_to_img <- function(sequence, base_colors) {
    # define some constants
    y_ <- 0.95          # starting point
    y_spacer <- 0.05    # spacer in y direction
    x_spacer <- 20      # spacer in x direction

    grid::pushViewport(grid::viewport(x = 0.5, y = 0.5, width = 0.99, height = 0.99))
    count <- 1
    for (i in sequence) {
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
