#' puts a text label on a specific position
#' @param vp_name Name of the viewport
#' @param x_ x coordinate
#' @param y_ y coordinate
#' @param w_ w coordinate
#' @param h_ h coordinate
#' @param label_txt Text to be printed on that position
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

#' plot alignment
#' @param mirna Sequence of miRNA sequence
#' @param target Sequence of target sequence
#' @param mirna_name Name of the miRNA
#' @param target_name Name of the target
#' @param alignment_type Type of alignment
#' @param match_color Color of matched nucleotides
#' @param mismatch_color Color of mismatches nucleotides
#' @param highlight_area Area to highlight
#' @param highlight_color Color of highlighted area
#' @param target_position_label Positional label of target sequence
#' @param target_position_label_rot Angle of positional label
#' @example
#' mirnali(mirna, target)

mirnali <- function(mirna,
                    target,
                    mirna_name = NULL,
                    target_name = NULL,
                    alignment_type = "|",
                    match_color = "black",
                    mismatch_color = "black",
                    highlight_area = NULL,
                    highlight_color = "gray90",
                    target_position_label = NULL,
                    target_position_label_rot = 90) {

    mirna <- toupper(unlist(strsplit(gsub("T", "U", mirna), split = "")))
    target <- toupper(unlist(strsplit(gsub("T", "U", target), split = "")))

    # check sequence properties
    if (length(mirna) != length(target)) {
        cat("Error:\n")
        cat("miRNA:\t", length(mirna), "\n")
        cat("Target:\t", length(target), "\n")
        stop("Unequal number of nucleotides. Please check.")
    }
    if (length(mirna) < 15 || length(mirna) > 40) {
        stop("Range of sequence should be between 15 and 40 nucleotides.")
    }
    if (identical(mirna,target)) {
        stop("miRNA and target sequences are equal, wrong orientation?")
    }

    # define and calculate some constants
    cells <- length(mirna) + 2
    x_step <- 1 / cells
    x_start <- x_step / 2
    positions <- seq(x_start, 1, x_step)
    s1_height <- 0.7 # Target
    s2_height <- 0.3 # miRNA
    s1_pos_height <- s1_height + 0.15
    s2_pos_height <- s2_height - 0.15

    # validation
    if (!is.null(target_position_label)) {
        if (length(target_position_label) != length(mirna)) {
            warning("Range of alignment exceeds given target locations.")
            cat("Warning:\n")
            cat("Range of label:", diff(range(target_position_label)) + 1, "\n")
            cat("Range of alignment ", length(mirna), "\n")
        }
    }

    # start plotting grids
    # create newpage to draw on
    grid::grid.newpage()

    # outer viewport
    grid::pushViewport(grid::viewport(name = "outer",
                        x = grid::unit(0.5, "npc"),
                        y = grid::unit(0.5, "npc"),
                        width = 1,
                        height = 1))

    # add sequence name label
    if (!is.null(mirna_name)) {
        text_label(x_ = 0.08,
                   y_ = 0.44,
                   w_ = 0.4,
                   h_ = 0.7,
                   label_txt = mirna_name)
    }
    if (!is.null(target_name)) {
        text_label(x_ = 0.08,
                   y_ = 0.56,
                   w_ = 0.4,
                   h_ = 0.7,
                   label_txt = target_name)
    }

    # main viewport
    grid::pushViewport(grid::viewport(name = "main",
                        x = grid::unit(0.5, "npc"),
                        y = grid::unit(0.5, "npc"),
                        width = 0.7,
                        height = 0.3,
                        just = c("center")))

    # define color and alignment type
    alignment_symbol <- vector()
    nucl_color <- vector()
    for (i in seq(1, length(mirna), 1)) {
        # Six cases:
        # 1. U-A | A-U      complement U-A
        # 2. G-C | C-G      complement G-C
        # 3. U-G | G-U
        # 4. U-C | C-U
        # 5. (ATGC)-(ATGC)   same nucleotide, no difference
        # 6. No complementary

        # 1. U - A | A -U -> true complementary
        if ((mirna[i] == "U"  && target[i] == "A") || (mirna[i] == "A"  && target[i] == "U")) {
            symb <- ifelse(alignment_type == ":", ":", "|")
            colr <- match_color
        }
        # 2. G - C | C -G -> true complementary
        else if ((mirna[i] == "G"  && target[i] == "C") || (mirna[i] == "C"  && target[i] == "G")) {
            symb <- ifelse(alignment_type == ":", ":", "|")
            colr <- match_color
        }
        # 3. U - G | G -U
        else if ((mirna[i] == "U"  && target[i] == "G") || (mirna[i] == "G"  && target[i] == "U")) {
            symb <- ifelse(alignment_type == ":", ".", " ")
            colr <- mismatch_color
        }
        # 4. U - C | C -U
        else if ((mirna[i] == "U"  && target[i] == "C") || (mirna[i] == "C"  && target[i] == "U")) {
            symb <- ifelse(alignment_type == ":", " ", " ")
            colr <- mismatch_color
        }
        # 5. same nucleotide, no difference
        else if (mirna[i] == target[i]) {
            symb <- ifelse(alignment_type == ":", " ", " ")
            colr <- mismatch_color
        }
        # 6. No complementary
        else {
            symb <- ifelse(alignment_type == ":", " ", " ")
            colr <- mismatch_color
        }
        alignment_symbol <- append(alignment_symbol, symb)
        nucl_color <- append(nucl_color, colr)
    }

    # add orientation label
    orientations <- c("5'", "3'")
    orientation_pos <- c(head(positions, 1), tail(positions, 1))
    for (i in c(s1_height, s2_height)) {
        if (i == s2_height) orientations <- rev(orientations)
        for (j in seq(1, 2, 1)) {
            text_label(x_ = orientation_pos[j],
                       y_ = i,
                       w_ = x_step,
                       h_ = 0.15,
                       label_txt = orientations[j],
                       gp_ = grid::gpar(fontsize = 10,
                                        fontface = "bold",
                                        col = "black"))
        }
    }

    # set position to actual sequence position
    positions <- positions[-1]

    # add highlight region
    if (!is.null(highlight_area)) {
        for (i in names(highlight_area)) {
            h_start <- positions[unname(unlist(highlight_area[i]))[1]] - x_start
            h_end <- positions[unname(unlist(highlight_area[i]))[2]] + x_start
            grid::grid.rect(x = h_start,
                            y = 0.5,
                            width = h_end - h_start,
                            height = 0.55,
                            gp = grid::gpar(fill = unname(unlist(highlight_color[i])),
                                            col = "white"),
                            just = c("left"))
        }
    }

    # add sequence label
    for (i in c(s1_height, s2_height)) {
        for (j in seq(1, length(mirna), 1)) {
            text_label(x_ = positions[j],
                       y_ = i,
                       w_ = x_step,
                       h_ = 0.15,
                       label_txt = ifelse(i == s2_height, target[j], mirna[j]),
                       gp_ = grid::gpar(fontsize = 12,
                                        fontface = "bold",
                                        col = nucl_color[j]))
        }
    }

    # add positional label
    if (!is.null(target_position_label)) {
        for (i in c(s1_pos_height, s2_pos_height)) {
            for (j in seq(1, length(mirna), 1)) {
                text_label(x_ = positions[j],
                           y_ = i,
                           w_ = x_step,
                           h_ = 0.15,
                           label_txt = ifelse(i == s1_pos_height,
                                              target_position_label[j],
                                              seq(1, length(mirna), 1)[j]),
                           gp_ = grid::gpar(fontsize = 8, col = "black"),
                           angle = ifelse(i == s1_pos_height,
                                          target_position_label_rot,
                                          0))
                }
        }
    }

    # add alignment label
    for (i in seq(1, length(mirna), 1)) {
        text_label(x_ = positions[i],
                   y_ = 0.5,
                   w_ = x_step,
                   h_ = 0.2,
                   label_txt = alignment_symbol[i],
                   gp_ = grid::gpar(fontsize = 16,
                                    col = "black"))
    }
    grid::popViewport(1)
}
