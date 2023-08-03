#' read fasta file and return list of entries
#' @param fastapath filepath of fasta file
#' @examples
#' fasta(test.fa)
fasta <- function(fastapath) {
    fastafile <- readLines(fastapath)
    # convert input to vector of header (0) and sequence lines (1) -> binary
    binary <- unname(sapply(fastafile, function(x) ifelse(substr(x, 1, 1) == ">", 0, 1)))

    # get groups of binary vector to classify individual sequences
    groups <- cumsum(diff(c(0, binary)) > 0) 
    groups[binary == 0] <- 0

    # get header
    header <- fastafile[which(groups == 0)]

    # concatenate sequences of individual lines to single string
    sequences <- lapply(seq_along(1:max(groups)),
                        function(x) {
                            paste0(fastafile[which(groups == x)], collapse = "")
                        })

    # name list of sequences with header
    names(sequences) <- header
    return(sequences)
}
