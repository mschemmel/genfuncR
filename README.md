# genfuncR
A home for R functions I use myself to visualize genomic features. Will grow over time. Should serve as a starting point for own representations and compilations (multi panel plots).

## gene_arrow
Returns a `grid.polygon` object. Can be combined to a series of genes. 
<p align="center">
<img src="./img/genearrow.png" width="50%">
</p>

## dna_to_img
Returns a colored DNA sequence for further editing.
<br>
Example:
```r
library(grid)
source("/PATH/OF/R_FILE.r")

# create random dna sequence
seq <- paste0(sample(c("A", "T", "G", "C"), 76, replace = TRUE), collapse = "")
seq <- unlist(strsplit(seq, split = ""))

# set colors
base_col <- list(A = "#9EE362", T = "#00C0D0", G = "#FFD403", C = "#FF9356")

# print colored dna sequences
dna_to_img(seq, base_col)
```

<p align="center">
<img src="./img/dna_to_img.png" width="60%">
</p>

## mirnali
```r
library(grid)
library(GetoptLong)

mirna <-  "AUUCUCAGCGCCUGUAUAGAGG"
target <- "UAGGAGUCGCGGAAUUAACACU"

# basic example
mirnali(mirna, target)
```

<p align="center">
<img src="./img/mirnali_basic.png" width="60%">
</p>

```r
# more advanced example
mirnali(mirna,
        target,
        mirna_name = "miRNA 1",
        target_name = "mRNA 1",
        alignment_type = ":",
        match_color = "forestgreen",
        mismatch_color = "firebrick",
        highlight_area = list("A" = c(5, 10), "B" = c(13, 17), "C" = c(20, 20)),
        highlight_color = list("A" = "gray90", "B" = "navajowhite", "C" = "gray90"),
        target_position_label = c(200:(199 + nchar(target))),
        target_position_label_rot = 90)
```

<p align="center">
<img src="./img/mirnali_adv.png" width="60%">
</p>