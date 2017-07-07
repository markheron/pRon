<!-- README.md is generated from README.Rmd. Please edit that file -->

pretty Representations of oligo-nucleotides
-------------------------------------------

The **pRon** package contains functions to help handle oligo-nucleotides (k-mers), including creating frequency/energy profile figures.

### Summary

This package provides some useful functions that I frequently used:

-   quick summary of large objects in memory
-   ruler style axis for plots
-   nicer (and more convinient) heatmap
-   smearing and running mean functions for fast smoothing

### Installation

``` r
devtools::install_bitbucket(repo="markheron/pRon", auth_user="user_name", password="your_password")
```

or install all of Marks packages with the **maRs** package:

``` r
install_marks_package_suite(auth_user="user_name")
```

### Usage

``` r
library(pRon)

oligo_names(2)
#>  [1] "AA" "AC" "AG" "AT" "CA" "CC" "CG" "CT" "GA" "GC" "GG" "GT" "TA" "TC"
#> [15] "TG" "TT"
```

``` r

fastas <- cut_out_fasta_multiple(chr=chr_vec, pos = pos_vec, strand = strand_vec, size = 100, order = 2, genome_dir = "path/to/genome")
num_fasta_matrix <- fasta2num(fastas=fastas, oligo_length = 2)
plotOligoFreqs(num_fasta_matrix)
```

### License

The **pRon** package is licensed under the GPLv3 (<http://www.gnu.org/licenses/gpl.html>).
