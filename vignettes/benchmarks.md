<!--
%\VignetteIndexEntry{pRon-benchmark}
%\VignetteEngine{knitr::knitr}
-->

Benchmarks: pROn
========================================================


```
## Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## 
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## 
## The following object is masked from 'package:stats':
## 
##     xtabs
## 
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, as.data.frame, cbind, colnames, duplicated,
##     eval, Filter, Find, get, intersect, lapply, Map, mapply,
##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
##     Position, rank, rbind, Reduce, rep.int, rownames, sapply,
##     setdiff, sort, table, tapply, union, unique, unlist
## 
## Loading required package: IRanges
```



Several benchmarks to check which way of implementing the functions is faster.


oligo_names
--------------------------------------------------------

```r
oligo_nuc_names <- c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", 
    "GG", "GT", "TA", "TC", "TG", "TT")
tmp <- microbenchmark(LAPPLY = unlist(lapply(oligo_nuc_names, function(x) paste0(x, 
    c("A", "C", "G", "T")))), SAPPLY = sapply(oligo_nuc_names, function(x) paste0(x, 
    c("A", "C", "G", "T"))), REP = paste0(rep(oligo_nuc_names, each = 4), c("A", 
    "C", "G", "T")), times = 1000)
boxplot(tmp, unit = "ms", log = TRUE, main = "oligo_names")
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2.png) 


fasta2num
--------------------------------------------------------

```r

base2num <- c(1:4, NA)
names(base2num) <- c("A", "C", "G", "T", "N")

fasta_seqs <- DNAStringSet(rep(paste0(sample(c("A", "C", "G", "T"), replace = TRUE, 
    3e+05), collapse = ""), 10))
```

```
## Note: no visible binding for global variable 'NAMES' 
## Note: no visible binding for global variable 'NAMES'
```

```r
tmp <- microbenchmark(lFISV_rowSums_vecmult = {
    tmp <- sapply(fasta_seqs, function(x) rowSums(letterFrequencyInSlidingView(x, 
        1, c("A", "C", "G", "T")) %*% (1:4)))
    tmp[tmp == 0] <- NA
}, lFISV_colSums_t = {
    tmp <- sapply(fasta_seqs, function(x) colSums(t(letterFrequencyInSlidingView(x, 
        1, c("A", "C", "G", "T"))) * (1:4)))
    tmp[tmp == 0] <- NA
}, strsplit = sapply(fasta_seqs, function(fasta) base2num[strsplit(as.character(fasta), 
    "")[[1]]]), times = 10)
```

```
## Note: no visible binding for global variable 'xp' 
## Note: no visible binding for global variable '.link_to_cached_object'
```

```r
boxplot(tmp, unit = "ms", log = TRUE, main = "fasta2num few long seqs")
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-31.png) 

```r

fasta_seqs <- DNAStringSet(rep(paste0(sample(c("A", "C", "G", "T"), replace = TRUE, 
    3000), collapse = ""), 10))
tmp <- microbenchmark(lFISV_rowSums_vecmult = {
    tmp <- sapply(fasta_seqs, function(x) rowSums(letterFrequencyInSlidingView(x, 
        1, c("A", "C", "G", "T")) %*% (1:4)))
    tmp[tmp == 0] <- NA
}, lFISV_colSums_t = {
    tmp <- sapply(fasta_seqs, function(x) colSums(t(letterFrequencyInSlidingView(x, 
        1, c("A", "C", "G", "T"))) * (1:4)))
    tmp[tmp == 0] <- NA
}, strsplit = sapply(fasta_seqs, function(fasta) base2num[strsplit(as.character(fasta), 
    "")[[1]]]), times = 1000)
boxplot(tmp, unit = "ms", log = TRUE, main = "fasta2num few short seqs")
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-32.png) 

```r

fasta_seqs <- DNAStringSet(rep(paste0(sample(c("A", "C", "G", "T"), replace = TRUE, 
    3000), collapse = ""), 1000))
tmp <- microbenchmark(lFISV_rowSums_vecmult = {
    tmp <- sapply(fasta_seqs, function(x) rowSums(letterFrequencyInSlidingView(x, 
        1, c("A", "C", "G", "T")) %*% (1:4)))
    tmp[tmp == 0] <- NA
}, lFISV_colSums_t = {
    tmp <- sapply(fasta_seqs, function(x) colSums(t(letterFrequencyInSlidingView(x, 
        1, c("A", "C", "G", "T"))) * (1:4)))
    tmp[tmp == 0] <- NA
}, strsplit = sapply(fasta_seqs, function(fasta) base2num[strsplit(as.character(fasta), 
    "")[[1]]]), times = 10)
boxplot(tmp, unit = "ms", log = TRUE, main = "fasta2num many short seqs")
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-33.png) 


