
context("mers")

test_that("Reverscomplement oligonucleotide list twice produces original", {
  
  for(oli_length in c(1:4)) {
    
    oligos <- oligo_names(oli_length)
    revcomp_order <- complementary_oligo_positions(oli_length)
    
    expect_equal(revcomp_order[revcomp_order] , 1:length(revcomp_order))
    expect_equal(oligos[revcomp_order][revcomp_order] , oligos)
  }
})


