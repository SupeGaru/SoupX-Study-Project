# This R function is designed to remove the background contamination from a count matrix.
#  Given an estimated level of background contamination for a channel, it calculates the resulting corrected 
#  count matrix with background contamination removed. The function subtracts the mean expected background counts for
#   each gene, and then redistributes any "unused" counts. A count is unused if its subtraction has no effect. For example, 
#   subtracting a count from a gene that has zero counts to begin with. The function requires clustering information to provide
#    an accurate removal of contaminating counts, either by setting clustering on the SoupChannel object or explicitly passing 
#    the clusters parameter. The function has three methods for removing the counts, which should almost always be left at the 
#    default ('subtraction'), which iteratively subtracts counts from all genes as described above. The 'soupOnly' method uses
#     a p-value based estimation procedure to identify those genes that can be confidently identified as having endogenous 
#     expression and removes everything else. The 'multinomial' method explicitly maximizes the multinomial likelihood for 
#     each cell. The function also has parameters that control the level of verbosity, the allowed deviation from the expected
#      number of soup counts, and the p-value cut-off used when the method is 'soupOnly'.