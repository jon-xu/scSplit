a <- unlist(read.table('barcodes_A.tsv'))
c <- unlist(read.table('barcodes_C.tsv'))
s1 <- unlist(read.table('barcodes_0.csv'))
s2 <- unlist(read.table('barcodes_1.csv'))
if (length(s1) > length(s2)) {
  temp <- s1
  s1 <- s2
  s2 <- temp
}
print (length(intersect(a,s1)) / length(a))
print (length(intersect(c,s2)) / length(c))
