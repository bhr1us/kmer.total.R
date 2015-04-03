#remove header from fasta file

fasta2raw <- function(input_file, output_file) {
  library(seqinr)
  x <- read.fasta(file = input_file, seqonly=TRUE, strip.desc=TRUE)
  y <- toString(x)
  z <- gsub(",", "", y)
  zz <- gsub(" ", "", z)
  write(zz, file = output_file)
  out <- paste(output_file, "has been saved.", collapse=" ")
  print(out)
}