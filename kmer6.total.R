#Compare 2 lists of kmers

options(warn=-1)

#define function
kmer.total <- function (kmerfile1, kmerfile2, output.file) {

cat ("Importing files...", "\n")
  
# read kmer files 
kmer1 <- read.table(kmerfile1)
kmer2 <- read.table(kmerfile2)

cat ("Selecting k-mers...", "\n")

#Select unique kmers
kmer1.single <- kmer1[kmer1$V2 == 1,]
kmer2.single <- kmer2[kmer2$V2 == 1,]

#Convert to vectors
kmer1.vec <- as.vector (kmer1.single$V1)
kmer2.vec <- as.vector (kmer2.single$V1)

cat ("Extracting unique k-mers...", "\n")

#Omit kmers with polynucleotides greater than 6
remove.A.1 <- grep ("A{6,}", kmer1.vec, value=TRUE)
remove.T.1 <- grep ("T{6,}", kmer1.vec, value=TRUE)
remove.C.1 <- grep ("C{6,}", kmer1.vec, value=TRUE)
remove.G.1 <- grep ("G{6,}", kmer1.vec, value=TRUE)
remove.kmer1 <- c(remove.A.1, remove.T.1, remove.C.1, remove.G.1)
kmer1.vec.uni <- unique(remove.kmer1)
kmer1.vec.new <- kmer1.vec[kmer1.vec != kmer1.vec.uni]

remove.A.2 <- grep ("A{6,}", kmer2.vec, value=TRUE)
remove.T.2 <- grep ("T{6,}", kmer2.vec, value=TRUE)
remove.C.2 <- grep ("C{6,}", kmer2.vec, value=TRUE)
remove.G.2 <- grep ("G{6,}", kmer2.vec, value=TRUE)
remove.kmer2 <- c(remove.A.2, remove.T.2, remove.C.2, remove.G.2)
kmer2.vec.uni <- unique(remove.kmer2)
kmer2.vec.new <- kmer2.vec[kmer2.vec != kmer2.vec.uni]

cat ("Counting k-mers...", "\n")

#Calculate number of kmers
kmer1.len <- length (kmer1.vec.new)
kmer2.len <- length (kmer2.vec.new)
mean.vec <- c(kmer1.len, kmer2.len)
mean.length <- mean(mean.vec)

#count shared kmers
total <- length (intersect (kmer1.vec.new, kmer2.vec.new))

#Calculate basic stats
per.shared <- (total/mean.length)*100
num.diff <- mean.length - total
ratio.dif <- num.diff/mean.length

#Print Summary
Print.L <- c("KMER ANALYSIS SUMMARY","======================", "Input Files:", 
"   kmer1 file ", kmerfile1, "      length ", kmer1.len, "   kmer2 file ", kmerfile2,
"      length ", kmer2.len,  "Analysis:", "   Mean kmer length ", mean.length,
"   # shared kmers ", total, "   % kmers shared ", per.shared,"   # kmers different ", num.diff,
"   Ratio of differences ", ratio.dif)

write (Print.L, file = output.file)

cat ("Output saved as: ", output.file)
}