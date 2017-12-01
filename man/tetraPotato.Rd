\name{tetraPotato}
\alias{tetraPotato}
\docType{data}
\title{tetraploid potato genotype data }

\description{ 
  Tetraploid potato (\emph{Solanum tuberosum L.}) genotype data.
}

\usage{data(tetraPotato)}
\format{
  The format is a matrix containing 1972 single-nucleotide polymorphism (SNP) markers for 156 individuals.
}  

\details{
The full-sib mapping population MSL603 consists of 156 F1 plants resulting from a cross
between female parent "Jacqueline Lee" and male parent "MSG227-2". The obtained genotype
data contain 1972 SNP markers with five allele dosages. This genotype data can be used to construct linkage map for tetraploid potato (see below example).
}

\author{ 
  Pariya Behrouzi and Ernst C. Wit \cr
  Maintainer: Pariya Behrouzi \email{pariya.behrouzi@gmail.com}\cr
}

\references{
1. Behrouzi, P., and Wit, E. C. (2017c). netgwas: An R Package for Network-Based Genome-Wide Association Studies. arXiv preprint, arXiv:1710.01236. \cr
2. Massa, Alicia N., Norma C. Manrique-Carpintero, Joseph J. Coombs, Daniel G. Zarka, Anne E. Boone, William W. Kirk, Christine A. Hackett, Glenn J. Bryan, and David S. Douches. "Genetic linkage mapping of economically important traits in cultivated tetraploid potato (Solanum tuberosum L.)." G3: Genes, Genomes, Genetics 5, no. 11 (2015): 2357-2364.
}

\examples{
\donttest{
data(tetraPotato)
#Shuffle the order of markers
potato <- tetraPotato[ ,sample(1:ncol(tetraPotato), ncol(tetraPotato), replace=FALSE)] 
#Constructing linkage map for tetraploid potato
out <- netmap(potato, method="npn", cross.typ= "outbred"); out
potato.map <- out$map; potato.map
plot(out)
}
}

\keyword{datasets}
