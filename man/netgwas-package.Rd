\name{netgwas-package}
\alias{netgwas-package}
\docType{package}
\title{
  Network Based Genome Wide Association Studies
}

\description{
The \code{R} package \pkg{netgwas} provides a set of tools based
		on undirected graphical models for accomplishing three important
		and interrelated goals in genetics: (1) linkage map construction, 
		(2) reconstructing intra- and inter-chromosomal conditional
		interactions (linkage disequilibrium) networks, and (3) exploring 
		high-dimensional genotype-phenotype network and genotype-phenotype- 
		environment interactions network. The netgwas can deal with biparental 
		species with any  ploidy level.
		The package implemented the recent improvements both for construction 
		of linkage maps in diploid and polyploid species in Behrouzi and Wit(2017b),
		and in reconstructing networks for non-Gaussian data, ordinal data, and 
		mixed continuous and discrete data in Behrouzi and Wit (2017a). One 
		application is to uncover epistatic interactions network, where the network
		captures the conditionally dependent short- and long-range linkage disequilibrium
		structure of a genomes and reveals aberrant marker-marker associations.
		In addition, Behrouzi and Wit(2017c) implemented their proposed method to explore 
		genotype-phenotype networks where nodes are either phenotypes or genotypes, and each
		phenotype is connected by an edge to a genotype or a group of genotypes if 
		there is a direct association between them, given the rest of the variables. 
		Different phenotypes may also interconnect. The conditionally dependent relationships
		between markers on a genome and phenotypes is determined through Gaussian copula graphical model. 
		We remark that environmental variables can also be included along with genotype-phenotype
		input data to reconstruct networks between genotypes, phenotypes, and environment 
		variables. Beside, the package contains functions for simulation and visualization,
		as well as three multivariate datasets taken from literature.
}

\author{Pariya Behrouzi and Ernst C. Wit \cr
Maintainers: Pariya Behrouzi \email{pariya.behrouzi@gmail.com} \cr
}

\references{
1. Behrouzi, P. and Wit, E.C., (2017a). Detecting epistatic selection with partially observed genotype data by using copula graphical models. Journal of the Royal Statistical Society: Series C (Applied Statistics).\cr
2. Behrouzi, P., and Wit, E. C. (2017b). De novo construction of q-ploid linkage maps using discrete graphical models. arXiv preprint, arXiv:1710.01063. \cr
3. Behrouzi, P., and Wit, E. C. (2017c). netgwas: An R Package for Network-Based Genome-Wide Association Studies. arXiv preprint, arXiv:1710.01236.
}

\examples{
\dontrun{
# Notice: first install 'RBGL' from Bioconductor for the proper functionality of 'netgwas'
source("https://bioconductor.org/biocLite.R")
biocLite("RBGL")


install.packages("netgwas")
library(netgwas)
}
}
\keyword{package}