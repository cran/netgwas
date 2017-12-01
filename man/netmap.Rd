\name{netmap}
\alias{netmap}

\title{Constructing linkage map for diploids and polyploids}

\description{
This is one of the main functions of \pkg{netgwas} package. This function reconstructs linkage maps for biparental diploid and polyploid organisms using three methods.
}

\usage{
netmap(data, method = "gibbs", rho = NULL, n.rho = NULL, rho.ratio = NULL, 
       cross.typ=c("inbred", "outbred"), iso.m= NULL, use.comu= FALSE,
	   ncores = "all", em.iter = 5, verbose = TRUE)
}

\arguments{
	\item{data}{ An (\eqn{n \times p}) matrix or a \code{data.frame} corresponding to a genotype data matrix (\eqn{n} is the sample size and \eqn{p} is the number of variables). Input data can contain missing values.}
	\item{method}{Three available methods to construct linkage map: "gibbs", "approx", and "npn". For a medium (~500) and a large number of variables we recommend to choose "gibbs" and "approx", respectively. Choosing "npn" for a very large number of variables (> 2000) is computationally efficient. The default method is "gibbs".}
	\item{rho}{ A decreasing sequence of non-negative numbers that control the sparsity level. Leaving the input as \code{rho = NULL}, the program automatically computes a sequence of \code{rho} based on \code{n.rho} and \code{rho.ratio}. Users can also supply a decreasing sequence values to override this.}
	\item{n.rho}{ The number of regularization parameters. The default value is \code{6}.}
	\item{rho.ratio}{ Determines distance between the elements of \code{rho} sequence. A small value of \code{rho.ratio} results in a large distance between the elements of \code{rho} sequence. And a large value of \code{rho.ratio} results into a small distance between elements of \code{rho}. If keep it as NULL the program internally chooses a value.}
	\item{cross.typ}{To be specified either "inbred" or "outbred".}
	\item{iso.m}{ Optional. Number of isolated markers.}
	\item{use.comu}{Use community detection algorithm to detect linkage groups. Default is FALSE.}
	\item{ncores}{ The number of cores to use for the calculations. Using \code{ncores = "all"} automatically detects number of available cores and runs the computations in parallel on (available cores - 1).}
	\item{em.iter}{The number of EM iterations. The default value is 5.}
	\item{verbose}{Providing a detail message for tracing output. The default value is \code{TRUE}.}
}

\details{
Constructing linkage maps for diploid and polyploid organisms. Diploid organisms contain two sets of chromosomes, one from each parent, whereas polyploids contain more than two sets of chromosomes. Inbreeding is mating between two parental lines where they have recent common biological ancestors. If they have no common ancestors up to roughly e.g. 4-6 generations, this is called outcrossing. In both cases the genomes of the derived progenies are random mosaics of the genome of the parents. However, in the case of inbreeding parental alleles are distinguishable in the genome of the progeny; in outcrossing this does not hold.
}

\value{
An object with S3 class \code{"netgwas"} is returned:  
 \item{map}{Constructed linkage map.}
 \item{opt.index}{ The index of selected graph using model selection.}
 \item{allres}{ A list containing results for different regularization parameter. Belongs to class "netgwas". To visualize a path of different 3D maps consider function \code{\link{plot.netgwas}}. Note that the input data is reordered based on the estimated linkage map and is saved as data in this argument.}
}

\author{ 
  Pariya Behrouzi and Ernst C. Wit \cr
  Maintainers: Pariya Behrouzi \email{pariya.behrouzi@gmail.com}\cr
}

\references{
1. Behrouzi, Pariya, and Ernst C. Wit. "De novo construction of q-ploid linkage maps using discrete graphical models." arXiv preprint arXiv:1710.01063 (2017). \cr
2. Behrouzi, Pariya, and Ernst C. Wit. "netgwas: An R Package for Network-Based Genome-Wide Association Studies." arXiv preprint arXiv:1710.01236 (2017). \cr
3. Guo, Jian, Elizaveta Levina, George Michailidis, and Ji Zhu. "Graphical models for ordinal data." Journal of Computational and Graphical Statistics 24, no. 1 (2015): 183-204. \cr
4. Liu, Han, Fang Han, Ming Yuan, John Lafferty, and Larry Wasserman. "High-dimensional semiparametric Gaussian copula graphical models." The Annals of Statistics 40, no. 4 (2012): 2293-2326. \cr 
5. Witten, Daniela M., Jerome H. Friedman, and Noah Simon. "New insights and faster computations for the graphical lasso." Journal of Computational and Graphical Statistics 20, no. 4 (2011): 892-900.
}

\examples{
\dontrun{
data(CviCol)
#Randomly change the order of markers across the genome
cvicol <- CviCol[ ,sample(1:ncol(CviCol), ncol(CviCol), replace=FALSE)]

#Constructing linkage map using gibbs method
out <- netmap(cvicol, cross.typ= "inbred", vis= FALSE, ncores=1); out
#Estimated linkage map
map <- out$map; map
#Plot the associated three-dimensional map
plot(out)
#Visualize other three-dimensional maps
plot(out$allres)
#Build a linkage map for 5th three-dimensional map
bm <- buildMap(out, opt.index=5); bm
####################

#Constructing linkage map using approx method
out2 <- netmap(cvicol, method="approx", cross.typ= "inbred", vis= FALSE, ncores=1); out2
#Estimated linkage map
map2 <- out2$map; map2
#Plot the associated 3-dimensional map
plot(out2)
#Visualize other three-dimensional maps
plot(out2$allres)
#Build a linkage map for 5th three-dimensional map
bm2 <- buildMap(out2, opt.index=5); bm2

#Constructing linkage map using npn method
out3 <- netmap(cvicol, method="npn", cross.typ= "inbred", vis= FALSE, ncores=1); out2
#Estimated linkage map
map3 <- out3$map; map3
#Plot the associated 3-dimensional map
plot(out3)
}
}
