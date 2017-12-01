\name{plot.netgwasmap}
\alias{plot.netgwasmap}

\title{ plot for S3 class "netgwasmap" }

\description{
Plot the graph associated with constructed linkage map via function \code{\link{netmap}}.
}

\usage{
\method{plot}{netgwasmap}( x,  layout, vertex.size, vertex.label, label.size, vertex.color, ... )
}

\arguments{
	\item{x}{An object from "netgwasmap" class.}
	\item{layout}{The vertex placement algorithm which is according to \pkg{igraph} package}
	\item{vertex.size}{ Optional integer to adjust vertex size in graph G. Default is 5.}
	\item{vertex.label}{To show the label of vertices. Default is FALSE.}
	\item{label.size}{Optional integer to adjust the size of node's label in graph G. Applicable when vertex.label is TRUE. Default is 0.8.}
	\item{vertex.color}{Optional integer vectors giving colors to the vertices. }
	\item{\dots}{System reserved (no specific usage)}
}

\author{ 
  Pariya Behrouzi and Ernst C. Wit \cr
  Maintainer: Pariya Behrouzi \email{pariya.behrouzi@gmail.com}
}

\references{
1. Behrouzi, P., and Wit, E. C. (2017b). De novo construction of q-ploid linkage maps using discrete graphical models. arXiv preprint, arXiv:1710.01063. \cr
2. Behrouzi, P., and Wit, E. C. (2017c). netgwas: An R Package for Network-Based Genome-Wide Association Studies. arXiv preprint, arXiv:1710.01236.
}


\seealso{\code{\link{netmap}}.}

