\name{plot.select}
\alias{plot.select}

\title{ Plot function for S3 class "select" }

\description{
Plot the optimal graph by model selection  
}

\usage{
\method{plot}{select}(x, vis= NULL, xlab= NULL, ylab= NULL, n.mem= NULL, vertex.label= FALSE
, ..., layout= NULL, label.vertex= "all", vertex.size= NULL, vertex.color= "red",
edge.color= "gray29", sel.label= NULL, label.size = NULL, w.btw= 800, w.within = 10)
}

\arguments{
	\item{x}{An object with S3 class "select" }
	\item{vis}{Visualizing the selected graph using two options: (i) "CI": plotting conditional independence (CI) relationships between variables, and (ii) "interactive": plotting the conditional independence network, where opens a new windows with interactive graph drawing facility. Default is "CI".}
	\item{xlab}{ONLY applicable when \code{vis}= "CI".}
	\item{ylab}{ONLY applicable when \code{vis}= "CI". }
	\item{n.mem}{A vector of memberships. For example, the CviCol dataset, which is provided in the package, contain 5 chromosomes which the total number of markers is \eqn{p = 90}, where the first 24 markers belong into chromosome 1, the next 14 markers into chromosome 2, ..., and chromosome 5 contains 19 markers. Thus, n.mem = c(24,14,17,16,19). If \code{n.mem = NULL} and \code{vis = "CI"} all vertices are coloured the same.}
	\item{vertex.label}{ONLY applicable when \code{vis}= "CI". Assign names to the vertices. Default is FALSE.}
	\item{\dots}{ONLY applicable when \code{vis}= "CI". System reserved (No specific usage)}
	\item{layout}{ONLY applicable when \code{vis}= "interactive". The layout specification. Some graph layouts examples: layout_with_fr, layout_in_circle, layout_as_tree, and layout.fruchterman.reingold. The default layout is \code{layout_with_fr}.}
	\item{label.vertex}{ONLY applicable when \code{vis}= "interactive". Assign names to the vertices. There are three options: "none", "some", "all". Specify "none" to omit vertex labels in the graph; using \code{label.vertex} = "some" you must provide a vector of vertex IDs or a single vertex ID to the \code{sel.label} argument, which you would like to be shown in the graph. Specify "all" to include all vertex labels in the graph. Default is "all".}
	\item{vertex.size}{ Optional. The size of vertices in the graph visualization. The default value is 7.}
	\item{vertex.color}{ONLY applicable when \code{vis}= "interactive". Optional vector (or a color name) giving the colors of the vertices. The default is "red" }
	\item{edge.color}{ONLY applicable when \code{vis}= "interactive". Optional. The default is "gray".}
	\item{sel.label}{ONLY applicable when \code{vis}= "interactive". A vector of vertex IDs or a single vertex ID, which you would like to be shown in the graph. ONLY applicable when \code{label.vertex}="some".}
	\item{label.size}{ONLY applicable when \code{vis}= "interactive". The font size of the vertex labels.}
	\item{w.btw}{Distance between different memberships of \code{n.mem} in layout.}
	\item{w.within}{Distance of nodes within a membership of \code{n.mem} in layout.}
}

\value{
An object with S3 class \code{"select"} is returned:  
 \item{network}{ Plot of a selected graph,  when \code{vis}= "CI".}
 \item{adjacency}{ Conditional independence (CI) relationships between variables, when \code{vis}= "CI"}
 \item{network}{Interactive plot of a selected graph with .eps format, when \code{vis}= "interactive"}
}

\author{ 
  Pariya Behrouzi \cr
  Maintainer: Pariya Behrouzi \email{pariya.behrouzi@gmail.com} \cr
}

\references{
Behrouzi, P., and Wit, E. C. (2017c). netgwas: An R Package for Network-Based Genome-Wide Association Studies. arXiv preprint, arXiv:1710.01236.
}

\seealso{\code{\link{select}}}