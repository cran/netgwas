\name{bp}
\alias{bp}
\docType{data}
\title{genotype-phenotype data in mice }

\description{ 
  Data from an intercross between BALB/cJ and CBA/CaJ mouse strains
}

\usage{data(bp)}
\format{
  The format is a matrix containing 93 SNP markers across the genome, and 4 phenotypes: blood pressure (bp), heart rate (hr), body weight (bw), and heart weight (heart-wt), as measured for 163 individuals.
}  

\details{
This genotype data can be used to reconstruct genotype-phenotype networks in mice (see below example) to identify genomic regions that regulate blood pressure, heart rate, and heart weight. 
}

\source{
Sugiyama, F., Churchill, G.A., Li, R., Libby, L.J., Carver, T., Yagami, K.I., John, S.W. and Paigen, B., 2002. QTL associated with blood pressure, heart rate, and heart weight in CBA/CaJ and BALB/cJ mice. Physiological genomics, 10(1), pp.5-12.
}

\examples{
\donttest{
data(bp)
#Constructing genotype-phenotype networks in mice
out <- netphenogeno(bp)
sel <- selectnet(out)
plot(sel, vis= "interactive", vertex.color = c(rep("red",4), rep( "white", 93) ))
}
}

\keyword{datasets}
