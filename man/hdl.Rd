\name{hdl}
\alias{hdl}
\docType{data}
\title{ Mus Musculus HDL data in mice}
\description{ 
 HDL QTL data was obtained from a F2 inner-cross between inbred MRL/MpJ and SM/J strains of mice. 
}

\usage{data(hdl)}
\format{
The format is a matrix consists of 280 observations for 15 variables: genotype data (genotype
states at 5 SNP markers) and phenotype data (HDL levels and normalized expression values of 10
genes). Three possible genotype states MM (homozygous) are denoted by 1, H (heterozygous) by 2, and SS (homozygous) by 3 and phenotypes are of class numeric.
}  

\details{
The \emph{Mus Musculus} HDL data were obtained from an F2 inner-cross between inbred MRL/MpJ and SM/J strains of mice.
}

\source{
Leduc MS, Blair RH, Verdugo RA, Tsaih SW, Walsh K, Churchill GA, Paigen B.(2012). ``Using
bioinformatics and systems genetics to dissect HDL-cholesterol genetics in an MRL/MpJ x SM/J
intercross." J Lipid Res., 6, 1163-75.
}

\examples{
\donttest{
data(hdl)
#Constructing genotype-phenotype networks in mice
out <- netphenogeno(hdl)
sel <- selectnet(out)
plot(sel, vis= "CI", vertex.label= TRUE, xlab= "Genotype-Phenotype", ylab= "Genotype-Phenotype")
cl <- c(rep("white", 5), rep("red", 10))
plot(sel, vis="interactive", n.mem=c(5, 10), vertex.color= cl, w.btw= 10, w.within = 8  )

}
}

\keyword{datasets}