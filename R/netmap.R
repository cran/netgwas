#-------------------------------------------------------------------------------#
# Package: Network-Based Genome-Wide Association Studies                        #
# netmap(): Reconstruct conditional dependence networks among markers on genome #
# Authors: Pariya Behrouzi, Ernst Wit                                           #
# maintainer: <pariya.Behrouzi@gmail.com>                                       #
# Date: Nov 21th 2017                                                           #
# Version: 0.0.1-1                                                              #
#-------------------------------------------------------------------------------#

netmap = function(data, method = "gibbs", rho = NULL, n.rho = NULL, rho.ratio = NULL, cross.typ=c("inbred", "outbred"), iso.m= NULL, use.comu= FALSE, ncores = "all", em.iter = 5, verbose = TRUE) 
{
	if(is.null(cross.typ)) stop("Please fill \"inbred\" or \"outbred\" as a cross-type! \n")
	if(!is.matrix(data)) data <- as.matrix(data)
	if(ncores == "all") ncores <- detectCores() - 1
	if(is.null(em.iter)) em.iter = 5
	em.tol = 0.001
	
	
	n = nrow(data)
	p = ncol(data)
	result = list()
	data <- cleaning.dat(data)
	
if( method == "gibbs" ||  method== "approx" ) 
{
	if( method == "gibbs")
	{
		if((is.null(rho)) && (is.null(n.rho)) ) n.rho = 6
		if(! is.null(rho)) n.rho  = length(rho)
		if(is.null(rho.ratio)) rho.ratio = 0.45
		
		est = vector("list", n.rho)
		for(chain in 1 : n.rho ) 
		{
			if(verbose)
			{
				m <- paste(c("Constructing linkage map is in progress:", floor(100 * chain/n.rho), "%"), collapse="")
				cat(m, "\r")
				flush.console()
			}
			if( chain == 1)
			{
				est[[chain]] = vector("list", n.rho)
				Theta = sparseMatrix(i = 1:ncol(data), j = 1:ncol(data), x = 1)
				est[[chain]] = Gibbs_method(data, rho=rho, n_rho=n.rho, rho_ratio=rho.ratio, Theta = Theta, ncores = ncores, chain = chain, max.elongation = em.iter, em.tol=em.tol)
			}else{
				est[[chain]] = vector("list", n.rho)
				Theta = est[[(chain - 1)]]$Theta
				Theta = as(Theta, "dgTMatrix") 
				Theta = as(Theta, "sparseMatrix")
				est[[chain]] = Gibbs_method(data, rho=rho, n_rho=n.rho, rho_ratio=rho.ratio, Theta= Theta, ncores = ncores, chain = chain, max.elongation = em.iter, em.tol=em.tol)
			}
		}
		rm(Theta)
		gc()
	}
	
	if(method == "approx")
	{

		if( !is.null(rho) ) n.rho = length(rho) 
		if( is.null(n.rho) ) n.rho = 6
		if(is.null(rho.ratio)) rho.ratio = 0.65
		
		ini = initialize(data, rho = rho, n_rho = n.rho, rho_ratio = rho.ratio, ncores=ncores )
		rho = ini$rho
		Z	= ini$Z
		ES	= ini$ES
		lower_upper = ini$lower_upper
		
		rm(ini)
		gc()
		
		est <- vector("list", n.rho)
		for(chain in 1 : n.rho) 
		{
		if(verbose)
			{
				m <- paste(c("Constructing linkage map is in progress:", floor(100 * chain/n.rho), "%"), collapse="")
				cat(m, "\r")
				flush.console()
			}
			est[[chain]] <- vector("list", n.rho)
			est[[(chain)]] <- approx_method(data, Z, ES=ES, rho=rho, lower_upper=lower_upper, chain = chain, ncores = ncores, em_tol=em.tol, em_iter=em.iter)
		}

		rm(lower_upper)
		gc()
	}

	result$Theta  = vector("list", n.rho)
	result$path   = vector("list", n.rho)
	result$ES	  = vector("list", n.rho)
	result$Z	  = vector("list", n.rho)
	result$rho	  = vector()
	result$loglik = vector()
	result$data	  = data
	rm(data)
	
	for(chain in 1:n.rho)
	{
		result$Theta[[chain]]	= est[[chain]]$Theta
		if(!is.null(colnames(result$data))) colnames(result$Theta[[chain]]) = colnames(result$data)
		result$path[[chain]]	= abs(sign(result$Theta[[chain]])) - Diagonal(p)
		result$rho[chain]		= est[[chain]]$rho
		result$loglik[chain]	= est[[chain]]$loglik
	}
	rm(est)
	class(result) = "netgwas"
}else{	
	if(method == "npn")
	{
		if(verbose)
		{
			m <- paste("Constructing linkage map is in progress ... \n")
			cat(m, "\r")
			flush.console()
		}
		if((is.null(rho)) && (is.null(n.rho)) ) n.rho = 6
		if(! is.null(rho)) n.rho  = length(rho)
		if(is.null(rho.ratio)) rho.ratio = 0.45	
		if( any(is.na(data)) ) npn.func = "shrinkage" else npn.func = "skeptic"
		tdata <- npn(data, npn.func= npn.func)		
		est <- huge(tdata, lambda= rho, nlambda=n.rho, lambda.min.ratio=rho.ratio, method="glasso", verbose=FALSE)
			
		result$Theta  = est$icov
		result$path   = est$path
		result$rho	  = est$lambda
		result$loglik = n/2 * (est$loglik)
		result$data	  = data
		
		rm(data, est)		
		if(is.null(colnames(result$Theta[[1]]))) {for(i in 1:length(result$rho)) colnames(result$Theta[[i]]) <- colnames(result$data)}
		if(is.null(colnames(result$path[[1]]))) {for(i in 1:length(result$rho)) colnames(result$path[[i]]) <- colnames(result$data)}
		class(result) = "netgwas"
	}	
}	
	#Detecting LGs and ordering within each LG
	if(is.null(use.comu)) use.comu <- FALSE
	sel.net <- selectnet(result, criteria = "ebic", ebic.gamma = 0.5, ncores = ncores, verbose=FALSE)
	if(method == "npn") colnames(sel.net$opt.theta) <- colnames(result$data)
	map <- buildMap.internal(sel.net$opt.theta, cross= cross.typ, num.iso.m= iso.m, use.comu=use.comu)

	result$data <- result$data[ , c(as.character(map[,1]))]
	results <- list( map= map, opt.index= sel.net$opt.index, cross.typ=cross.typ, allres=result)
	class(results$allres) = "netgwas"
	class(results) = "netgwasmap"
	rm(result)
	
	cat("Constructing linkage map is done.             \r\n")
	return(results)
}

#-----------------------------------------------------#
#   		Plot for class "netgwasmap"      	      #
#-----------------------------------------------------#
plot.netgwasmap = function( x, layout=NULL, vertex.size=NULL, vertex.label= NULL, label.size= NULL, vertex.color=NULL, ...)
{
	if(class(x) != "netgwasmap") stop("netgwas.object should belong to the netgwas class. \n ")
	if(is.null( vertex.color))  vertex.color <- "red"
	if(is.null( label.size)) label.size <- 0.8
	if(is.null( vertex.size)) vertex.size <- 5
	if(is.null( vertex.label)) vertex.label <- FALSE

	opt.theta <- as.matrix(x$allres$Theta[[x$opt.index]]) 
	p <- ncol(opt.theta)
	path <-  as.matrix(abs(sign(opt.theta)) - diag(rep(1,p)))
	adj <- graph.adjacency(path, mode="undirected")
	adj$label.cex <- label.size
	if(vertex.label == TRUE) {vertex.label <- colnames(opt.theta)}else{vertex.label <- NA}
	if(is.null(layout)) layout <- layout.fruchterman.reingold

	split.screen( figs = c( 1, 2 ) )
	split.screen( figs = c( 1, 1 ), screen = 1 )
	split.screen( figs = c( 2, 1 ), screen = 2 )
	screen(1)
	plot(adj, layout=layout, edge.curved = F, vertex.label= vertex.label, vertex.color=vertex.color, edge.color="gray40", vertex.size=vertex.size, vertex.label.dist=0, vertex.label.color="darkblue", main="Three-dimensional map")
	screen(4)
	image(path, col = gray.colors(256), xlab= "markers", ylab="markers" ,main= "Conditional dependence relationships \nbefore ordering", cex.main=.8, cex.lab=.8, cex.axis=.8)

	index <- as.character(x$map[ ,1])
	rownames(path) <- colnames(path)
	path.After <- path[c(index), c(index)] 
	screen(5)
	image(path.After, col = gray.colors(256), xlab= "markers", ylab="markers" ,main="Conditional dependence relationships \nafter ordering", cex.main=0.8, cex.lab=.8, cex.axis=.8)
}

#-----------------------------------------------------#
#   		Summary for class "netgwasmap"            #
#-----------------------------------------------------#
print.netgwasmap = function(x, ...){
	cat("Number of linkage groups: ", length(unique(sort(x$map[,2]))), "\n")
	cat("Number of markers per linkage group: ", table(x$map[, 2]), "\n")
	cat("Total number of markers in the linkage map:", length(x$map[,1]),".", " (", (ncol(x$allres$Theta[[1]]) - length(x$map[,1])), " markers removed from the input genotype data) \n")
	cat("Number of sample size: n =", nrow(x$allres$data), "\n")
	cat("Number of categories in data:", length(unique(sort(as.matrix(x$allres$data))))," (" , unique(sort(as.matrix(x$allres$data))), ")", "\n")
	cat("The estimated linkage map is inserted in <YOUR OUTPUT NAME>$map \n")
	cat("To visualize the associated 3-dimensional map consider plot(<YOUR OUTPUT NAME>) \n")
	cat("----------------------- \n")
	cat("To visualize the other possible 3-dimensional maps consider plot(<YOUR OUTPUT NAME>$allres) \n")
	cat("To build a linkage map for your desired 3-dimensional map consider builMap() function \n")
}