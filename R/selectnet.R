#-------------------------------------------------------------------------------#
# Package: Network-Based Genome-Wide Association Studies                        #
# selectnet(): Select optimal network using three methods: Bayesian information #
#criterion(BIC), extended BIC (eBIC), and Akaike information criterion (AIC)    #
# Author: Pariya Behrouzi                                                       #
# Email: <pariya.Behrouzi@gmail.com>                                            #
# Date: Nov 21th 2017                                                           #
# Version: 0.0.1-1                                                              #
#-------------------------------------------------------------------------------#

selectnet = function(netgwas.obj, opt.index = NULL, criteria = NULL, ebic.gamma = 0.5, ncores = NULL, verbose=TRUE)
{
	if(! is.null(opt.index)){
		if(is.null(ncores)) ncores <- 1
		sel <- list( )
		sel$opt.adj		<-  netgwas.obj$path[[opt.index]]
		sel$opt.theta	<- netgwas.obj$Theta[[opt.index]]
		sel$opt.rho		<- netgwas.obj$rho[opt.index]
		sel$opt.index	<- opt.index
	}else{	
		if(is.null(ncores)) ncores <- detectCores() - 1
		if(is.null(criteria)) criteria <- "ebic"
		if(!ncores) ncores = 1		
		sel	= model.selection( netgwas.obj, criterion = criteria, lower.upper=lower.upper, ebic.gamma=ebic.gamma, ncores = ncores, verbose = verbose)
	}		

	theta <- as.matrix(sel$opt.theta)
	rownames(theta) <- colnames(theta)
	par.cor <- calculate.strength.theta(theta)
	par.cor[upper.tri(par.cor)] <- 0
	par.cor <- par.cor + t(par.cor)
	diag(par.cor ) <- 1
	par.cor[par.cor > 1] <- 1
	sel$par.cor <- Matrix(par.cor)
	
	rm(par.cor, theta, netgwas.obj)
	class(sel) = "select"
	return(sel)
}

readkey <- function()
{
    cat ("Press [enter] to continue")
    line <- readline()
}

#-----------------------------------------------------#
#different plots for class "select"                   #
# Author: Pariya Behrouzi                             #
# Email: <pariya.Behrouzi@gmail.com>                  #
#-----------------------------------------------------#
plot.select = function(x, vis= NULL, xlab= NULL, ylab= NULL, n.var = NULL, vertex.label = FALSE, ..., layout = NULL, label.vertex = "all", vertex.size = NULL, vertex.color = "red" , edge.color = "gray29",  sel.label = NULL)
{
	if(class(x) != "select") stop("The input of this plot function should be from \"select\" class (More info in: selectnet( ) ). \n")
	if(is.null(vis)) vis <-  "CI"
	
	if(vis == "CI" ){

		if(! vertex.label) {
			vertex.label = NA
		}else{
			if(!is.null(colnames(x$opt.adj) )) 
				{
					vertex.label = colnames(x$opt.adj)
			}else{
					vertex.label= NA
			}
		}
		if(is.null(vertex.size)) vertex.size = 7
		
		adj = graph.adjacency(as.matrix(x$opt.adj), mode="undirected", diag=FALSE)
		if(is.null(n.var)) 
			{
				memberships = 1
				vertex.color = "red"
		}else{
		LG = length(n.var)
		memberships = NULL
		i = 1
		while( i <= LG)
			{
				grp <- rep(i, n.var[i])
				memberships = c(memberships, grp)
				i = i + 1
			}
			color <- sample(rainbow(max(memberships)+10, alpha=0.3), max(memberships))
			vertex.color = color[memberships]
		}
		adj$layout	= layout.fruchterman.reingold 
		
		plot(adj, vertex.color= vertex.color , edge.color='gray40', vertex.size = vertex.size, vertex.label = vertex.label, vertex.label.dist = 0, main= "Selected graph")	 	  
		if(length(memberships) > 1) legend("bottomright", paste("group", 1:length(n.var)), cex=0.7, col= color, pch=rep(20,10))
		readkey()
		
		if(is.null(xlab)) xlab <- ""
		if(is.null(ylab)) ylab <- ""
		image(as.matrix(x$opt.adj), xlab= xlab, ylab= ylab, col = gray.colors(256) ,main="Conditional dependence relationships" , cex=0.8) 
	}	
	
	if(vis == "interactive"){
		adj <- as.matrix(x$opt.adj)
		  
		if(is.null(vertex.size)) vertex.size <- 7
		if(is.null(label.vertex)) label.vertex <- "all"
		if(is.null(vertex.color)) vertex.color <- "red"
		if(is.null(edge.color)) edge.color <- "gray29"
		if(is.null(sel.label)) sel.label <- NULL
		if((label.vertex == "some") && (is.null(sel.label )) ) stop("Please select some vertex label(s) or fix label.vertex to either none or all.")

		p <- ncol(adj)
		A <- graph.adjacency(adj, mode= "undirected")
		if(is.null(layout)) layout <- layout_with_fr(A)
		V(A)$label.cex <- 1
	 
		if(label.vertex == "none")
		{
			V(A)$label <- NA
			tkplot(A, layout=layout, vertex.color=vertex.color, edge.color=edge.color, vertex.size=vertex.size, vertex.label.dist=0)  
		}

		if(label.vertex == "some") 
		{
			V(A)$label <- colnames(adj)
			tkplot(A, vertex.label=ifelse(V(A)$label %in% sel.label, V(A)$label, NA ), layout=layout, vertex.color=vertex.color, edge.color=edge.color, vertex.size=vertex.size, vertex.label.dist=0)
		}
		if(label.vertex == "all") 
		{	
			V(A)$label <- colnames(adj)
			tkplot(A, vertex.label=colnames(adj) , layout=layout, vertex.color=vertex.color, edge.color=edge.color, vertex.size=vertex.size, vertex.label.dist=0)  
		}
	}
}

calculate.strength.theta <- function(theta){
  p <- ncol(theta)
  cond.cor <- matrix(NA, ncol=p, nrow=p)
  for(i in 1:nrow(theta))
  {
    for(j in 1:ncol(theta))
    {
      cond.cor[i,j] <- - theta[i,j]/ sqrt(theta[i,i])* sqrt(theta[j,j])
    }
  }
  id <- colnames(theta)
  rownames(cond.cor) <- id
  colnames(cond.cor) <- id
  return(cond.cor)
}


#-----------------------------------------------------#
#   		Summary for class "select"                #
#-----------------------------------------------------#
print.select = function(x, ...){
	cat("To plot selected graph: plot(<YOUR OUTPUT NAME>) \n")
	cat("To visualize an interactive network: plot(<YOUR OUTPUT NAME>, vis= \"interactive\") \n")
	cat("To plot partial correlations: image(<YOUR OUTPUT NAME>$par.cor) \n")
}