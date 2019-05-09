\name{GSM_QP}
\alias{GSM_QP}
\title{GSM_QP function}
\description{
  Function to estimate the true signals for each tissue type with quadratic programming.
}
\usage{
GSM_QP(ob, weight, l = 0, u = 2^34, meq =0)
}
\arguments{
	\item{ob}{data matrix of the mixture signals (in anti-log scale), with genes in row, cell type in column.}

	\item{weight}{weight matrix, with cell type in column, tissue types in row.}

	\item{l,u}{values for the lower- (l) and upper- (u) bound used in setting the vector for values of b0 (bvec for solve.QP) in solving for quadratic programming. Defaults to 0 and 2^34 respectively.}
	
	\item{meq}{default to zero (used to set meq for solve.QP)}

}
\details{
	This functions depends on the solve.QP function from quadprog package.  Its major job is to set up values for the needed parameters call solve.QP to obtain the solution.
}
\value{
	\item{sol}{estimated true signals for each gene (in row) in each cell type (in column).}
}
\examples{
	## load package DSA
	library(DSA)
	
	## load sample data
	data('mix.signals')
	data('cell.gene')
	data('weight')
	
	# In this data, thr first three samples contain signal from only one cell	
	pure <- mix[, 1:3]
	mix <- mix[, 4:14]

	weight <- weight[4:14, ]

	data <- as.matrix(2^mix)
	estimate_weight <- as.matrix(weight)
	estimate_weight <- estimate_weight/rowSums(estimate_weight)

	# Obtain deconvolution
	paraM <- GSM_QP(data, estimate_weight, l = min(data), u = max(data) , meq =0) 
	
	# check the estimated cell-type specific signals
	for(i in 1:3){
		print(cor(2^(pure[,i]),paraM[,i]))
	}
}
