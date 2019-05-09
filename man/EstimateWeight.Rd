\name{EstimateWeight}
\alias{EstimateWeight}
\title{EstimateWeight function}
\description{
  Function to estimate the weight matrix.  Based on the set of marker genes for each cell type, this function estimate the cell-specific proportions (weight) for each sample.
}
\usage{
EstimateWeight(mix_ob, gene_list, method = "LM", l = 0, u = 2^34)
}
\arguments{
	\item{data}{data matrix of the mixture signals, with genes in row, cell type in column, in anti-log scale.}
	
	\item{mix_ob}{data matrix of the mixture signals, with genes in row, cell type in column, in anti-log scale.}
	
	\item{gene_list}{list of the length in the number of tissue types. Each list element contains gene symbols representing the tissue type.}

	\item{method}{methods used in estimating true signals for each tissue type.  Default to 'LM' for linear regression.  
						Other methods included 'QP_LM' (quadradic programming with constraint on the estimated parameter on linear regression)}

	\item{l,u}{values for the lower- (l) and upper- (u) bound used in setting the vector for values of b0 (bvec for solve.QP) in solving for quadratic programming. Defaults to 0 and 2^34 respectively.}

}
\value{
A list object containing two elements:
	\item{weight}{matrix of estimated weight, with cell type in row, tissue types in column.}

	\item{mse}{means square error of the fitted linear model.  mse is 'NULL' if method is 'QP_LM'.}
}
\seealso{
	\code{\link{GSM_QP}}, 
	\code{solve.QP}
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

	weight.orig <- weight
	weight <- weight[4:14, ]

	# Generate a list object set of marker genes for each cell type
	unique_cell_type <- unique(cell.gene[,2]) 
	gene_list <- list()
	for( i in 1 : length(unique_cell_type)){
		gene_list[[i]] <- cell.gene[cell.gene[,2] == unique_cell_type[i],1]
	} 
	names(gene_list) <- unique_cell_type

	# Estimate weight matrix
	estimated_weight <- EstimateWeight(2^mix, gene_list, method="LM")
	
	t.weight <- as.numeric(weight)/100
	e.weight <- as.numeric(t(estimated_weight$weight))
	limits <- range(t.weight, e.weight)

	plot(t.weight, e.weight, xlim=limits, pch=19, ylim=limits)
	cor(t.weight, e.weight)
	abline(a=0, b=1)
}
